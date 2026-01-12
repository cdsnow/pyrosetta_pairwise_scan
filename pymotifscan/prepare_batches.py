#!/usr/bin/env python3
"""
Splits candidate PDBs into batches for stable parallel processing.
"""

import os
import argparse
from pathlib import Path
import math

def discover_candidates(candidates_dir):
    candidates_path = Path(candidates_dir)
    return sorted(candidates_path.glob('**/*.pdb'))

def main():
    parser = argparse.ArgumentParser(description="Split candidates into batches")
    parser.add_argument('candidates_dir', help="Directory with candidate PDBs")
    parser.add_argument('--num-batches', '-n', type=int, default=20, help="Number of batches")
    parser.add_argument('--output-dir', '-o', default='batches', help="Output directory for batch lists")
    args = parser.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)

    print("Discovering candidates...")
    candidates = discover_candidates(args.candidates_dir)
    print(f"Found {len(candidates)} candidates.")

    # Sort by size (LPT - largest processing time first approximation)
    print("Sorting by file size...")
    candidates.sort(key=lambda p: os.path.getsize(str(p)), reverse=True)

    # Distribute round-robin to balance load
    batches = [[] for _ in range(args.num_batches)]
    for i, candidate in enumerate(candidates):
        batches[i % args.num_batches].append(str(candidate))

    print(f"Writing {args.num_batches} batch files to {args.output_dir}...")
    
    run_script_lines = [
        "#!/bin/bash",
        "# Run this script to execute all batches in parallel.",
        "",
        "# Initialize conda",
        "__conda_setup=\"$('/home/csnow/miniconda3/bin/conda' 'shell.bash' 'hook' 2> /dev/null)\"",
        "if [ $? -eq 0 ]; then",
        "    eval \"$__conda_setup\"",
        "else",
        "    if [ -f \"/home/csnow/miniconda3/etc/profile.d/conda.sh\" ]; then",
        "        . \"/home/csnow/miniconda3/etc/profile.d/conda.sh\"",
        "    else",
        "        export PATH=\"/home/csnow/miniconda3/bin:$PATH\"",
        "    fi",
        "fi",
        "unset __conda_setup",
        "",
        "conda activate pyrosetta",
        "",
        "mkdir -p results_batches",
        ""
    ]

    for i, batch in enumerate(batches):
        batch_file = os.path.join(args.output_dir, f"batch_{i}.txt")
        with open(batch_file, 'w') as f:
            for path in batch:
                f.write(f"{path}\n")
        
        # Add command to run script
        result_file = f"results_batches/result_{i}.json"
        log_file = f"results_batches/log_{i}.txt"
        
        # Using a python script for each batch.
        # We assume runscan.sh environment setup is handled by the user or this script.
        python_path = "/home/csnow/miniconda3/envs/pyrosetta/bin/python3"
        cmd = (f"{python_path} pymotifscan/energy_scan.py "
               f"--file-list {batch_file} "
               f"--target flagorigin_and_loop.pdb "
               f"--output {result_file} "
               f"--workers 1 "  # Each batch is serial internally to avoid multiprocessing issues!
               f"> {log_file} 2>&1 &")
        
        run_script_lines.append(f"echo 'Starting batch {i}...'")
        run_script_lines.append(cmd)
        run_script_lines.append(f"PID_{i}=$!")

    run_script_lines.append("")
    run_script_lines.append("echo 'All batches started. Waiting for completion...'")
    run_script_lines.append("wait")
    run_script_lines.append("echo 'All batches complete.'")
    run_script_lines.append("")
    run_script_lines.append("# Merge results")
    run_script_lines.append(f"{python_path} pymotifscan/merge_results.py results_batches merged_results.json")

    with open("run_batches.sh", "w") as f:
        f.write("\n".join(run_script_lines))
    
    os.chmod("run_batches.sh", 0o755)
    print("Created run_batches.sh")

if __name__ == "__main__":
    main()
