#!/usr/bin/env python3
"""
Parallelized PyRosetta energy scanner for candidate-target binding evaluation.

Computes delta_E = E(complex) - E(target) - E(motif) for each candidate model.
Supports graceful restart from interrupted scans via checkpoint files.
Optimized to initialize PyRosetta once per worker process.
"""

import os
import sys
import json
import glob

# Debug info
# print(f"DEBUG: sys.version={sys.version}", file=sys.stderr)
# print(f"DEBUG: sys.executable={sys.executable}", file=sys.stderr)

import argparse
import tempfile
from pathlib import Path
import multiprocessing
from datetime import datetime
import signal
from contextlib import contextmanager

# Global variables for worker processes
_worker_sfxn = None
_worker_target_pdb_string = None
_worker_target_E = None


class TimeoutException(Exception):
    pass


@contextmanager
def silence_stdout_stderr():
    """Context manager to suppress stdout and stderr."""
    try:
        # Save original file descriptors
        saved_stdout_fd = os.dup(sys.stdout.fileno())
        saved_stderr_fd = os.dup(sys.stderr.fileno())

        # Open /dev/null
        with open(os.devnull, 'w') as devnull:
            # Redirect stdout/stderr to /dev/null
            os.dup2(devnull.fileno(), sys.stdout.fileno())
            os.dup2(devnull.fileno(), sys.stderr.fileno())

        try:
            yield
        finally:
            # Restore original file descriptors
            os.dup2(saved_stdout_fd, sys.stdout.fileno())
            os.dup2(saved_stderr_fd, sys.stderr.fileno())
            
            # Close duplicated fds
            os.close(saved_stdout_fd)
            os.close(saved_stderr_fd)
    except Exception:
        # Fallback if redirection fails (e.g. if stdout/stderr are not real files)
        yield


@contextmanager
def time_limit(seconds):
    def signal_handler(signum, frame):
        raise TimeoutException("Timed out!")
    signal.signal(signal.SIGALRM, signal_handler)
    signal.alarm(seconds)
    try:
        yield
    finally:
        signal.alarm(0)


def split_multi_model_pdb(pdb_path):
    """Split a multi-model PDB file into individual model strings.

    Returns list of (model_num, pdb_string, metadata) tuples.
    """
    models = []
    current_model = []
    model_num = None
    metadata = {}

    with open(pdb_path, 'r') as f:
        for line in f:
            if line.startswith('MODEL'):
                model_num = int(line.split()[1])
                current_model = []
                metadata = {}
            elif line.startswith('ENDMDL'):
                if current_model and model_num is not None:
                    models.append((model_num, ''.join(current_model), metadata))
            elif line.startswith('REMARK') and 'Score:' in line:
                parts = line.split()
                for i, p in enumerate(parts):
                    if p == 'Score:' and i + 1 < len(parts):
                        try:
                            metadata['original_score'] = float(parts[i + 1])
                        except ValueError:
                            pass
                    if p == 'ID:' and i + 1 < len(parts):
                        metadata['original_id'] = parts[i + 1]
            elif line.startswith('ATOM') or line.startswith('HETATM') or line.startswith('TER'):
                current_model.append(line)

    return models


def get_energies(pose, sfxn):
    """Extract key energy terms from a scored pose."""
    from pyrosetta.rosetta.core.scoring import ScoreType

    sfxn(pose)
    energies = pose.energies()

    return {
        'total': energies.total_energy(),
        'fa_rep': energies.total_energies()[ScoreType.fa_rep],
        'hbond': sum([
            energies.total_energies()[ScoreType.hbond_sc],
            energies.total_energies()[ScoreType.hbond_bb_sc],
            energies.total_energies()[ScoreType.hbond_sr_bb],
            energies.total_energies()[ScoreType.hbond_lr_bb],
        ])
    }


def extract_chain(pose, chain_id):
    """Extract a single chain from a pose, returning a new pose."""
    import pyrosetta
    from pyrosetta.rosetta.core.pose import Pose

    new_pose = Pose()
    pdb_info = pose.pdb_info()

    for i in range(1, pose.total_residue() + 1):
        if pdb_info.chain(i) == chain_id:
            res = pose.residue(i).clone()
            if new_pose.total_residue() == 0:
                new_pose.append_residue_by_jump(res, 1)
            else:
                # If it's a polymer and following the previous residue, try to bond
                # Otherwise, jump. For motifs, jump is generally safer for ligands.
                if res.is_polymer() and new_pose.residue(new_pose.total_residue()).is_polymer():
                    new_pose.append_residue_by_bond(res)
                else:
                    new_pose.append_residue_by_jump(res, new_pose.total_residue())

    return new_pose


def evaluate_single_model(model_pdb_string, target_pdb_string, sfxn, target_E_cache=None):
    """Evaluate a single model's binding energy.

    Args:
        model_pdb_string: PDB content as string for this model
        target_pdb_string: PDB content as string for the target
        sfxn: PyRosetta score function
        target_E_cache: Pre-computed target energies (optional, for efficiency)

    Returns:
        Dictionary with energy values, or None if model couldn't be loaded
    """
    import pyrosetta
    from pyrosetta.rosetta.core.import_pose import pose_from_pdbstring

    # Combine target and motif PDB strings
    combined_pdb_string = target_pdb_string + model_pdb_string

    # Load the complex (model contains both target fragment and motif)
    complex_pose = pyrosetta.Pose()
    pose_from_pdbstring(complex_pose, combined_pdb_string)

    # Check what loaded
    if complex_pose.total_residue() == 0:
        return {'skip_reason': 'no_residues_loaded'}

    # Check chains present
    pdb_info = complex_pose.pdb_info()
    chains = set()
    for i in range(1, complex_pose.total_residue() + 1):
        chains.add(pdb_info.chain(i))

    has_target_chain = 'F' in chains
    has_motif_chain = 'Z' in chains

    if not has_motif_chain:
        return {'skip_reason': 'no_motif_chain_Z'}

    # Get target energy (use cache if provided)
    if target_E_cache is not None:
        target_E = target_E_cache
    else:
        # Create pose from string
        target_pose = pyrosetta.Pose()
        pose_from_pdbstring(target_pose, target_pdb_string)
        target_E = get_energies(target_pose, sfxn)

    # Score the complex
    complex_E = get_energies(complex_pose, sfxn)

    # Try to extract and score motif separately
    motif_pose = extract_chain(complex_pose, 'Z')
    if motif_pose.total_residue() > 0:
        motif_E = get_energies(motif_pose, sfxn)
    else:
        # Can't extract motif, use zero as approximation
        motif_E = {'total': 0.0, 'fa_rep': 0.0, 'hbond': 0.0}

    # Compute deltas: delta = complex - target - motif
    delta = {f'delta_{k}': complex_E[k] - target_E[k] - motif_E[k] for k in complex_E}

    return {
        'complex': complex_E,
        'target': target_E,
        'motif': motif_E,
        'chains_found': list(chains),
        'n_residues': complex_pose.total_residue(),
        **delta
    }


def worker_init(target_path, params_files):
    """Initialize the worker process.

    Sets up PyRosetta, loads params, and pre-computes target energy.
    """
    global _worker_sfxn, _worker_target_pdb_string, _worker_target_E

    import pyrosetta
    from pyrosetta.rosetta.core.import_pose import pose_from_pdbstring
    import time
    import random

    # Stagger initialization to avoid potential race conditions
    time.sleep(random.random() * 0.5)

    # Build init flags
    params_str = " ".join(params_files)
    init_flags = f"-extra_res_fa {params_str} -ignore_waters false -mute all"

    # Initialize PyRosetta
    # Note: With spawn context, this is always a fresh process, so we just init.
    # If init fails, we want it to raise an exception.
    try:
        with silence_stdout_stderr():
            pyrosetta.init(init_flags, set_logging_handler=None)
    except Exception as e:
        # If it says already initialized, that's fine (though unlikely with spawn)
        # But if it's another error, we should know.
        print(f"[{os.getpid()}] Warning: PyRosetta init raised: {e}")

    # Create score function
    _worker_sfxn = pyrosetta.get_fa_scorefxn()

    # Read target PDB content
    with open(target_path, 'r') as f:
        _worker_target_pdb_string = f.read()

    # Pre-compute target energy
    target_pose = pyrosetta.Pose()
    pose_from_pdbstring(target_pose, _worker_target_pdb_string)
    _worker_target_E = get_energies(target_pose, _worker_sfxn)


def process_candidate_pdb(pdb_path):
    """Worker function to process a single candidate PDB file.

    Uses globally initialized PyRosetta and target data.
    """
    global _worker_sfxn, _worker_target_pdb_string, _worker_target_E

    results = {
        'pdb_file': str(pdb_path),
        'models': [],
        'skipped': 0,
        'error': None
    }

    try:
        # Set a 5-minute timeout for processing the file
        with time_limit(300):
            models = split_multi_model_pdb(pdb_path)
            results['total_models'] = len(models)

            for model_num, pdb_string, metadata in models:
                try:
                    energies = evaluate_single_model(
                        pdb_string,
                        _worker_target_pdb_string,
                        _worker_sfxn,
                        _worker_target_E
                    )
                    if energies:
                        # Check if model was skipped
                        if 'skip_reason' in energies:
                            results['skipped'] += 1
                            continue

                        model_result = {
                            'model_num': model_num,
                            **metadata,
                            **energies
                        }
                        results['models'].append(model_result)
                except Exception as e:
                    results['models'].append({
                        'model_num': model_num,
                        'error': str(e)
                    })

    except TimeoutException:
        results['error'] = "Processing timed out (300s)"
    except Exception as e:
        results['error'] = str(e)

    return results


def discover_candidates(candidates_dir):
    """Find all candidate PDB files."""
    candidates_path = Path(candidates_dir)
    pdb_files = sorted(candidates_path.glob('**/*.pdb'))
    return pdb_files


def load_checkpoint(checkpoint_path):
    """Load checkpoint data if it exists."""
    if os.path.exists(checkpoint_path):
        with open(checkpoint_path, 'r') as f:
            return json.load(f)
    return {'completed': {}, 'started': None}


def save_checkpoint(checkpoint_path, data):
    """Save checkpoint data atomically."""
    temp_path = checkpoint_path + '.tmp'
    with open(temp_path, 'w') as f:
        json.dump(data, f)
    os.replace(temp_path, checkpoint_path)


def run_energy_scan(candidates_dir, target_path, output_path, checkpoint_path=None,
                    num_workers=None, batch_size=10, file_list=None):
    """Run the parallelized energy scan with checkpoint support.

    Args:
        candidates_dir: Path to directory containing candidate PDB files
        target_path: Path to the target PDB file
        output_path: Path to save final results JSON
        checkpoint_path: Path for checkpoint file (auto-generated if None)
        num_workers: Number of parallel workers (default: CPU count - 1)
        batch_size: Number of candidates to process before checkpointing
        file_list: Optional list of PDB files to process (overrides directory scan)

    Returns:
        Dictionary with all results
    """
    if num_workers is None:
        # Use half of available CPUs to reduce memory pressure and risk of OOM
        num_workers = max(1, multiprocessing.cpu_count() // 2)

    if checkpoint_path is None:
        checkpoint_path = output_path + '.checkpoint'

    # Discover all candidates
    if file_list:
        with open(file_list, 'r') as f:
            # Assume one path per line
            all_candidates = [Path(line.strip()) for line in f if line.strip()]
        # Resolve relative paths against current dir if needed, but assuming user provides correct paths
        print(f"Loaded {len(all_candidates)} candidates from {file_list}")
    else:
        all_candidates = discover_candidates(candidates_dir)
        print(f"Found {len(all_candidates)} candidate PDB files")

    # Load checkpoint
    checkpoint = load_checkpoint(checkpoint_path)
    completed = checkpoint['completed']

    # Filter out already completed candidates
    pending = [c for c in all_candidates if str(c) not in completed]
    print(f"Pending: {len(pending)}, Already completed: {len(completed)}")

    if not pending:
        print("All candidates already processed!")
        return completed

    # Sort pending candidates by file size (largest first) to optimize parallel utilization
    # This prevents the "straggler problem" where a large file at the end blocks completion.
    pending.sort(key=lambda p: os.path.getsize(str(p)), reverse=True)
    print("Sorted candidates by size (descending) for optimal scheduling.")

    # Prepare arguments for workers (just paths, rest is in init)
    work_items = [str(pdb) for pdb in pending]

    # Find params files for initialization
    params_dir = "params"
    # Ensure params_dir path is absolute or correct relative to execution
    if not os.path.isabs(params_dir):
         params_dir = os.path.abspath(params_dir)
    
    params_files = glob.glob(os.path.join(params_dir, "*.params"))
    
    # Process in batches with checkpointing
    if checkpoint.get('started') is None:
        checkpoint['started'] = datetime.now().isoformat()

    print("STARK_MARKER_99999")
    print(f"Starting scan (num_workers={num_workers})...")
    
    # We need to make target_path absolute just in case workers change dir (unlikely but safe)
    abs_target_path = os.path.abspath(target_path)

    processed = 0
    
    if num_workers > 1:
        # Use spawn to ensure clean PyRosetta state in workers
        ctx = multiprocessing.get_context('spawn')
        # maxtasksperchild=1 ensures a fresh process for every single file
        # This is slower but prevents memory accumulation and isolates crashes better.
        with ctx.Pool(processes=num_workers, initializer=worker_init, initargs=(abs_target_path, params_files), maxtasksperchild=1) as pool:
            try:
                # Use imap_unordered for better progress tracking
                for result in pool.imap_unordered(process_candidate_pdb, work_items):
                    pdb_file = result['pdb_file']
                    completed[pdb_file] = result
                    processed += 1

                    # Progress update
                    total_done = len(completed)
                    total_all = len(all_candidates)
                    n_models = len(result['models'])
                    n_skipped = result.get('skipped', 0)
                    n_total = result.get('total_models', n_models + n_skipped)
                    
                    # Extract parent dir for clearer logging
                    parent_dir = os.path.basename(os.path.dirname(pdb_file))
                    filename = os.path.basename(pdb_file)
                    
                    print(f"[{total_done}/{total_all}] {parent_dir}/{filename}: "
                          f"{n_models}/{n_total} models evaluated, {n_skipped} skipped")

                    # Checkpoint periodically
                    if processed % batch_size == 0:
                        checkpoint['completed'] = completed
                        save_checkpoint(checkpoint_path, checkpoint)
                        print(f"  Checkpoint saved ({total_done} candidates)")
            except Exception as e:
                print(f"\nCRITICAL ERROR IN POOL: {e}")
                print("Saving checkpoint and exiting. Some workers may have crashed.")
                checkpoint['completed'] = completed
                save_checkpoint(checkpoint_path, checkpoint)
                raise e
    else:
        # Serial execution for num_workers=1
        # This avoids multiprocessing overhead and environment issues.
        # Initialize locally
        worker_init(abs_target_path, params_files)
        
        for pdb_path_str in work_items:
            result = process_candidate_pdb(pdb_path_str)
            pdb_file = result['pdb_file']
            completed[pdb_file] = result
            processed += 1

            # Progress update
            total_done = len(completed)
            total_all = len(all_candidates)
            n_models = len(result['models'])
            n_skipped = result.get('skipped', 0)
            n_total = result.get('total_models', n_models + n_skipped)
            
            parent_dir = os.path.basename(os.path.dirname(pdb_file))
            filename = os.path.basename(pdb_file)
            
            print(f"[{total_done}/{total_all}] {parent_dir}/{filename}: "
                  f"{n_models}/{n_total} models evaluated, {n_skipped} skipped")

            if processed % batch_size == 0:
                checkpoint['completed'] = completed
                save_checkpoint(checkpoint_path, checkpoint)
                print(f"  Checkpoint saved ({total_done} candidates)")

    # Final save
    checkpoint['completed'] = completed
    checkpoint['finished'] = datetime.now().isoformat()
    save_checkpoint(checkpoint_path, checkpoint)

    # Save final results
    final_results = {
        'metadata': {
            'target': target_path,
            'candidates_dir': str(candidates_dir),
            'total_candidates': len(all_candidates),
            'started': checkpoint.get('started'),
            'finished': checkpoint['finished']
        },
        'results': completed
    }

    with open(output_path, 'w') as f:
        json.dump(final_results, f, indent=2)

    print(f"\nScan complete! Results saved to {output_path}")

    # Clean up checkpoint file on successful completion
    if os.path.exists(checkpoint_path):
        os.remove(checkpoint_path)
        print(f"Checkpoint file removed: {checkpoint_path}")

    return completed


def summarize_results(results_path):
    """Print a summary of scan results."""
    with open(results_path, 'r') as f:
        data = json.load(f)

    results = data['results']

    total_models = 0
    total_skipped = 0
    favorable_models = 0
    clashing_models = 0

    best_binding = None
    best_binding_info = None

    for pdb_file, result in results.items():
        total_skipped += result.get('skipped', 0)

        for model in result.get('models', []):
            if 'error' in model:
                continue

            total_models += 1
            delta_total = model.get('delta_total', 0)
            delta_fa_rep = model.get('delta_fa_rep', 0)

            if delta_total < 0:
                favorable_models += 1
            if delta_fa_rep > 10:
                clashing_models += 1

            if best_binding is None or delta_total < best_binding:
                best_binding = delta_total
                best_binding_info = {
                    'pdb': os.path.basename(pdb_file),
                    'model': model.get('model_num'),
                    'delta_total': delta_total,
                    'delta_hbond': model.get('delta_hbond'),
                    'delta_fa_rep': delta_fa_rep
                }

    print(f"\n{'='*60}")
    print("SCAN SUMMARY")
    print(f"{ '='*60}")
    print(f"Total PDB files: {len(results)}")
    print(f"Total models evaluated: {total_models}")
    print(f"Models skipped (incomplete fragments): {total_skipped}")
    print(f"Favorable binding (delta_total < 0): {favorable_models}")
    print(f"Clashing (delta_fa_rep > 10): {clashing_models}")

    if best_binding_info:
        print(f"\nBest binding candidate:")
        print(f"  PDB: {best_binding_info['pdb']}")
        print(f"  Model: {best_binding_info['model']}")
        print(f"  delta_total: {best_binding_info['delta_total']:.2f}")
        print(f"  delta_hbond: {best_binding_info['delta_hbond']:.2f}")
        print(f"  delta_fa_rep: {best_binding_info['delta_fa_rep']:.2f}")


def main():
    print(f"DEBUG: sys.version={sys.version}")
    print(f"DEBUG: sys.executable={sys.executable}")
    print(f"DEBUG: sys.path={sys.path}")
    
    parser = argparse.ArgumentParser(
        description="Parallelized PyRosetta energy scanner for candidate binding evaluation"
    )
    parser.add_argument(
        '--candidates', '-c',
        default='candidates',
        help='Directory containing candidate PDB files (default: candidates)'
    )
    parser.add_argument(
        '--target', '-t',
        default='flagorigin_and_loop.pdb',
        help='Target PDB file (default: flagorigin_and_loop.pdb)'
    )
    parser.add_argument(
        '--output', '-o',
        default='energy_scan_results.json',
        help='Output JSON file (default: energy_scan_results.json)'
    )
    parser.add_argument(
        '--workers', '-w',
        type=int,
        default=None,
        help='Number of parallel workers (default: CPU count - 1)'
    )
    parser.add_argument(
        '--batch-size', '-b',
        type=int,
        default=10,
        help='Checkpoint every N candidates (default: 10)'
    )
    parser.add_argument(
        '--file-list', '-f',
        default=None,
        help='File containing list of PDB paths to process (overrides candidates dir)'
    )
    parser.add_argument(
        '--summarize', '-s',
        action='store_true',
        help='Only summarize existing results file'
    )

    args = parser.parse_args()

    if args.summarize:
        if os.path.exists(args.output):
            summarize_results(args.output)
        else:
            print(f"Results file not found: {args.output}")
            sys.exit(1)
        return

    # Validate inputs
    if not os.path.exists(args.target):
        print(f"Error: Target file not found: {args.target}")
        sys.exit(1)

    if not os.path.isdir(args.candidates):
        print(f"Error: Candidates directory not found: {args.candidates}")
        sys.exit(1)

    # Run the scan
    run_energy_scan(
        candidates_dir=args.candidates,
        target_path=args.target,
        output_path=args.output,
        num_workers=args.workers,
        batch_size=args.batch_size,
        file_list=args.file_list
    )

    # Print summary
    if os.path.exists(args.output):
        summarize_results(args.output)


if __name__ == '__main__':
    # Ensure workers use the same python executable as the parent
    try:
        multiprocessing.set_executable(sys.executable)
    except Exception:
        pass
    main()