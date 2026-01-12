#!/usr/bin/env python3
"""
Merges JSON results from multiple batch files.
"""

import json
import os
import sys
import glob

def merge_results(batch_dir, output_file):
    merged_results = {}
    metadata = {}
    
    files = glob.glob(os.path.join(batch_dir, "result_*.json"))
    print(f"Found {len(files)} result files in {batch_dir}")

    for fpath in files:
        try:
            with open(fpath, 'r') as f:
                data = json.load(f)
                
                # Merge results
                if 'results' in data:
                    merged_results.update(data['results'])
                
                # Keep last metadata found (should be mostly same)
                if 'metadata' in data:
                    metadata = data['metadata']
        except Exception as e:
            print(f"Error reading {fpath}: {e}")

    final_output = {
        'metadata': metadata,
        'results': merged_results
    }
    
    # Update total count
    final_output['metadata']['total_processed'] = len(merged_results)
    
    with open(output_file, 'w') as f:
        json.dump(final_output, f, indent=2)
    
    print(f"Merged results saved to {output_file}")
    print(f"Total candidates processed: {len(merged_results)}")

if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage: python3 merge_results.py <batch_dir> <output_file>")
        sys.exit(1)
    
    merge_results(sys.argv[1], sys.argv[2])
