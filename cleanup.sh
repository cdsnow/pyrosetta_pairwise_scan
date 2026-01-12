#!/bin/bash
# Clean up all generated files from the pipeline

rm -rf candidates_transformed
rm -rf candidates_pruned
rm -rf batches/ results_batches/ batches_prune/
rm -f merged_results.json pruned_results.json
rm -f energy_scan_results.json
rm -f *.checkpoint
rm -f run_batches.sh run_prune_batches.sh