#!/bin/bash
echo "Restarting remaining batches..."
python prune_candidates.py --results batches_prune/batch_1.json --output-dir candidates_v3 --verbose > batches_prune/log_1.txt 2>&1 &
python prune_candidates.py --results batches_prune/batch_10.json --output-dir candidates_v3 --verbose > batches_prune/log_10.txt 2>&1 &
python prune_candidates.py --results batches_prune/batch_11.json --output-dir candidates_v3 --verbose > batches_prune/log_11.txt 2>&1 &
python prune_candidates.py --results batches_prune/batch_18.json --output-dir candidates_v3 --verbose > batches_prune/log_18.txt 2>&1 &
wait
echo "Done."
