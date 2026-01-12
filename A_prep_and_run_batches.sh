# Step 1: Transform partial sidechains to small molecules
python pymotifscan/transform_candidates.py --indir candidates --outdir candidates_transformed

# Step 2: Prepare and run energy scan batches
python pymotifscan/prepare_batches.py candidates_transformed --num-batches 20
./run_batches.sh