#!/bin/bash
# Run this script to execute all batches in parallel.

# Initialize conda
__conda_setup="$('/home/csnow/miniconda3/bin/conda' 'shell.bash' 'hook' 2> /dev/null)"
if [ $? -eq 0 ]; then
    eval "$__conda_setup"
else
    if [ -f "/home/csnow/miniconda3/etc/profile.d/conda.sh" ]; then
        . "/home/csnow/miniconda3/etc/profile.d/conda.sh"
    else
        export PATH="/home/csnow/miniconda3/bin:$PATH"
    fi
fi
unset __conda_setup

conda activate pyrosetta

mkdir -p results_batches

echo 'Starting batch 0...'
/home/csnow/miniconda3/envs/pyrosetta/bin/python3 energy_scan.py --file-list batches/batch_0.txt --target flagorigin_and_loop.pdb --output results_batches/result_0.json --workers 1 > results_batches/log_0.txt 2>&1 &
PID_0=$!
echo 'Starting batch 1...'
/home/csnow/miniconda3/envs/pyrosetta/bin/python3 energy_scan.py --file-list batches/batch_1.txt --target flagorigin_and_loop.pdb --output results_batches/result_1.json --workers 1 > results_batches/log_1.txt 2>&1 &
PID_1=$!
echo 'Starting batch 2...'
/home/csnow/miniconda3/envs/pyrosetta/bin/python3 energy_scan.py --file-list batches/batch_2.txt --target flagorigin_and_loop.pdb --output results_batches/result_2.json --workers 1 > results_batches/log_2.txt 2>&1 &
PID_2=$!
echo 'Starting batch 3...'
/home/csnow/miniconda3/envs/pyrosetta/bin/python3 energy_scan.py --file-list batches/batch_3.txt --target flagorigin_and_loop.pdb --output results_batches/result_3.json --workers 1 > results_batches/log_3.txt 2>&1 &
PID_3=$!
echo 'Starting batch 4...'
/home/csnow/miniconda3/envs/pyrosetta/bin/python3 energy_scan.py --file-list batches/batch_4.txt --target flagorigin_and_loop.pdb --output results_batches/result_4.json --workers 1 > results_batches/log_4.txt 2>&1 &
PID_4=$!
echo 'Starting batch 5...'
/home/csnow/miniconda3/envs/pyrosetta/bin/python3 energy_scan.py --file-list batches/batch_5.txt --target flagorigin_and_loop.pdb --output results_batches/result_5.json --workers 1 > results_batches/log_5.txt 2>&1 &
PID_5=$!
echo 'Starting batch 6...'
/home/csnow/miniconda3/envs/pyrosetta/bin/python3 energy_scan.py --file-list batches/batch_6.txt --target flagorigin_and_loop.pdb --output results_batches/result_6.json --workers 1 > results_batches/log_6.txt 2>&1 &
PID_6=$!
echo 'Starting batch 7...'
/home/csnow/miniconda3/envs/pyrosetta/bin/python3 energy_scan.py --file-list batches/batch_7.txt --target flagorigin_and_loop.pdb --output results_batches/result_7.json --workers 1 > results_batches/log_7.txt 2>&1 &
PID_7=$!
echo 'Starting batch 8...'
/home/csnow/miniconda3/envs/pyrosetta/bin/python3 energy_scan.py --file-list batches/batch_8.txt --target flagorigin_and_loop.pdb --output results_batches/result_8.json --workers 1 > results_batches/log_8.txt 2>&1 &
PID_8=$!
echo 'Starting batch 9...'
/home/csnow/miniconda3/envs/pyrosetta/bin/python3 energy_scan.py --file-list batches/batch_9.txt --target flagorigin_and_loop.pdb --output results_batches/result_9.json --workers 1 > results_batches/log_9.txt 2>&1 &
PID_9=$!
echo 'Starting batch 10...'
/home/csnow/miniconda3/envs/pyrosetta/bin/python3 energy_scan.py --file-list batches/batch_10.txt --target flagorigin_and_loop.pdb --output results_batches/result_10.json --workers 1 > results_batches/log_10.txt 2>&1 &
PID_10=$!
echo 'Starting batch 11...'
/home/csnow/miniconda3/envs/pyrosetta/bin/python3 energy_scan.py --file-list batches/batch_11.txt --target flagorigin_and_loop.pdb --output results_batches/result_11.json --workers 1 > results_batches/log_11.txt 2>&1 &
PID_11=$!
echo 'Starting batch 12...'
/home/csnow/miniconda3/envs/pyrosetta/bin/python3 energy_scan.py --file-list batches/batch_12.txt --target flagorigin_and_loop.pdb --output results_batches/result_12.json --workers 1 > results_batches/log_12.txt 2>&1 &
PID_12=$!
echo 'Starting batch 13...'
/home/csnow/miniconda3/envs/pyrosetta/bin/python3 energy_scan.py --file-list batches/batch_13.txt --target flagorigin_and_loop.pdb --output results_batches/result_13.json --workers 1 > results_batches/log_13.txt 2>&1 &
PID_13=$!
echo 'Starting batch 14...'
/home/csnow/miniconda3/envs/pyrosetta/bin/python3 energy_scan.py --file-list batches/batch_14.txt --target flagorigin_and_loop.pdb --output results_batches/result_14.json --workers 1 > results_batches/log_14.txt 2>&1 &
PID_14=$!
echo 'Starting batch 15...'
/home/csnow/miniconda3/envs/pyrosetta/bin/python3 energy_scan.py --file-list batches/batch_15.txt --target flagorigin_and_loop.pdb --output results_batches/result_15.json --workers 1 > results_batches/log_15.txt 2>&1 &
PID_15=$!
echo 'Starting batch 16...'
/home/csnow/miniconda3/envs/pyrosetta/bin/python3 energy_scan.py --file-list batches/batch_16.txt --target flagorigin_and_loop.pdb --output results_batches/result_16.json --workers 1 > results_batches/log_16.txt 2>&1 &
PID_16=$!
echo 'Starting batch 17...'
/home/csnow/miniconda3/envs/pyrosetta/bin/python3 energy_scan.py --file-list batches/batch_17.txt --target flagorigin_and_loop.pdb --output results_batches/result_17.json --workers 1 > results_batches/log_17.txt 2>&1 &
PID_17=$!
echo 'Starting batch 18...'
/home/csnow/miniconda3/envs/pyrosetta/bin/python3 energy_scan.py --file-list batches/batch_18.txt --target flagorigin_and_loop.pdb --output results_batches/result_18.json --workers 1 > results_batches/log_18.txt 2>&1 &
PID_18=$!
echo 'Starting batch 19...'
/home/csnow/miniconda3/envs/pyrosetta/bin/python3 energy_scan.py --file-list batches/batch_19.txt --target flagorigin_and_loop.pdb --output results_batches/result_19.json --workers 1 > results_batches/log_19.txt 2>&1 &
PID_19=$!

echo 'All batches started. Waiting for completion...'
wait
echo 'All batches complete.'

# Merge results
/home/csnow/miniconda3/envs/pyrosetta/bin/python3 merge_results.py results_batches merged_results.json
