#!/bin/bash

# initialize conda
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

# activate environment
conda activate pyrosetta

# run scan
python3 pairwise_energy_scan.py \
    --candidates candidates_v2 \
    --target flagorigin_and_loop.pdb \
    --output energy_scan_results_v2.json