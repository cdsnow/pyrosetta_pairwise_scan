
import sys
import os
from pyrosetta import init, Pose, get_fa_scorefxn
from pyrosetta.rosetta.core.import_pose import pose_from_pdbstring
from pyrosetta.rosetta.core.scoring import ScoreType

def get_energies(pose, sfxn):
    sfxn(pose)
    energies = pose.energies()
    return {
        'total': energies.total_energy(),
        'fa_rep': energies.total_energies()[ScoreType.fa_rep],
    }

def test_scoring():
    init("-mute all")
    sfxn = get_fa_scorefxn()
    
    target_file = "flagorigin_and_loop.pdb"
    candidate_file = "candidates_v2/F6_ASP_SC/ARG_N_NE_NH2.pdb"
    
    with open(target_file) as f:
        target_str = f.read()
        
    with open(candidate_file) as f:
        cand_lines = f.readlines()
        
    # Extract Model 1 only
    model_1 = []
    in_model = False
    for line in cand_lines:
        if line.startswith("MODEL"):
            in_model = True
        elif line.startswith("ENDMDL"):
            if in_model: break
        elif in_model or not any(l.startswith("MODEL") for l in cand_lines):
            # Handle both multi-model and single-model cases broadly
            if line.startswith("ATOM") or line.startswith("HETATM"):
                model_1.append(line)
    
    cand_str = "".join(model_1)
    
    # CASE 1: Original (Collision)
    combined = target_str + cand_str
    pose = Pose()
    try:
        pose_from_pdbstring(pose, combined)
        e = get_energies(pose, sfxn)
        print(f"Original Score (Combined): Total={e['total']:.2f}, fa_rep={e['fa_rep']:.2f}")
    except Exception as ex:
        print(f"Original Failed: {ex}")

    # CASE 2: Filtered (Keep only Chain Z from candidate)
    filtered_cand_lines = [l for l in model_1 if l[21] == 'Z']
    filtered_cand_str = "".join(filtered_cand_lines)
    
    combined_filt = target_str + filtered_cand_str
    pose_filt = Pose()
    try:
        pose_from_pdbstring(pose_filt, combined_filt)
        e = get_energies(pose_filt, sfxn)
        print(f"Filtered Score (Only Z):   Total={e['total']:.2f}, fa_rep={e['fa_rep']:.2f}")
    except Exception as ex:
        print(f"Filtered Failed: {ex}")

if __name__ == "__main__":
    test_scoring()
