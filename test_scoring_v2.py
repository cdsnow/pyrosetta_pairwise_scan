#!/usr/bin/env python3
"""
Improved test for scoring candidates_v2 PDB files with PyRosetta.
Investigates water and MAM loading issues.
"""

import pyrosetta
from pyrosetta import rosetta
import os
import glob

def extract_first_model(pdb_file):
    """Extract first model from a multi-model PDB file."""
    with open(pdb_file, 'r') as f:
        lines = f.readlines()

    model_lines = []
    in_model = False
    found_model = False

    for line in lines:
        if line.startswith('MODEL'):
            in_model = True
            found_model = True
            continue
        elif line.startswith('ENDMDL'):
            break
        elif line.startswith('ATOM') or line.startswith('HETATM'):
            if found_model and in_model:
                model_lines.append(line)
            elif not found_model:  # No MODEL records
                model_lines.append(line)

    return model_lines


def analyze_pdb_content(lines):
    """Analyze what residues are in the PDB."""
    residues = {}
    for line in lines:
        if line.startswith('ATOM') or line.startswith('HETATM'):
            resname = line[17:20].strip()
            chain = line[21]
            resnum = line[22:26].strip()
            key = (chain, resnum, resname)
            if key not in residues:
                residues[key] = []
            atom_name = line[12:16].strip()
            residues[key].append(atom_name)
    return residues


def main():
    # Initialize PyRosetta with params files
    params_dir = "params"
    params_files = glob.glob(os.path.join(params_dir, "*.params"))

    params_str = " ".join(params_files)

    # Key flags:
    # -ignore_unrecognized_res false: fail if residue not recognized
    # -load_PDB_components false: don't try to auto-load from CCD
    # -ignore_waters false: load water molecules
    init_flags = f"-extra_res_fa {params_str} -ignore_waters false -mute all"

    print(f"Initializing PyRosetta with flags: {init_flags}\n")
    pyrosetta.init(init_flags)

    sfxn = pyrosetta.get_fa_scorefxn()

    # Test files
    test_files = [
        ("candidates_v2/F5_ASP_SC/ARG_NE_NH1.pdb", "ACT + MGN"),
        ("candidates_v2/F3_LYS_SC/ARG_wNE.pdb", "MAM + MGN + HOH"),
        ("candidates_v2/F1_ASP_BB_N/SER_wOG.pdb", "ALA + MEO + HOH"),
        ("candidates_v2/F5_ASP_SC/ASN_ND2.pdb", "ACT + ACM"),
        ("candidates_v2/F2_TYR_SC/ARG_NE_NH1.pdb", "MEO + MGN"),
    ]

    for pdb_file, expected in test_files:
        if not os.path.exists(pdb_file):
            print(f"\nFile not found: {pdb_file}")
            continue

        print(f"\n{'='*70}")
        print(f"Testing: {pdb_file}")
        print(f"Expected: {expected}")
        print('='*70)

        try:
            # Extract first model
            model_lines = extract_first_model(pdb_file)

            # Analyze content
            residues = analyze_pdb_content(model_lines)
            print(f"\nPDB contains {len(residues)} residue(s):")
            for (chain, resnum, resname), atoms in sorted(residues.items()):
                print(f"  Chain {chain} #{resnum} {resname}: {len(atoms)} atoms - {atoms[:5]}{'...' if len(atoms)>5 else ''}")

            # Write temp file
            temp_pdb = "/tmp/test_model.pdb"
            with open(temp_pdb, 'w') as f:
                f.writelines(model_lines)
                f.write("END\n")

            # Load the pose
            pose = pyrosetta.pose_from_pdb(temp_pdb)

            print(f"\nLoaded pose with {pose.total_residue()} residue(s):")
            for i in range(1, pose.total_residue() + 1):
                res = pose.residue(i)
                print(f"  {i}. {res.name3()} chain={pose.pdb_info().chain(i)} "
                      f"pdb#{pose.pdb_info().number(i)} atoms={res.natoms()}")

            # Check what's missing
            loaded_keys = set()
            for i in range(1, pose.total_residue() + 1):
                chain = pose.pdb_info().chain(i)
                resnum = str(pose.pdb_info().number(i))
                resname = pose.residue(i).name3()
                loaded_keys.add((chain, resnum, resname))

            pdb_keys = set(residues.keys())
            missing = pdb_keys - loaded_keys
            if missing:
                print(f"\n  WARNING: Missing residues: {missing}")

            # Score
            total_score = sfxn(pose)
            print(f"\nTotal score: {total_score:.3f}")

            # Breakdown
            print("Per-residue scores:")
            for i in range(1, pose.total_residue() + 1):
                res_score = pose.energies().residue_total_energy(i)
                print(f"  {pose.residue(i).name3()}: {res_score:.3f}")

        except Exception as e:
            print(f"ERROR: {e}")
            import traceback
            traceback.print_exc()

    print("\n" + "="*70)
    print("Testing complete!")


if __name__ == '__main__':
    main()
