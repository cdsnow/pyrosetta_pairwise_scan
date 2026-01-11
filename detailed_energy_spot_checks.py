#!/usr/bin/env python3
"""
Detailed energy spot-check script for validating H-bond analysis.

Selects diverse candidates, performs detailed PyRosetta H-bond analysis,
and generates a PyMOL visualization showing interfacial interactions.
"""

import os
import sys
import json
import argparse
from pathlib import Path

# Global PyRosetta objects (initialized once)
_pyrosetta_initialized = False
_sfxn = None


def init_pyrosetta(params_dir="params"):
    """Initialize PyRosetta with ligand params."""
    global _pyrosetta_initialized, _sfxn

    if _pyrosetta_initialized:
        return

    import pyrosetta
    import glob

    params_files = glob.glob(os.path.join(params_dir, "*.params"))
    params_str = " ".join(params_files)
    init_flags = f"-extra_res_fa {params_str} -ignore_waters false -mute all"

    pyrosetta.init(init_flags, set_logging_handler=None)
    _sfxn = pyrosetta.get_fa_scorefxn()
    _pyrosetta_initialized = True


def load_results(results_path):
    """Load energy scan results from JSON."""
    with open(results_path, 'r') as f:
        return json.load(f)


def select_diverse_candidates(results, num_candidates=10, max_fa_rep=50.0, min_hbond=-1.0):
    """
    Select diverse candidates for spot-checking.

    Filters by:
    - delta_fa_rep < max_fa_rep (skip severe clashes)
    - delta_hbond < min_hbond (must have some H-bonding)

    Then selects candidates spanning the range of delta_hbond values.
    """
    candidates = []

    for pdb_path, entry in results.get('results', {}).items():
        for model in entry.get('models', []):
            if 'error' in model:
                continue
            if 'delta_fa_rep' not in model or 'delta_hbond' not in model:
                continue

            delta_fa_rep = model['delta_fa_rep']
            delta_hbond = model['delta_hbond']
            delta_total = model.get('delta_total', 0)

            # Filter criteria
            if delta_fa_rep > max_fa_rep:
                continue
            if delta_hbond > min_hbond:
                continue

            candidates.append({
                'pdb_path': pdb_path,
                'model_num': model['model_num'],
                'delta_total': delta_total,
                'delta_fa_rep': delta_fa_rep,
                'delta_hbond': delta_hbond,
                'original_id': model.get('original_id', 'unknown')
            })

    if not candidates:
        print(f"No candidates found with fa_rep < {max_fa_rep} and hbond < {min_hbond}")
        return []

    # Sort by delta_hbond (most negative = strongest H-bonding)
    candidates.sort(key=lambda x: x['delta_hbond'])

    # Select diverse candidates spanning the range
    if len(candidates) <= num_candidates:
        return candidates

    # Take evenly spaced candidates
    step = len(candidates) / num_candidates
    selected = []
    for i in range(num_candidates):
        idx = int(i * step)
        selected.append(candidates[idx])

    return selected


def split_multi_model_pdb(pdb_path):
    """Split a multi-model PDB file into individual model strings."""
    models = {}
    current_model = []
    model_num = None

    with open(pdb_path, 'r') as f:
        for line in f:
            if line.startswith('MODEL'):
                model_num = int(line.split()[1])
                current_model = []
            elif line.startswith('ENDMDL'):
                if current_model and model_num is not None:
                    models[model_num] = ''.join(current_model)
            elif line.startswith('ATOM') or line.startswith('HETATM') or line.startswith('TER'):
                current_model.append(line)

    return models


def get_interfacial_hbonds(pose, target_chain='F', motif_chain='Z'):
    """
    Extract interfacial H-bonds between target and motif chains.

    Returns list of dicts with H-bond details including atom coordinates.
    """
    import pyrosetta
    from pyrosetta.rosetta.core.scoring import ScoreType

    global _sfxn

    # Score to populate H-bonds
    _sfxn(pose)
    hbond_set = pose.get_hbonds()

    pdb_info = pose.pdb_info()
    interfacial_hbonds = []

    for i in range(1, hbond_set.nhbonds() + 1):
        hb = hbond_set.hbond(i)

        # Get donor info
        don_res_idx = hb.don_res()
        don_res = pose.residue(don_res_idx)
        don_h_idx = hb.don_hatm()  # Hydrogen atom index
        don_heavy_idx = don_res.atom_base(don_h_idx)  # Heavy atom bonded to H

        don_chain = pdb_info.chain(don_res_idx)
        don_resnum = pdb_info.number(don_res_idx)
        don_resname = don_res.name3()
        don_heavy_name = don_res.atom_name(don_heavy_idx).strip()
        don_h_name = don_res.atom_name(don_h_idx).strip()

        # Get acceptor info
        acc_res_idx = hb.acc_res()
        acc_res = pose.residue(acc_res_idx)
        acc_atm_idx = hb.acc_atm()

        acc_chain = pdb_info.chain(acc_res_idx)
        acc_resnum = pdb_info.number(acc_res_idx)
        acc_resname = acc_res.name3()
        acc_atom_name = acc_res.atom_name(acc_atm_idx).strip()

        # Check if interfacial
        is_interfacial = False
        if (don_chain == target_chain and acc_chain == motif_chain):
            is_interfacial = True
            direction = 'target_donates'
        elif (don_chain == motif_chain and acc_chain == target_chain):
            is_interfacial = True
            direction = 'motif_donates'

        if not is_interfacial:
            continue

        # Get coordinates
        don_h_xyz = don_res.xyz(don_h_idx)
        don_heavy_xyz = don_res.xyz(don_heavy_idx)
        acc_xyz = acc_res.xyz(acc_atm_idx)

        energy = hb.energy()

        interfacial_hbonds.append({
            'energy': energy,
            'direction': direction,
            'donor': {
                'chain': don_chain,
                'resnum': don_resnum,
                'resname': don_resname,
                'heavy_atom': don_heavy_name,
                'h_atom': don_h_name,
                'heavy_xyz': (don_heavy_xyz.x, don_heavy_xyz.y, don_heavy_xyz.z),
                'h_xyz': (don_h_xyz.x, don_h_xyz.y, don_h_xyz.z),
            },
            'acceptor': {
                'chain': acc_chain,
                'resnum': acc_resnum,
                'resname': acc_resname,
                'atom': acc_atom_name,
                'xyz': (acc_xyz.x, acc_xyz.y, acc_xyz.z),
            }
        })

    return interfacial_hbonds


def analyze_candidate(candidate, target_pdb_path):
    """
    Perform detailed H-bond analysis on a single candidate.

    Returns dict with analysis results and interfacial H-bonds.
    """
    import pyrosetta
    from pyrosetta.rosetta.core.import_pose import pose_from_pdbstring

    pdb_path = candidate['pdb_path']
    model_num = candidate['model_num']

    # Load target
    with open(target_pdb_path, 'r') as f:
        target_pdb_string = f.read()

    # Get specific model from candidate file
    models = split_multi_model_pdb(pdb_path)
    if model_num not in models:
        return None

    motif_pdb_string = models[model_num]

    # Combine into complex
    combined_pdb_string = target_pdb_string + motif_pdb_string

    # Load pose
    complex_pose = pyrosetta.Pose()
    pose_from_pdbstring(complex_pose, combined_pdb_string)

    if complex_pose.total_residue() == 0:
        return None

    # Get interfacial H-bonds
    interfacial_hbonds = get_interfacial_hbonds(complex_pose)

    # Sum interfacial H-bond energy
    total_interfacial_energy = sum(hb['energy'] for hb in interfacial_hbonds)

    return {
        'candidate': candidate,
        'interfacial_hbonds': interfacial_hbonds,
        'total_interfacial_energy': total_interfacial_energy,
        'num_interfacial_hbonds': len(interfacial_hbonds),
        'combined_pdb_string': combined_pdb_string
    }


def generate_pymol_script(analyses, output_pml, output_dir="spot_check_pdbs"):
    """
    Generate a comprehensive PyMOL script visualizing all candidates.

    Creates:
    - PDB files for each candidate complex
    - A .py script that loads all structures and draws H-bond visualizations
    - PyMOL scenes for each candidate for easy navigation
    """
    os.makedirs(output_dir, exist_ok=True)

    pml_lines = [
        "# Detailed H-bond spot-check visualization",
        "# Generated by detailed_energy_spot_checks.py",
        "",
        "from pymol.cgo import *",
        "from pymol import cmd",
        "from pymol import util",
        "",
        "# Settings",
        "cmd.set('stick_radius', 0.15)",
        "cmd.set('sphere_scale', 0.25)",
        "cmd.set('label_size', 14)",
        "cmd.set('label_color', 'white')",
        "cmd.set('float_labels', 'on')",
        "cmd.set('label_position', (0, 0, 3))",
        "",
    ]

    all_hbond_cgos = []  # Collect all CGO elements
    candidate_info = []  # Track info for scene creation

    for i, analysis in enumerate(analyses):
        if analysis is None:
            continue

        candidate = analysis['candidate']
        hbonds = analysis['interfacial_hbonds']

        # Create a unique name for this candidate
        pdb_basename = os.path.basename(candidate['pdb_path']).replace('.pdb', '')
        model_num = candidate['model_num']
        obj_name = f"spot_{i+1}_{pdb_basename}_m{model_num}"
        obj_name = obj_name.replace('-', '_').replace('.', '_')[:50]  # Sanitize

        # Save combined PDB
        pdb_filename = f"{obj_name}.pdb"
        pdb_path = os.path.join(output_dir, pdb_filename)
        with open(pdb_path, 'w') as f:
            f.write(analysis['combined_pdb_string'])

        # Track CGO indices for this candidate
        cgo_start_idx = len(all_hbond_cgos)

        # Add load command
        pml_lines.append(f"# Candidate {i+1}: {candidate['pdb_path']} model {model_num}")
        pml_lines.append(f"# delta_hbond: {candidate['delta_hbond']:.3f}, delta_fa_rep: {candidate['delta_fa_rep']:.3f}")
        pml_lines.append(f"# Interfacial H-bonds: {len(hbonds)}, Total energy: {analysis['total_interfacial_energy']:.3f}")
        pml_lines.append(f"cmd.load('{pdb_path}', '{obj_name}')")
        pml_lines.append(f"cmd.hide('everything', '{obj_name}')")
        pml_lines.append(f"cmd.show('cartoon', '{obj_name} and chain F')")
        pml_lines.append(f"cmd.show('sticks', '{obj_name} and chain Z')")
        pml_lines.append(f"cmd.color('lightblue', '{obj_name} and chain F')")
        pml_lines.append(f"cmd.color('orange', '{obj_name} and chain Z')")
        pml_lines.append("")

        # Track label names for this candidate
        label_names = []

        # Create CGO cylinders for H-bonds
        for j, hb in enumerate(hbonds):
            h_xyz = hb['donor']['h_xyz']
            acc_xyz = hb['acceptor']['xyz']
            energy = hb['energy']

            # Color based on energy strength (more negative = stronger = greener)
            # Range roughly -2 to 0
            intensity = min(1.0, max(0.0, -energy / 2.0))
            r = 1.0 - intensity
            g = intensity
            b = 0.2

            # Contract cylinder to central 80% of the H->acceptor vector
            dx = acc_xyz[0] - h_xyz[0]
            dy = acc_xyz[1] - h_xyz[1]
            dz = acc_xyz[2] - h_xyz[2]
            start_x = h_xyz[0] + dx * 0.1
            start_y = h_xyz[1] + dy * 0.1
            start_z = h_xyz[2] + dz * 0.1
            end_x = h_xyz[0] + dx * 0.9
            end_y = h_xyz[1] + dy * 0.9
            end_z = h_xyz[2] + dz * 0.9

            # CGO cylinder from H to acceptor (contracted to central 80%)
            cgo_name = f"hb_{obj_name}_{j}"
            cyl = (
                f"CYLINDER, {start_x:.3f}, {start_y:.3f}, {start_z:.3f}, "
                f"{end_x:.3f}, {end_y:.3f}, {end_z:.3f}, "
                f"0.08, {r:.2f}, {g:.2f}, {b:.2f}, {r:.2f}, {g:.2f}, {b:.2f}"
            )
            all_hbond_cgos.append(cyl)

            # Add label at midpoint
            mid_x = (h_xyz[0] + acc_xyz[0]) / 2
            mid_y = (h_xyz[1] + acc_xyz[1]) / 2
            mid_z = (h_xyz[2] + acc_xyz[2]) / 2

            # Create pseudoatom for label
            label_name = f"lbl_{obj_name}_{j}"
            label_names.append(label_name)
            donor_info = f"{hb['donor']['resname']}{hb['donor']['resnum']}"
            acc_info = f"{hb['acceptor']['resname']}{hb['acceptor']['resnum']}"
            label_text = f"{energy:.2f}"

            pml_lines.append(f"cmd.pseudoatom('{label_name}', pos=[{mid_x:.3f}, {mid_y:.3f}, {mid_z:.3f}])")
            pml_lines.append(f"cmd.hide('nonbonded', '{label_name}')")  # Hide the pseudoatom marker
            pml_lines.append(f"cmd.label('{label_name}', '\"{label_text}\"')")

        # Show interacting residues as sticks
        if hbonds:
            target_residues = set()
            motif_residues = set()
            for hb in hbonds:
                if hb['donor']['chain'] == 'F':
                    target_residues.add(hb['donor']['resnum'])
                else:
                    motif_residues.add(hb['donor']['resnum'])
                if hb['acceptor']['chain'] == 'F':
                    target_residues.add(hb['acceptor']['resnum'])
                else:
                    motif_residues.add(hb['acceptor']['resnum'])

            if target_residues:
                resi_sel = '+'.join(str(r) for r in target_residues)
                pml_lines.append(f"cmd.show('sticks', '{obj_name} and chain F and resi {resi_sel}')")

        pml_lines.append("")

        # Track info for scene creation
        cgo_end_idx = len(all_hbond_cgos)
        candidate_info.append({
            'index': i + 1,
            'obj_name': obj_name,
            'label_names': label_names,
            'cgo_indices': (cgo_start_idx, cgo_end_idx),
            'delta_hbond': candidate['delta_hbond'],
            'num_hbonds': len(hbonds),
            'pdb_basename': pdb_basename,
            'model_num': model_num
        })

    # Write all H-bond CGOs as a single object
    if all_hbond_cgos:
        pml_lines.append("# All interfacial H-bonds as CGO cylinders")
        pml_lines.append("hbond_cgo = [")
        for cgo in all_hbond_cgos:
            pml_lines.append(f"    {cgo},")
        pml_lines.append("]")
        pml_lines.append("cmd.load_cgo(hbond_cgo, 'interfacial_hbonds')")
        pml_lines.append("")
        pml_lines.append("# Color legend: Green = strong H-bond (< -1.5), Yellow = medium, Red = weak (> -0.5)")
        pml_lines.append("")

    # Create per-candidate CGO objects for scene visibility control
    pml_lines.append("# Per-candidate H-bond CGO objects")
    for info in candidate_info:
        start, end = info['cgo_indices']
        if end > start:
            cgo_obj_name = f"hb_{info['obj_name']}"
            pml_lines.append(f"{cgo_obj_name}_cgo = [")
            for idx in range(start, end):
                pml_lines.append(f"    {all_hbond_cgos[idx]},")
            pml_lines.append("]")
            pml_lines.append(f"cmd.load_cgo({cgo_obj_name}_cgo, '{cgo_obj_name}')")
    pml_lines.append("")

    # Create scenes for each candidate
    pml_lines.append("# Create scenes for each candidate")
    pml_lines.append("# Use PgUp/PgDn or scene buttons to navigate")
    pml_lines.append("")

    for info in candidate_info:
        scene_name = f"spot_{info['index']}"
        obj_name = info['obj_name']
        cgo_obj_name = f"hb_{obj_name}"

        pml_lines.append(f"# Scene {info['index']}: {info['pdb_basename']} model {info['model_num']}")
        pml_lines.append(f"# dHbond: {info['delta_hbond']:.3f}, {info['num_hbonds']} H-bonds")

        # Disable all objects first
        pml_lines.append("cmd.disable('spot_*')")
        pml_lines.append("cmd.disable('hb_spot_*')")
        pml_lines.append("cmd.disable('lbl_*')")
        pml_lines.append("cmd.disable('interfacial_hbonds')")

        # Enable the groups first, then the specific objects
        pml_lines.append("cmd.enable('spot_checks')")
        pml_lines.append("cmd.enable('hbond_cgos')")
        pml_lines.append("cmd.enable('hbond_labels')")

        # Enable this candidate's objects
        pml_lines.append(f"cmd.enable('{obj_name}')")
        pml_lines.append(f"cmd.enable('{cgo_obj_name}')")

        # Enable labels for this candidate
        for lbl in info['label_names']:
            pml_lines.append(f"cmd.enable('{lbl}')")

        # Zoom to this candidate
        pml_lines.append(f"cmd.zoom('{obj_name}', buffer=5)")

        # Ensure groups are enabled and coloring is applied before storing scene
        pml_lines.append("cmd.enable('spot_checks')")
        pml_lines.append("cmd.enable('hbond_cgos')")
        pml_lines.append("cmd.enable('hbond_labels')")
        pml_lines.append("util.cnc('all', _self=cmd)")

        # Store the scene
        pml_lines.append(f"cmd.scene('{scene_name}', 'store', message='{info['pdb_basename']} m{info['model_num']}: dHbond={info['delta_hbond']:.2f}')")
        pml_lines.append("")

    # Add summary table as comments
    pml_lines.append("# " + "="*70)
    pml_lines.append("# SUMMARY")
    pml_lines.append("# " + "="*70)
    pml_lines.append(f"# {'Candidate':<40} | {'dHbond':>8} | {'#HB':>4} | {'IntfE':>8}")
    pml_lines.append("# " + "-"*70)

    for i, analysis in enumerate(analyses):
        if analysis is None:
            continue
        c = analysis['candidate']
        short_name = os.path.basename(c['pdb_path'])[:35]
        pml_lines.append(
            f"# {short_name:<40} | {c['delta_hbond']:>8.3f} | "
            f"{analysis['num_interfacial_hbonds']:>4} | {analysis['total_interfacial_energy']:>8.3f}"
        )

    # Final display settings
    num_candidates = len([a for a in analyses if a])
    pml_lines.extend([
        "",
        "# Display settings",
        "cmd.bg_color('black')",
        "cmd.set('ray_shadows', 'off')",
        "cmd.set('specular', 0.2)",
        "cmd.set('scene_buttons', 1)",
        "cmd.set('scene_loop', 1)",
        "",
        "# Group all spot-check objects",
        "cmd.group('spot_checks', 'spot_*')",
        "cmd.group('hbond_labels', 'lbl_*')",
        "cmd.group('hbond_cgos', 'hb_spot_*')",
        "",
        "# Recall first scene",
        "cmd.scene('spot_1', 'recall')",
        "",
        f"print('Loaded {num_candidates} spot-check candidates')",
        "print('')",
        "print('NAVIGATION:')",
        "print('  PgUp/PgDn    - Previous/Next scene')",
        "print('  Scene buttons - Click numbered buttons at bottom')",
        "print('  scene spot_N  - Jump to specific scene (N=1-10)')",
        "print('')",
        "print('SHOW ALL:')",
        "print('  enable spot_checks; enable hbond_cgos; enable interfacial_hbonds')",
        "print('')",
    ])

    # Write PML file
    with open(output_pml, 'w') as f:
        f.write('\n'.join(pml_lines))

    print(f"Generated PyMOL script: {output_pml}")
    print(f"PDB files saved to: {output_dir}/")


def print_detailed_report(analyses):
    """Print a detailed text report of the H-bond analysis."""
    print("\n" + "="*80)
    print("DETAILED H-BOND SPOT-CHECK REPORT")
    print("="*80)

    for i, analysis in enumerate(analyses):
        if analysis is None:
            continue

        c = analysis['candidate']
        hbonds = analysis['interfacial_hbonds']

        print(f"\n--- Candidate {i+1} ---")
        print(f"File: {c['pdb_path']}")
        print(f"Model: {c['model_num']}")
        print(f"Delta Total: {c['delta_total']:.3f}")
        print(f"Delta fa_rep: {c['delta_fa_rep']:.3f}")
        print(f"Delta H-bond (from scan): {c['delta_hbond']:.3f}")
        print(f"Interfacial H-bond energy: {analysis['total_interfacial_energy']:.3f}")
        print(f"Number of interfacial H-bonds: {len(hbonds)}")

        if hbonds:
            print("\nInterfacial H-bonds:")
            print(f"  {'Donor':<25} -> {'Acceptor':<25} | {'Energy':>8}")
            print("  " + "-"*65)
            for hb in sorted(hbonds, key=lambda x: x['energy']):
                d = hb['donor']
                a = hb['acceptor']
                donor_str = f"{d['chain']}:{d['resname']}{d['resnum']}:{d['heavy_atom']}"
                acc_str = f"{a['chain']}:{a['resname']}{a['resnum']}:{a['atom']}"
                print(f"  {donor_str:<25} -> {acc_str:<25} | {hb['energy']:>8.3f}")
        else:
            print("\n  (No interfacial H-bonds detected)")


def main():
    parser = argparse.ArgumentParser(
        description="Detailed H-bond spot-check for validating energy analysis"
    )
    parser.add_argument(
        '--results', '-r',
        default='merged_results.json',
        help='Path to energy scan results JSON'
    )
    parser.add_argument(
        '--target', '-t',
        default='flagorigin_and_loop.pdb',
        help='Path to target PDB file'
    )
    parser.add_argument(
        '--output', '-o',
        default='spot_check.py',
        help='Output PyMOL Python script (.py)'
    )
    parser.add_argument(
        '--output-dir',
        default='spot_check_pdbs',
        help='Directory for output PDB files'
    )
    parser.add_argument(
        '--num-candidates', '-n',
        type=int,
        default=10,
        help='Number of candidates to spot-check'
    )
    parser.add_argument(
        '--max-fa-rep',
        type=float,
        default=50.0,
        help='Maximum delta_fa_rep to include (filter clashes)'
    )
    parser.add_argument(
        '--min-hbond',
        type=float,
        default=-1.0,
        help='Maximum delta_hbond to include (must be negative)'
    )
    parser.add_argument(
        '--params-dir',
        default='params',
        help='Directory containing .params files'
    )
    parser.add_argument(
        '--specific',
        nargs='+',
        metavar='PDB:MODEL',
        help='Analyze specific candidates (format: path/to/file.pdb:model_num)'
    )

    args = parser.parse_args()

    # Validate inputs
    if not os.path.exists(args.target):
        print(f"Error: Target file not found: {args.target}")
        sys.exit(1)

    # Handle specific candidates if provided
    if args.specific:
        candidates = []
        for spec in args.specific:
            if ':' not in spec:
                print(f"Error: Invalid format '{spec}'. Use path/to/file.pdb:model_num")
                sys.exit(1)
            pdb_path, model_str = spec.rsplit(':', 1)
            try:
                model_num = int(model_str)
            except ValueError:
                print(f"Error: Invalid model number '{model_str}' in '{spec}'")
                sys.exit(1)
            if not os.path.exists(pdb_path):
                print(f"Error: PDB file not found: {pdb_path}")
                sys.exit(1)
            candidates.append({
                'pdb_path': pdb_path,
                'model_num': model_num,
                'delta_total': 0.0,  # Unknown
                'delta_fa_rep': 0.0,
                'delta_hbond': 0.0,
                'original_id': 'manual'
            })
        print(f"Analyzing {len(candidates)} specific candidates")
    else:
        if not os.path.exists(args.results):
            print(f"Error: Results file not found: {args.results}")
            sys.exit(1)

        print(f"Loading results from: {args.results}")
        results = load_results(args.results)

        print(f"Selecting {args.num_candidates} diverse candidates...")
        print(f"  Filtering: delta_fa_rep < {args.max_fa_rep}, delta_hbond < {args.min_hbond}")
        candidates = select_diverse_candidates(
            results,
            num_candidates=args.num_candidates,
            max_fa_rep=args.max_fa_rep,
            min_hbond=args.min_hbond
        )

    if not candidates:
        print("No suitable candidates found. Try relaxing filter criteria.")
        sys.exit(1)

    print(f"Selected {len(candidates)} candidates:")
    for i, c in enumerate(candidates):
        if c.get('original_id') == 'manual':
            print(f"  {i+1}. {os.path.basename(c['pdb_path'])} model {c['model_num']} (manual)")
        else:
            print(f"  {i+1}. {os.path.basename(c['pdb_path'])} model {c['model_num']}: "
                  f"dHbond={c['delta_hbond']:.3f}, dFaRep={c['delta_fa_rep']:.3f}")

    print("\nInitializing PyRosetta...")
    init_pyrosetta(args.params_dir)

    print("\nAnalyzing candidates...")
    analyses = []
    for i, candidate in enumerate(candidates):
        print(f"  Analyzing {i+1}/{len(candidates)}: {os.path.basename(candidate['pdb_path'])} model {candidate['model_num']}")
        analysis = analyze_candidate(candidate, args.target)
        analyses.append(analysis)

        if analysis:
            print(f"    Found {analysis['num_interfacial_hbonds']} interfacial H-bonds "
                  f"(energy: {analysis['total_interfacial_energy']:.3f})")

    # Print detailed report
    print_detailed_report(analyses)

    # Generate PyMOL visualization
    print(f"\nGenerating PyMOL visualization...")
    generate_pymol_script(analyses, args.output, args.output_dir)

    print(f"\nDone! To visualize in PyMOL:")
    print(f"  pymol -r {args.output}")
    print(f"  Or from PyMOL: run {args.output}")


if __name__ == '__main__':
    main()
