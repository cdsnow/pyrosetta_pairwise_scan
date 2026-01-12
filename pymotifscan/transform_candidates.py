#!/usr/bin/env python3
"""
Transform candidate PDB files to use small molecule fragments.

For each amino acid fragment:
1. If it's a "tip fragment" (sidechain only, no backbone), convert to small molecule
2. If it's a full residue truncated at CB (has backbone but only CB as sidechain), rename to ALA
3. Keep water (HOH) and other residues as-is

Output goes to candidates_v2/ directory.
"""

import os
import glob
import re
import ligand_constants


def is_backbone_atom(name):
    """Check if atom is backbone or CB."""
    return name in ['N', 'CA', 'C', 'O', 'CB', 'H', 'HA', '1HA', '2HA', 'OXT',
                    '1H', '2H', '3H', 'H1', 'H2', 'H3']


def has_backbone(atoms):
    """Check if residue has backbone atoms."""
    atom_names = {a['name'] for a in atoms}
    return 'CA' in atom_names and 'N' in atom_names


def is_truncated_at_cb(resname, atoms):
    """Check if residue has backbone + only CB (no other sidechain)."""
    if resname in ['GLY', 'ALA']:
        return False
    atom_names = {a['name'] for a in atoms}
    sidechain_atoms = atom_names - {'N', 'CA', 'C', 'O', 'CB', 'H', 'HA', '1HA', '2HA',
                                     'OXT', '1H', '2H', '3H', 'H1', 'H2', 'H3',
                                     '1HB', '2HB', '3HB', 'HB', 'HB1', 'HB2', 'HB3'}
    return len(sidechain_atoms) == 0 and 'CB' in atom_names


def format_atom_name(name):
    """Format atom name for PDB (4 chars, proper alignment)."""
    if len(name) == 4:
        return name
    elif name[0].isdigit():
        return f"{name:<4}"
    else:
        return f" {name:<3}"


def transform_pdb_line(line, new_resname, atom_mapping):
    """Transform a PDB ATOM line with new residue name and atom name mapping."""
    atom_name = line[12:16].strip()

    # Check if atom should be dropped
    if atom_name in atom_mapping:
        new_atom = atom_mapping[atom_name]
        if new_atom is None:
            return None  # Drop this atom
    else:
        # Atom not in mapping - keep original name for tip fragments
        new_atom = atom_name

    # Format the new line
    atom_field = format_atom_name(new_atom)
    new_line = (
        line[:12] +
        atom_field +
        line[16:17] +
        f"{new_resname:>3}" +
        line[20:]
    )
    return new_line


def process_residue(resname, atoms, lines):
    """
    Process a residue and return transformed lines.

    Returns: list of transformed PDB lines
    """
    # Water - keep as is
    if resname in ['HOH', 'WAT']:
        return lines

    # Check if this is a backbone-containing residue
    if has_backbone(atoms):
        # Check if truncated at CB -> convert to ALA
        if is_truncated_at_cb(resname, atoms):
            return [transform_pdb_line(l, 'ALA', {}) for l in lines if transform_pdb_line(l, 'ALA', {}) is not None]
        else:
            # Full residue with sidechain - keep as is
            return lines

    # Tip fragment (no backbone) - convert to small molecule
    target_lig = ligand_constants.RESIDUE_TO_LIGAND.get(resname)
    atom_mapping = ligand_constants.ATOM_MAPPING.get(resname, {})

    if target_lig:
        transformed = []
        for line in lines:
            new_line = transform_pdb_line(line, target_lig, atom_mapping)
            if new_line is not None:
                transformed.append(new_line)
        return transformed
    else:
        # Unknown residue - keep as is
        return lines


def parse_target_from_signature(signature):
    """Parse target residue from signature string like 'F6_ASP_SC'."""
    match = re.match(r'^([A-Z])(\d+)_', signature)
    if match:
        return match.group(1), match.group(2)
    return None, None


def process_pdb_file(input_path, output_path, target_chain=None, target_resnum=None):
    """Process a single PDB file."""

    with open(input_path, 'r') as f:
        lines = f.readlines()

    output_lines = []
    current_res_key = None
    current_res_atoms = []
    current_res_lines = []
    current_resname = None

    for line in lines:
        if line.startswith('ATOM') or line.startswith('HETATM'):
            # Parse atom info
            atom_name = line[12:16].strip()
            resname = line[17:20].strip()
            chain = line[21]
            resnum = line[22:27]  # includes insertion code

            res_key = (chain, resnum)

            if res_key != current_res_key:
                # Process previous residue
                if current_res_key is not None:
                    # Skip target residue
                    if not (current_res_key[0] == target_chain and current_res_key[1].strip() == target_resnum):
                        transformed = process_residue(current_resname, current_res_atoms, current_res_lines)
                        output_lines.extend(transformed)

                # Start new residue
                current_res_key = res_key
                current_resname = resname
                current_res_atoms = []
                current_res_lines = []

            # Parse atom coordinates
            try:
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
            except:
                x, y, z = 0, 0, 0

            current_res_atoms.append({
                'name': atom_name,
                'x': x, 'y': y, 'z': z
            })
            current_res_lines.append(line)

        elif line.startswith('MODEL') or line.startswith('ENDMDL') or line.startswith('END'):
            # Process any pending residue before MODEL/ENDMDL/END
            if current_res_key is not None:
                # Skip target residue
                if not (current_res_key[0] == target_chain and current_res_key[1].strip() == target_resnum):
                    transformed = process_residue(current_resname, current_res_atoms, current_res_lines)
                    output_lines.extend(transformed)
                current_res_key = None
                current_res_atoms = []
                current_res_lines = []
                current_resname = None

            output_lines.append(line)
        else:
            # REMARK and other lines - keep as is
            output_lines.append(line)

    # Process any remaining residue
    if current_res_key is not None:
        if not (current_res_key[0] == target_chain and current_res_key[1].strip() == target_resnum):
            transformed = process_residue(current_resname, current_res_atoms, current_res_lines)
            output_lines.extend(transformed)

    # Write output
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    with open(output_path, 'w') as f:
        f.writelines(output_lines)


def main():
    import argparse
    parser = argparse.ArgumentParser(description='Transform candidate PDB files to use small molecule fragments')
    parser.add_argument('--indir', default='candidates', help='Input directory')
    parser.add_argument('--outdir', default='candidates_transformed', help='Output directory')
    args = parser.parse_args()

    # Find all PDB files
    pdb_files = glob.glob(os.path.join(args.indir, '**/*.pdb'), recursive=True)

    print(f"Found {len(pdb_files)} PDB files to process")

    for i, input_path in enumerate(pdb_files):
        if (i + 1) % 500 == 0:
            print(f"  Processed {i+1}/{len(pdb_files)} files...")

        # Construct output path
        rel_path = os.path.relpath(input_path, args.indir)
        output_path = os.path.join(args.outdir, rel_path)

        # Extract target from signature (first part of rel_path)
        signature = rel_path.split(os.sep)[0]
        target_chain, target_resnum = parse_target_from_signature(signature)

        process_pdb_file(input_path, output_path, target_chain, target_resnum)

    print(f"Done! Processed {len(pdb_files)} files.")
    print(f"Output written to: {args.outdir}")


if __name__ == '__main__':
    main()
