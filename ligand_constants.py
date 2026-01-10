
# Atom Mappings and Param Definitions for Candidate Transformation
# Updated to match MEO, ACT, ACM, MMA, MGN, IMD, IND params files

# Mapping Residue -> Ligand Identity
RESIDUE_TO_LIGAND = {
    'SER': 'MEO', 'THR': 'MEO', 'TYR': 'MEO',
    'ASP': 'ACT', 'GLU': 'ACT',
    'ASN': 'ACM', 'GLN': 'ACM',
    'LYS': 'MMA',
    'ARG': 'MGN',
    'HIS': 'IMD',
    'TRP': 'IND'
}

# Mapping Residue Atoms -> Ligand Atoms
# These mappings convert amino acid sidechain atoms to the corresponding small molecule atoms
# Atoms not in the mapping are dropped (e.g., hydrogens on carbons that become methyls)
ATOM_MAPPING = {
    # MEO (methanol): C1-O1-HO with methyl H's (H1,H2,H3)
    'SER': {
        'CB': 'C1', 'OG': 'O1', 'HG': 'HO',
        # Drop CB hydrogens - will use implicit H in Rosetta
        '1HB': None, '2HB': None, 'HB2': None, 'HB3': None
    },
    'THR': {
        'CB': 'C1', 'OG1': 'O1', 'HG1': 'HO',
        # Drop other THR atoms
        'CG2': None, '1HG2': None, '2HG2': None, '3HG2': None,
        'HG21': None, 'HG22': None, 'HG23': None, 'HB': None
    },
    'TYR': {
        'CZ': 'C1', 'OH': 'O1', 'HH': 'HO',
        # Drop aromatic ring atoms
        'CG': None, 'CD1': None, 'CD2': None, 'CE1': None, 'CE2': None,
        'HD1': None, 'HD2': None, 'HE1': None, 'HE2': None
    },

    # ACT (acetate): C1 (methyl) - C2 (carboxylate) - O1, O2
    'ASP': {
        'CB': 'C1', 'CG': 'C2', 'OD1': 'O1', 'OD2': 'O2',
        # Drop CB hydrogens
        '1HB': None, '2HB': None, 'HB2': None, 'HB3': None
    },
    'GLU': {
        'CG': 'C1', 'CD': 'C2', 'OE1': 'O1', 'OE2': 'O2',
        # Drop CB and CG hydrogens
        'CB': None, '1HB': None, '2HB': None, 'HB2': None, 'HB3': None,
        '1HG': None, '2HG': None, 'HG2': None, 'HG3': None
    },

    # ACM (acetamide): C1 (methyl) - C2 (carbonyl) - O1, N1 - HN1, HN2
    'ASN': {
        'CB': 'C1', 'CG': 'C2', 'OD1': 'O1', 'ND2': 'N1',
        '1HD2': 'HN1', '2HD2': 'HN2', 'HD21': 'HN2', 'HD22': 'HN1',
        # Drop CB hydrogens
        '1HB': None, '2HB': None, 'HB2': None, 'HB3': None
    },
    'GLN': {
        'CG': 'C1', 'CD': 'C2', 'OE1': 'O1', 'NE2': 'N1',
        '1HE2': 'HN1', '2HE2': 'HN2', 'HE21': 'HN2', 'HE22': 'HN1',
        # Drop CB and CG hydrogens
        'CB': None, '1HB': None, '2HB': None, 'HB2': None, 'HB3': None,
        '1HG': None, '2HG': None, 'HG2': None, 'HG3': None
    },

    # MMA (methylammonium): C1 (methyl) - N1 - HN1, HN2, HN3
    'LYS': {
        'CE': 'C1', 'NZ': 'N1',
        '1HZ': 'HN1', '2HZ': 'HN2', '3HZ': 'HN3',
        'HZ1': 'HN1', 'HZ2': 'HN2', 'HZ3': 'HN3',
        # Drop other atoms
        'CD': None, '1HD': None, '2HD': None, 'HD2': None, 'HD3': None,
        'CG': None, '1HG': None, '2HG': None, 'HG2': None, 'HG3': None,
        'CB': None, '1HB': None, '2HB': None, 'HB2': None, 'HB3': None,
        '1HE': None, '2HE': None, 'HE2': None, 'HE3': None
    },

    # MGN (methylguanidinium): C1 (methyl) - NE - CZ - NH1, NH2
    'ARG': {
        'CD': 'C1', 'NE': 'NE', 'CZ': 'CZ', 'NH1': 'NH1', 'NH2': 'NH2',
        'HE': 'HE',
        '1HH1': 'HH11', '2HH1': 'HH12', 'HH11': 'HH11', 'HH12': 'HH12',
        '1HH2': 'HH21', '2HH2': 'HH22', 'HH21': 'HH21', 'HH22': 'HH22',
        # Drop other atoms
        'CG': None, '1HG': None, '2HG': None, 'HG2': None, 'HG3': None,
        'CB': None, '1HB': None, '2HB': None, 'HB2': None, 'HB3': None,
        '1HD': None, '2HD': None, 'HD2': None, 'HD3': None
    },

    # IMD (imidazole): CG - ND1 - CE1 - NE2 - CD2 (ring)
    'HIS': {
        'CG': 'CG', 'ND1': 'ND1', 'CD2': 'CD2', 'CE1': 'CE1', 'NE2': 'NE2',
        'HE2': 'HE2', 'HE1': 'HE1', 'HD2': 'HD2', 'HD1': 'HD1',
        # Drop CB
        'CB': None, '1HB': None, '2HB': None, 'HB2': None, 'HB3': None
    },

    # IND (indole): full indole ring system
    'TRP': {
        'CG': 'CG', 'CD1': 'CD1', 'CD2': 'CD2', 'NE1': 'NE1', 'CE2': 'CE2',
        'CE3': 'CE3', 'CZ2': 'CZ2', 'CZ3': 'CZ3', 'CH2': 'CH2',
        'HE1': 'HE1', 'HD1': 'HD1', 'HE3': 'HE3', 'HZ2': 'HZ2', 'HZ3': 'HZ3', 'HH2': 'HH2',
        # Drop CB
        'CB': None, '1HB': None, '2HB': None, 'HB2': None, 'HB3': None
    }
}

# Param definitions are now in the separate .params files in the params/ directory
# These are kept here for reference but not used by the transformation script
LIGAND_PARAMS = {
    'MEO': [
        ('C1', 'CH3', -0.05), ('O1', 'OH', -0.66), ('HO', 'Hpol', 0.43),
        ('H1', 'Hapo', 0.09), ('H2', 'Hapo', 0.09), ('H3', 'Hapo', 0.10)
    ],
    'ACT': [
        ('C1', 'CH3', -0.27), ('C2', 'COO', 0.62), ('O1', 'OOC', -0.76), ('O2', 'OOC', -0.76),
        ('H1', 'Hapo', 0.09), ('H2', 'Hapo', 0.09), ('H3', 'Hapo', -0.01)
    ],
    'ACM': [
        ('C1', 'CH3', -0.10), ('C2', 'CNH2', 0.55), ('O1', 'ONH2', -0.55), ('N1', 'NH2O', -0.62),
        ('HN1', 'Hpol', 0.32), ('HN2', 'Hpol', 0.13),
        ('H1', 'Hapo', 0.09), ('H2', 'Hapo', 0.09), ('H3', 'Hapo', 0.09)
    ],
    'MMA': [
        ('C1', 'CH3', 0.13), ('N1', 'Nlys', -0.30),
        ('HN1', 'Hpol', 0.33), ('HN2', 'Hpol', 0.33), ('HN3', 'Hpol', 0.24),
        ('H1', 'Hapo', 0.09), ('H2', 'Hapo', 0.09), ('H3', 'Hapo', 0.09)
    ],
    'MGN': [
        ('C1', 'CH3', 0.05), ('NE', 'NtrR', -0.35), ('CZ', 'aroC', 0.64),
        ('NH1', 'Narg', -0.40), ('NH2', 'Narg', -0.40),
        ('HE', 'Hpol', 0.30), ('HH11', 'Hpol', 0.30), ('HH12', 'Hpol', 0.30),
        ('HH21', 'Hpol', 0.30), ('HH22', 'Hpol', 0.28),
        ('H1', 'Hapo', 0.06), ('H2', 'Hapo', 0.06), ('H3', 'Hapo', 0.06)
    ],
    'IMD': [
        ('CG', 'CH0', 0.05), ('ND1', 'Nhis', -0.55), ('CD2', 'aroC', -0.05),
        ('CE1', 'aroC', 0.20), ('NE2', 'Ntrp', -0.30),
        ('HE2', 'Hpol', 0.30), ('HE1', 'Hapo', 0.10), ('HD2', 'Hapo', 0.10), ('HG', 'Hapo', 0.15)
    ],
    'IND': [
        ('CG', 'CH0', -0.05), ('CD1', 'aroC', 0.05), ('CD2', 'CH0', 0.00),
        ('NE1', 'Ntrp', -0.40), ('CE2', 'CH0', 0.10),
        ('CE3', 'aroC', -0.10), ('CZ2', 'aroC', -0.10), ('CZ3', 'aroC', -0.10), ('CH2', 'aroC', -0.10),
        ('HE1', 'Hpol', 0.35), ('HD1', 'Haro', 0.10), ('HG', 'Haro', 0.10),
        ('HE3', 'Haro', 0.10), ('HZ2', 'Haro', 0.10), ('HZ3', 'Haro', 0.10), ('HH2', 'Haro', -0.05)
    ]
}

# Bond definitions (kept for reference)
LIGAND_BONDS = {
    'MEO': [('C1','O1'), ('O1','HO'), ('C1','H1'), ('C1','H2'), ('C1','H3')],
    'ACT': [('C1','C2'), ('C2','O1'), ('C2','O2'), ('C1','H1'), ('C1','H2'), ('C1','H3')],
    'ACM': [('C1','C2'), ('C2','O1'), ('C2','N1'), ('N1','HN1'), ('N1','HN2'), ('C1','H1'), ('C1','H2'), ('C1','H3')],
    'MMA': [('C1','N1'), ('N1','HN1'), ('N1','HN2'), ('N1','HN3'), ('C1','H1'), ('C1','H2'), ('C1','H3')],
    'MGN': [('C1','NE'), ('NE','CZ'), ('NE','HE'), ('CZ','NH1'), ('CZ','NH2'),
            ('NH1','HH11'), ('NH1','HH12'), ('NH2','HH21'), ('NH2','HH22'),
            ('C1','H1'), ('C1','H2'), ('C1','H3')],
    'IMD': [('CG','ND1'), ('CG','CD2'), ('ND1','CE1'), ('CD2','NE2'), ('CE1','NE2'),
            ('NE2','HE2'), ('CE1','HE1'), ('CD2','HD2'), ('CG','HG')],
    'IND': [('CG','CD1'), ('CG','CD2'), ('CD1','NE1'), ('NE1','CE2'), ('CE2','CD2'),
            ('CD2','CE3'), ('CE2','CZ2'), ('CE3','CZ3'), ('CZ2','CH2'), ('CZ3','CH2'),
            ('NE1','HE1'), ('CD1','HD1'), ('CG','HG'), ('CE3','HE3'), ('CZ2','HZ2'), ('CZ3','HZ3'), ('CH2','HH2')]
}
