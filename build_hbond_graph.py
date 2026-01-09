import sys
import os
import math
import networkx as nx
from collections import defaultdict
import pyrosetta
from pyrosetta import rosetta

# Residue Definitions defining potential Donors and Acceptors
# BB always implies N (donor) and O (acceptor), except PRO (no N donor)
# Sidechains defined below.
# Format: 'sc_donors': [atom_names], 'sc_acceptors': [atom_names]
# Note: For HIS, we will check protonation dynamically if possible, or just add both potentials if ambiguous.
# But better to check PDB for H atoms.
RESIDUE_SPECS = {
    'ALA': {},
    'GLY': {},
    'VAL': {},
    'LEU': {},
    'ILE': {},
    'PHE': {},
    'PRO': {'no_bb_donor': True},
    'SER': {'sc_donors': ['OG'], 'sc_acceptors': ['OG']},
    'THR': {'sc_donors': ['OG1'], 'sc_acceptors': ['OG1']},
    'TYR': {'sc_donors': ['OH'], 'sc_acceptors': ['OH']},
    'ASP': {'sc_acceptors': ['OD1', 'OD2']},
    'GLU': {'sc_acceptors': ['OE1', 'OE2']},
    'ASN': {'sc_donors': ['ND2'], 'sc_acceptors': ['OD1']},
    'GLN': {'sc_donors': ['NE2'], 'sc_acceptors': ['OE1']},
    'LYS': {'sc_donors': ['NZ']},
    'ARG': {'sc_donors': ['NE', 'NH1', 'NH2']},
    'TRP': {'sc_donors': ['NE1']},
    'CYS': {'sc_donors': ['SG'], 'sc_acceptors': ['SG']}, # Weak, but possible
    'MET': {'sc_acceptors': ['SD']},
    'HIS': {'sc_donors': ['ND1', 'NE2'], 'sc_acceptors': ['ND1', 'NE2']}, # Dynamic handling needed
    'HOH': {'sc_donors': ['O'], 'sc_acceptors': ['O']} # Special handling for water
}

def normalize_resnum(r):
    # Try to convert to int to handle "   3" vs "0003"
    try:
        return str(int(r))
    except ValueError:
        return r.strip()

class PDBMap:
    def __init__(self, pdb_file):
        self.atoms = {} # (chain, resnum, atomname) -> (x, y, z)
        self.residues = [] # List of (chain, resnum, resname)
        self.polar_hydrogens = defaultdict(list) # (chain, resnum, heavy_atom) -> [h_atom_coords]
        self.load(pdb_file)
        
    def load(self, pdb_file):
        with open(pdb_file, 'r') as f:
            prev_res = None
            for line in f:
                if line.startswith("ATOM") or line.startswith("HETATM"):
                    chain = line[21]
                    resnum_raw = line[22:27].strip() # 4 digits + insertion code
                    resnum = normalize_resnum(resnum_raw)
                    resname = line[17:20].strip()
                    name = line[12:16].strip()
                    try:
                        x = float(line[30:38])
                        y = float(line[38:46])
                        z = float(line[46:54])
                    except ValueError:
                        continue
                    
                    self.atoms[(chain, resnum, name)] = (x, y, z)
                    
                    current_res = (chain, resnum, resname)
                    if current_res != prev_res:
                        self.residues.append(current_res)
                        prev_res = current_res
                    
                    # Store Hydrogens associated with heavy atoms
                    # Logic: Hydrogens usually start with H.
                    # We need to map H to its parent. PDB convention:
                    
                    is_hydrogen = False
                    if name.startswith("H"): is_hydrogen = True
                    # PDB convention for hydrogens: often " HA "," HB " etc. or "1H  "
                    # If stripped name starts with H calculation above is correct.
                    # Some hydrogens start with digit e.g. "1HG1" -> stripped "1HG1"
                    if not is_hydrogen and name[0].isdigit() and (len(name)>1 and name[1]=='H'):
                        is_hydrogen = True
                        
                    if is_hydrogen:
                         self._map_hydrogen(chain, resnum, name, x, y, z, resname)
                         # Debug first few
                         if len(self.polar_hydrogens) < 2:
                             print(f"DEBUG: Mapped H {name} to residue {resname} {resnum}")

    def _map_hydrogen(self, chain, resnum, h_name, x, y, z, resname):
        # Heuristic mapping of H to heavy atom
        parent = None
        # Backbone Amide
        if h_name == 'H': parent = 'N'
        
        # Sidechain heuristics based on name
        elif resname == 'LYS' and 'HZ' in h_name: parent = 'NZ'
        elif resname == 'ARG':
            if 'HE' in h_name: parent = 'NE'
            if 'HH1' in h_name: parent = 'NH1'
            if 'HH2' in h_name: parent = 'NH2'
        elif resname == 'ASN' and 'HD2' in h_name: parent = 'ND2'
        elif resname == 'GLN' and 'HE2' in h_name: parent = 'NE2'
        elif resname == 'SER' and 'HG' in h_name: parent = 'OG'
        elif resname == 'THR' and 'HG1' in h_name: parent = 'OG1'
        elif resname == 'TYR' and 'HH' in h_name: parent = 'OH'
        elif resname == 'TRP' and 'HE1' in h_name: parent = 'NE1'
        elif resname == 'HIS':
            if 'HD1' in h_name: parent = 'ND1'
            if 'HE2' in h_name: parent = 'NE2'
        elif resname == 'HOH': parent = 'O'
        elif resname == 'CYS' and 'HG' in h_name: parent = 'SG'
        
        if parent:
            self.polar_hydrogens[(chain, resnum, parent)].append((x, y, z))
        else:
            # Fallback: Distance-based assignment
            # Assign to closest heavy atom in the same residue
            # This handles cases like 1HH1, 2HH1 etc where naming is messy
            self.polar_hydrogens[(chain, resnum, 'GENERIC')].append((x, y, z))

    def get_closest_heavy_atom(self, chain, resnum, h_xyz):
        # Helper to find closest heavy atom in residue
        # Not efficiently implemented (requires iterating atoms), skipping for now.
        # Instead, relying on 'GENERIC' bucket in get_closest_h
        pass

    def get_atom(self, chain, resnum, atomname):
        return self.atoms.get((chain, resnum, atomname))

    def get_closest_h(self, chain, resnum, donor_heavy, acceptor_xyz):
        # Return the H attached to donor_heavy that is closest to acceptor_xyz
        
        # 1. Try specific parent mapping
        candidates = self.polar_hydrogens.get((chain, resnum, donor_heavy), [])
        
        # 2. Add generic candidates (unmapped hydrogens)
        candidates.extend(self.polar_hydrogens.get((chain, resnum, 'GENERIC'), []))
        
        # 3. Special case for Backbone: if donor is N, include 'H'
        if donor_heavy == 'N':
             candidates.extend(self.polar_hydrogens.get((chain, resnum, 'N'), []))

        if not candidates:
             # Fallback: No hydrogens known. Return heavy atom coord.
             xyz = self.get_atom(chain, resnum, donor_heavy)
             if not xyz:
                 # Debug why heavy atom lookup failed
                 print(f"DEBUG: Heavy Atom Lookup Failed for Key: {(chain, resnum, donor_heavy)}")
             return xyz
        
        # Filter candidates: Only those within bonding distance of donor_heavy (approx 1.0A)
        # This ensures we don't pick a random H from the other side of the residue
        donor_xyz = self.atoms.get((chain, resnum, donor_heavy))
        valid_candidates = []
        if donor_xyz:
            dx, dy, dz = donor_xyz
            for (hx, hy, hz) in candidates:
                d2 = (hx-dx)**2 + (hy-dy)**2 + (hz-dz)**2
                if d2 < 1.5**2: # generous 1.5A cutoff for H-X bond
                    valid_candidates.append((hx, hy, hz))
        else:
            valid_candidates = candidates # If we can't check distance, take them all
        
        if not valid_candidates:
            # Fallback 1: Use any candidate (ignore distance check)
            if candidates:
                valid_candidates = candidates
            else:
                 # Fallback 2: No hydrogens found mappable to this donor.
                 # Return the heavy atom coordinates themselves so at least a cylinder appears (zero length start, or just start from heavy)
                 # User wants cylinders. Better to start from Heavy than nothing.
                 return self.get_atom(chain, resnum, donor_heavy)
            
        best_h = None
        min_d2 = float('inf')
        ax, ay, az = acceptor_xyz
        for (hx, hy, hz) in valid_candidates:
            d2 = (hx-ax)**2 + (hy-ay)**2 + (hz-az)**2
            if d2 < min_d2:
                min_d2 = d2
                best_h = (hx, hy, hz)
        return best_h

def build_graph(pdb_map, hb2_file='2ij2.opt.hb2'):
    G = nx.Graph()
    
    # 1. Create Nodes
    # Mapping: (chain, resnum, atomname, role) -> NodeID
    # Role: 'donor' or 'acceptor'
    
    node_lookup = {} # (chain, resnum_int_str, atom, role) -> node_id

    # Let's iterate residues
    for i, (chain, resnum, resname) in enumerate(pdb_map.residues):
        spec = RESIDUE_SPECS.get(resname, {})
        norm_resnum = normalize_resnum(resnum)
        
        # Helper to add node and index it
        def add_node_to_graph(node_id, **kwargs):
            G.add_node(node_id, **kwargs)
            # Index it
            key = (kwargs['chain'], normalize_resnum(kwargs['resnum']), kwargs['atom'], kwargs['role'])
            node_lookup[key] = node_id
            
            # Special aliases (BB N/O)
            if kwargs['type'] == 'BB_Donor' and kwargs['atom'] == 'N':
                 node_lookup[(kwargs['chain'], norm_resnum, 'N', 'donor')] = node_id
            if kwargs['type'] == 'BB_Acceptor' and kwargs['atom'] == 'O':
                 node_lookup[(kwargs['chain'], norm_resnum, 'O', 'acceptor')] = node_id
        
        # Backbone Nodes
        if not spec.get('no_bb_donor') and resname not in ['HOH']:
            node_id = f"{chain}_{resnum}_BB_Donor"
            add_node_to_graph(node_id, 
                       chain=chain, resnum=resnum, resname=resname, 
                       type="BB_Donor", atom="N", role="donor")
            
        if resname not in ['HOH']:
            node_id = f"{chain}_{resnum}_BB_Acceptor"
            add_node_to_graph(node_id,
                       chain=chain, resnum=resnum, resname=resname,
                       type="BB_Acceptor", atom="O", role="acceptor")
            
        # Sidechain Nodes
        sc_donors = spec.get('sc_donors', [])
        sc_acceptors = spec.get('sc_acceptors', [])
        
        # Refine HIS
        if resname == 'HIS':
            # Check if HD1 exists -> ND1 is donor
            if pdb_map.polar_hydrogens.get((chain, resnum, 'ND1')):
                add_node_to_graph(f"{chain}_{resnum}_ND1_Donor", 
                           chain=chain, resnum=resnum, resname=resname,
                           type="SC_Donor", atom="ND1", role="donor")
            else:
                 add_node_to_graph(f"{chain}_{resnum}_ND1_Acceptor", 
                           chain=chain, resnum=resnum, resname=resname,
                           type="SC_Acceptor", atom="ND1", role="acceptor")

            if pdb_map.polar_hydrogens.get((chain, resnum, 'NE2')):
                add_node_to_graph(f"{chain}_{resnum}_NE2_Donor",
                           chain=chain, resnum=resnum, resname=resname,
                           type="SC_Donor", atom="NE2", role="donor")
            else:
                add_node_to_graph(f"{chain}_{resnum}_NE2_Acceptor",
                           chain=chain, resnum=resnum, resname=resname,
                           type="SC_Acceptor", atom="NE2", role="acceptor")
        elif resname == 'HOH':
            # Water
            add_node_to_graph(f"{chain}_{resnum}_Water_Donor", 
                       chain=chain, resnum=resnum, resname=resname,
                       type="Water_Donor", atom="O", role="donor")
            add_node_to_graph(f"{chain}_{resnum}_Water_Acceptor", 
                       chain=chain, resnum=resnum, resname=resname,
                       type="Water_Acceptor", atom="O", role="acceptor")
        else:
            # Standard SC
            for atom in sc_donors:
                node_id = f"{chain}_{resnum}_{atom}_{'Donor'}"
                add_node_to_graph(node_id, chain=chain, resnum=resnum, resname=resname,
                           type="SC_Donor", atom=atom, role="donor")
            for atom in sc_acceptors:
                node_id = f"{chain}_{resnum}_{atom}_{'Acceptor'}"
                add_node_to_graph(node_id, chain=chain, resnum=resnum, resname=resname,
                           type="SC_Acceptor", atom=atom, role="acceptor")

    # 2. Add Intra-residue and Contiguous Edges
    # Iterate nodes to link them
    # Group by residue
    nodes_by_res = defaultdict(list)
    for n, data in G.nodes(data=True):
        nodes_by_res[(data['chain'], data['resnum'])].append(n)
        
    for (chain, resnum), nodes in nodes_by_res.items():
        # Intra-residue edges
        for i in range(len(nodes)):
            for j in range(i+1, len(nodes)):
                G.add_edge(nodes[i], nodes[j], type="intra_residue")
    
    # Contiguous residues
    for i in range(len(pdb_map.residues) - 1):
        c1, r1, n1 = pdb_map.residues[i]
        c2, r2, n2 = pdb_map.residues[i+1]
        
        if c1 == c2:
            u = f"{c1}_{r1}_BB_Acceptor"
            v = f"{c2}_{r2}_BB_Donor"
            if G.has_node(u) and G.has_node(v):
                G.add_edge(u, v, type="polymer")

    # 3. Parse HBPLUS and add Hydrogen Bond Edges
    lines = open(hb2_file).readlines()
    
    # Debug counter
    failures = 0
    success = 0
    
    for line_idx, line in enumerate(lines[8:]): 
        if len(line) < 35: continue
        
        # Parse Donor
        d_chain = line[0]        
        d_resnum_raw = line[1:5]
        d_resnum = normalize_resnum(d_resnum_raw)
        
        d_resname = line[6:9] # Adjusted from 7:10
        d_atom = line[10:14].strip() # Adjusted from 11:15
        
        # Parse Acceptor
        a_chain = line[14] # From 20
        a_resnum_raw = line[15:19] # From 21:25
        a_resnum = normalize_resnum(a_resnum_raw)
        
        a_resname = line[20:23] # From 27:30
        a_atom = line[24:28].strip() # From 31:35
        
        dist_da = float(line[29:33]) # From 35:39
        cat = line[34:36] # From 40:42 
        
        # Find Nodes
        key_d = (d_chain, d_resnum, d_atom, "donor")
        key_a = (a_chain, a_resnum, a_atom, "acceptor")
        
        donor_node = node_lookup.get(key_d)
        acceptor_node = node_lookup.get(key_a)
        
        if donor_node and acceptor_node:
            success += 1
            # Need coords for cylinder
            # Donor Heavy XYZ
            d_xyz = pdb_map.get_atom(d_chain, d_resnum, d_atom)
            
            # Acceptor Heavy XYZ
            a_xyz = pdb_map.get_atom(a_chain, a_resnum, a_atom)
            
            # Find H position
            h_xyz = None
            if d_xyz and a_xyz:
                h_xyz = pdb_map.get_closest_h(d_chain, d_resnum, d_atom, a_xyz)
            
            if not h_xyz or not d_xyz or not a_xyz:
                # Debug why coord missing
                # only print if really broken
                 pass
            
            # Geometric Sanity Check (Fix for Phantom Edges)
            if d_xyz and a_xyz:
                 dx, dy, dz = d_xyz[0]-a_xyz[0], d_xyz[1]-a_xyz[1], d_xyz[2]-a_xyz[2]
                 geo_dist = math.sqrt(dx*dx + dy*dy + dz*dz)
                 if geo_dist > 4.0:
                      failures += 1
                      if failures <= 10: # Print first few
                          print(f"DEBUG: Skipping Phantom Edge (d={geo_dist:.1f}A > 4.0A). {d_chain}:{d_resnum}:{d_atom} -> {a_chain}:{a_resnum}:{a_atom}")
                      continue

            G.add_edge(donor_node, acceptor_node, 
                       type="hbond", 
                       category=cat, 
                       distance=dist_da,
                       h_xyz=h_xyz,
                       a_xyz=a_xyz) 
        else:
            failures += 1
            if failures < 10:
                print(f"DEBUG: Missed Edge. D: {key_d} -> {donor_node}, A: {key_a} -> {acceptor_node}. Line: {line.strip()}")
                # Is matching partial?
                if not donor_node:
                    print(f"  -> Donor Node Mismatch. Keys available for {d_chain}:{d_resnum}: {[k for k in node_lookup.keys() if k[0]==d_chain and k[1]==d_resnum]}")

    print(f"HBPLUS Parsing: {success} edges added, {failures} lines failed (or header/irrelevant).")
    return G


def normalize_resnum(r):
    # Try to convert to int to handle "   3" vs "0003"
    try:
        return str(int(r))
    except ValueError:
        return r.strip()

def find_node(G, chain, resnum, atom, role):
    # Search for node matching the signature
    target_resnum = normalize_resnum(resnum)
    
    # 1. Try exact atomic match (for SC atoms, atom name is in ID or attr)
    for n, d in G.nodes(data=True):
        if d['chain'] == chain and normalize_resnum(d['resnum']) == target_resnum and d['role'] == role:
            # Check atom match
            if d['atom'] == atom:
                return n
            
            # Special Handling for Backbone
            if role == 'donor' and atom == 'N' and d['type'] == 'BB_Donor':
                return n
            if role == 'acceptor' and atom == 'O' and d['type'] == 'BB_Acceptor':
                return n
                
            # Special Handling for Waters
            if d['resname'] == 'HOH' and atom == 'O':
                 # Determine if donor or acceptor node based on role arg
                 if role == 'donor' and d['type'] == 'Water_Donor': return n
                 if role == 'acceptor' and d['type'] == 'Water_Acceptor': return n
                 
            # Dynamic HIS mapping: hbplus says "ND1" but we might have "SC_Donor" or "SC_Acceptor" 
            # assigned to "ND1". If the atom attr matches, we are good.
            # d['atom'] should match 'ND1'
    
    return None

def get_hbond_type(u_data, v_data):
    # Determine type: wat_wat, aa_wat, aa_aa_bb_bb, aa_aa_sc_sc, aa_aa_bb_sc
    u_wat = (u_data['resname'] == 'HOH')
    v_wat = (v_data['resname'] == 'HOH')
    
    if u_wat and v_wat: return 'wat_wat'
    if u_wat or v_wat: return 'aa_wat'
    
    # Both AA
    # Check roles BB vs SC
    # Node types: BB_Donor, BB_Acceptor vs SC_...
    u_bb = 'BB' in u_data['type']
    v_bb = 'BB' in v_data['type']
    
    if u_bb and v_bb: return 'aa_aa_bb_bb'
    if not u_bb and not v_bb: return 'aa_aa_sc_sc'
    return 'aa_aa_bb_sc' # Mixed

def add_rosetta_energies(G, pdb_file):
    print("Initializing PyRosetta...")
    # Silence output if possible or keep std
    pyrosetta.init("-ex1 -ex2 -use_input_sc -flip_HNQ -ignore_unrecognized_res -ignore_waters 0")
    
    print(f"Loading Pose from {pdb_file}...")
    pose = pyrosetta.pose_from_pdb(pdb_file)
    
    # Score pose to populate energies (and hbonds)
    scorefxn = pyrosetta.get_fa_scorefxn()
    total_score = scorefxn(pose)
    print(f"Total Rosetta Score: {total_score:.3f} REU")
    
    # Analyze HBonds
    hbond_set = pose.get_hbonds()
    print(f"Rosetta found {hbond_set.nhbonds()} hydrogen bonds.")
    
    total_hb_energy = 0.0
    matched_edges = 0
    
    # We iterate Rosetta HBonds and match to Graph Edges
    # This handles the case where G has edge.
    # Note: G might track 'donated' H. Rosetta tracks Donor Atom and Acc Atom.
    # We need to map Rosetta (ResI, AtomName) to G Node ID.
    
    # Pre-calculate mapping from (Chain, ResNum, AtomName, Role) -> NodeID
    atom_to_node = {}
    for n, d in G.nodes(data=True):
        # We need to match based on how Rosetta names atoms vs PDB
        # Rosetta usually matches PDB atom names.
        # Key: (Chain, ResNum_Normalized, AtomName, Role)
        # Note: A single atom can be both Donor and Acceptor (e.g. water, His).
        # Rosetta Hbond has donor_atom and acceptor_atom.
        # donor_atom is the Heavy Atom!
        
        # d['resnum'] in G is normalized string (e.g. "341")
        key = (d['chain'], d['resnum'], d['atom'], d['role'])
        atom_to_node[key] = n

    pdb_info = pose.pdb_info()
    
    for i in range(1, hbond_set.nhbonds() + 1):
        hb = hbond_set.hbond(i)
        energy = hb.energy()
        total_hb_energy += energy
        
        # Identify Donor
        d_res = hb.don_res()
        donor_res_obj = pose.residue(d_res)
        d_h_idx = hb.don_hatm()
        d_heavy_idx = donor_res_obj.atom_base(d_h_idx)
        d_atom_name = donor_res_obj.atom_name(d_heavy_idx).strip()
        
        # Identify Acceptor
        a_res = hb.acc_res()
        acceptor_res_obj = pose.residue(a_res)
        a_atom_idx = hb.acc_atm()
        a_atom_name = acceptor_res_obj.atom_name(a_atom_idx).strip()
        
        # Get PDB Info
        d_pdb_chain = pdb_info.chain(d_res)
        d_pdb_res = pdb_info.number(d_res) # Integer
        # ICD code? pdb_info.icode(d_res)
        d_resnum = str(d_pdb_res) # normalization likely matches str(int)
        
        a_pdb_chain = pdb_info.chain(a_res)
        a_pdb_res = pdb_info.number(a_res)
        a_resnum = str(a_pdb_res)
        
        # Handle Waters: Rosetta might define chain differently if multiple waters allow?
        # Usually waters have chain and resnum.
        
        # Lookup in Graph
        # Donor Node: role='donor'
        d_key = (d_pdb_chain, d_resnum, d_atom_name, 'donor')
        # Acceptor Node: role='acceptor'
        a_key = (a_pdb_chain, a_resnum, a_atom_name, 'acceptor')
        
        u = atom_to_node.get(d_key)
        v = atom_to_node.get(a_key)
        
        if u and v:
            if G.has_edge(u, v):
                # Update energy
                G[u][v]['energy'] = energy
                G[u][v]['rosetta_found'] = True
                matched_edges += 1
            else:
                 # Rosetta found edge, but HBPLUS didn't?
                 # User cares about graph from HBPLUS. We can optionally add it?
                 # "Compute energies for each of the interactions IN OUR HYDROGEN BOND GRAPH"
                 # So we only annotate existing edges.
                 pass
        else:
             # Debug partial match?
             pass

    print(f"Total Rosetta H-bond Energy (Model): {total_hb_energy:.3f} REU")
    print(f"Matched {matched_edges} Rosetta interactions to Graph Edges.")

def write_pymol_script(G, filename="visualize_hbonds.pml"):
    print(f"Writing PyMOL script to {filename}...")
    with open(filename, 'w') as f:
        f.write("from pymol.cgo import *\n")
        f.write("from pymol import cmd\n\n")
        
        # Categories
        obj_aa_aa_bb_bb = []
        obj_aa_aa_sc_sc = []
        obj_aa_aa_bb_sc = []
        obj_aa_wat = []
        obj_wat_wat = []
        obj_weak = [] # New category for REU >= 0.0
        
        count = 0
        hbond_edges = 0
        
        count_red = 0
        
        for u, v, data in G.edges(data=True):
            if data.get('type') == 'hbond':
                hbond_edges += 1
                h_xyz = data.get('h_xyz')
                a_xyz = data.get('a_xyz')
                
                if h_xyz and a_xyz:
                    x1, y1, z1 = h_xyz
                    x2, y2, z2 = a_xyz
                    
                    # Contract 80%
                    dx, dy, dz = x2 - x1, y2 - y1, z2 - z1
                    dist = math.sqrt(dx*dx + dy*dy + dz*dz)
                    if dist > 0.001:
                        sx, sy, sz = x1 + dx*0.1, y1 + dy*0.1, z1 + dz*0.1
                        ex, ey, ez = x1 + dx*0.9, y1 + dy*0.9, z1 + dz*0.9
                    else:
                        sx, sy, sz, ex, ey, ez = x1, y1, z1, x2, y2, z2
                    
                    # Determine Color and Category
                    energy = data.get('energy', 0.0)
                    
                    if energy >= 0.0:
                        count_red += 1
                        col = "1.0, 0.0, 0.0" # Red
                        res_list = obj_weak
                    else:
                        hb_type = get_hbond_type(G.nodes[u], G.nodes[v])
                        
                        # Colors
                        # sc_sc: Yellow [1,1,0]
                        # sc_bb: Magenta [1,0,1]
                        # bb_bb: Green [0,1,0]
                        # aa_wat: Cyan [0,1,1]
                        # wat_wat: Density (Light Blue) [0.6, 0.6, 1.0]
                        
                        if hb_type == 'aa_aa_sc_sc':
                            col = "1.0, 1.0, 0.0"
                            res_list = obj_aa_aa_sc_sc
                        elif hb_type == 'aa_aa_bb_sc':
                            col = "1.0, 0.0, 1.0"
                            res_list = obj_aa_aa_bb_sc
                        elif hb_type == 'aa_aa_bb_bb':
                            col = "0.0, 1.0, 0.0"
                            res_list = obj_aa_aa_bb_bb
                        elif hb_type == 'aa_wat':
                            col = "0.0, 1.0, 1.0"
                            res_list = obj_aa_wat
                        else: # wat_wat
                            col = "0.6, 0.6, 1.0"
                            res_list = obj_wat_wat
                    
                    # CYLINDER
                    cyl_str = f"CYLINDER, {sx:.3f}, {sy:.3f}, {sz:.3f}, {ex:.3f}, {ey:.3f}, {ez:.3f}, 0.1, {col}, {col}"
                    res_list.append(cyl_str)
                    count += 1

        # Write lists
        def write_list(name, items):
            if not items: return
            f.write(f"{name} = [\n")
            for item in items:
                f.write(f"   {item},\n")
            f.write("]\n")
            f.write(f"cmd.load_cgo({name}, '{name}')\n\n")

        write_list("aa_aa_bb_bb", obj_aa_aa_bb_bb)
        write_list("aa_aa_sc_sc", obj_aa_aa_sc_sc)
        write_list("aa_aa_bb_sc", obj_aa_aa_bb_sc)
        write_list("aa_wat", obj_aa_wat)
        write_list("wat_wat", obj_wat_wat)
        write_list("hbond_weak", obj_weak)
        
        print(f"Wrote {count} H-bond cylinders to {filename}")
        print(f"Number of cylinders with non-favorable energy (REU >= 0): {count_red} (available object 'hbond_weak')")

import argparse
import subprocess
import os.path

def run_hbplus(pdb_file):
    # Determine expected output name
    # HBPLUS logic: usually replaces extension or appends?
    # Actually, hbplus output is based on PDB root name.
    # If file is '2ij2.opt.pdb', hbplus creates '2ij2.opt.hb2'?
    # Or '2ij2.hb2'?
    # It seems to handle dot correctly? 
    # Let's verify by just running it.
    # Assuming hbplus is in PATH.
    
    print(f"Running HBPLUS on {pdb_file}...")
    print(f"Running HBPLUS on {pdb_file}...")
    try:
        # Determine hbplus path relative to this script
        script_dir = os.path.dirname(os.path.abspath(__file__))
        local_hbplus = os.path.join(script_dir, "hbplus")
        
        if os.path.exists(local_hbplus):
            hbplus_cmd = local_hbplus
        else:
            # Fallback
            hbplus_cmd = "hbplus"
            
        subprocess.run([hbplus_cmd, pdb_file], check=True)
    except Exception as e:
        print(f"Error running hbplus ({hbplus_cmd}): {e}")
        # If hbplus failed, maybe we can proceed if .hb2 already exists?
        # Check below logic.
        pass # Let the next check see if file exists.
        
    # Predict hb2 filename
    root, ext = os.path.splitext(pdb_file)
    hb2_file = root + ".hb2"
    
    # HBPLUS sometimes truncates names or behaves weirdly with long names?
    # Ideally checking if file exists.
    if not os.path.exists(hb2_file):
        # Maybe it used '2ij2.hb2' if input was '2ij2.opt.pdb'?
        # HBPLUS convention is complex.
        # Simple check: 
        candidates = [root + ".hb2", root.split('.')[0] + ".hb2"]
        for c in candidates:
            if os.path.exists(c):
                return c
        print(f"Warning: Expected output {hb2_file} not found.")
        return None
    return hb2_file

def main():
    parser = argparse.ArgumentParser(description="Build H-bond Network Graph from PDB.")
    parser.add_argument("pdb_file", help="Input PDB file")
    
    args = parser.parse_args()
    pdb_file = args.pdb_file
    
    if not os.path.exists(pdb_file):
        print(f"Error: {pdb_file} not found.")
        sys.exit(1)
    
    # Run HBPLUS
    hb2_file = run_hbplus(pdb_file)
    if not hb2_file:
        print("Could not find HBPLUS output. Exiting.")
        sys.exit(1)
        
    print(f"Using HBPLUS output: {hb2_file}")

    print("Parsing PDB...")
    pdb_map = PDBMap(pdb_file)
    print(f"Loaded {len(pdb_map.atoms)} atoms and {len(pdb_map.residues)} residues.")
    
    print("Building Graph...")
    G = build_graph(pdb_map, hb2_file)
    print(f"Graph stats: {G.number_of_nodes()} nodes, {G.number_of_edges()} edges")
    
    print("Computing Energies...")
    # NOTE: PyRosetta init is called inside add_rosetta_energies. 
    # Since this is a script, safe to call.
    add_rosetta_energies(G, pdb_file)
    
    root, _ = os.path.splitext(pdb_file)
    output_pkl = f"{root}.graph.pkl"
    output_pml = f"{root}_hbonds.py"
    
    write_pymol_script(G, output_pml)
    
    import pickle
    with open(output_pkl, "wb") as f:
        pickle.dump(G, f)
    print(f"Saved graph to {output_pkl}")

if __name__ == "__main__":
    main()
