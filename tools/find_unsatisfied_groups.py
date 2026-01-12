import argparse
import networkx as nx
import pickle
import json
import sys

def main():
    parser = argparse.ArgumentParser(description="Find unsatisfied H-bond donors and acceptors in graph.")
    parser.add_argument("graph_file", help="Path to .graph.pkl file")
    parser.add_argument("--output", default="unsatisfied.json", help="Output JSON file")
    parser.add_argument("--pml", help="Output PyMOL script (.pml)")
    parser.add_argument("--chain", help="Filter by Chain ID (e.g. 'F')")
    args = parser.parse_args()

    # Load graph
    try:
        with open(args.graph_file, "rb") as f:
            G = pickle.load(f)
    except Exception as e:
        print(f"Error loading graph: {e}")
        sys.exit(1)
        
    unsatisfied = {
        "donors": [],
        "acceptors": []
    }
    
    # Iterate nodes
    for n, data in G.nodes(data=True):
        role = data.get('role')
        if not role: continue
        
        # Filter by Chain if requested
        if args.chain and data.get('chain') != args.chain:
            continue
            
        # Omit waters from unsatisfied reporting
        if data.get('resname') == 'HOH':
            continue
        
        # Check edges
        has_hbond = False
        for neighbor in G.neighbors(n):
            edge_data = G[n][neighbor]
            if edge_data.get('type') == 'hbond':
                has_hbond = True
                break
        
        if not has_hbond:
            # Add to list
            item = {
                "chain": data.get('chain'),
                "resnum": data.get('resnum'),
                "resname": data.get('resname'),
                "atom": data.get('atom'),
                "type": data.get('type') # BB_Donor etc
            }
            if role == 'donor':
                unsatisfied["donors"].append(item)
            elif role == 'acceptor':
                unsatisfied["acceptors"].append(item)
                
    # Save
    with open(args.output, "w") as f:
        json.dump(unsatisfied, f, indent=2)
    
    print(f"Found {len(unsatisfied['donors'])} unsatisfied donors and {len(unsatisfied['acceptors'])} unsatisfied acceptors.")
    print(f"Saved to {args.output}")

    if args.pml:
        write_pymol_script(unsatisfied, args.pml)
        print(f"Saved PyMOL script to {args.pml}")

def write_pymol_script(unsatisfied, filename):
    with open(filename, 'w') as f:
        # Donors
        selections = []
        for d in unsatisfied['donors']:
            # Construct atom selection. 
            # Note: resnum in graph is string. PyMOL handles 'resi 100' or 'resi 100A' fine.
            sel = f"(chain {d['chain']} and resi {d['resnum']} and name {d['atom']})"
            selections.append(sel)
        
        if selections:
            # Join with ' or '
            # Chunking to avoid extremely long lines if needed, but PyMOL usually OK with distinct or blocks.
            full_sel = " or ".join(selections)
            f.write(f"select unsat_donors, {full_sel}\n")
        else:
            f.write("select unsat_donors, none\n")
            
        # Acceptors
        selections = []
        for a in unsatisfied['acceptors']:
            sel = f"(chain {a['chain']} and resi {a['resnum']} and name {a['atom']})"
            selections.append(sel)
            
        if selections:
            full_sel = " or ".join(selections)
            f.write(f"select unsat_acceptors, {full_sel}\n")
        else:
            f.write("select unsat_acceptors, none\n")
            
        f.write("show spheres, unsat_donors\n")
        f.write("show spheres, unsat_acceptors\n")
        f.write("set sphere_scale, 0.4\n")
        # Visual styling
        f.write("color white, unsat_donors\n")
        f.write("color pink, unsat_acceptors\n")
        f.write("deselect\n")

if __name__ == "__main__":
    main()
