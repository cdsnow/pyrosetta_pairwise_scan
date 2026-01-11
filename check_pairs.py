import json

with open('pruned_results.json') as f:
    data = json.load(f)

entry = data['results'].get('candidates_v3/F6_ASP_SC/ARG_NE.pdb')
if not entry:
    print("Entry not found")
else:
    models = {m['model_num']: m for m in entry['models']}
    
    pairs = [(8, 9), (12, 14)]
    
    for m1_id, m2_id in pairs:
        m1 = models.get(m1_id)
        m2 = models.get(m2_id)
        
        print(f"--- Pair {m1_id} vs {m2_id} ---")
        if not m1 or not m2:
            print("Model missing")
            continue
            
        print(f"Model {m1_id}: Total={m1['delta_total']:.3f}, HBond={m1['delta_hbond']:.3f}")
        print(f"Model {m2_id}: Total={m2['delta_total']:.3f}, HBond={m2['delta_hbond']:.3f}")
        
        # Check domination
        m1_dom_m2 = (m1['delta_total'] <= m2['delta_total'] and m1['delta_hbond'] <= m2['delta_hbond'])
        m2_dom_m1 = (m2['delta_total'] <= m1['delta_total'] and m2['delta_hbond'] <= m1['delta_hbond'])
        
        print(f"{m1_id} dominates {m2_id}? {m1_dom_m2}")
        print(f"{m2_id} dominates {m1_id}? {m2_dom_m1}")
        print("")
