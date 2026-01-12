import json
import os
import random

def select_diverse_candidates(json_path, n=10):
    with open(json_path, 'r') as f:
        data = json.load(f)
    
    results = data.get('results', {})
    print(f"Total entries in results: {len(results)}")
    
    # Group by category (parent folder)
    categories = {}
    first = True
    for pdb_path, entry in results.items():
        if not entry:
            print(f"Empty entry for {pdb_path}")
            continue
            
        if 'error' in entry:
            # print(f"Error in entry {pdb_path}: {entry['error']}")
            continue
            
        if not entry.get('models'):
            # print(f"No models for {pdb_path}")
            continue

        if first:
            print(f"Sample model: {entry['models'][0]}")
            first = False
            
        # Extract category (e.g., candidates_v2/F1_ASP_BB_N/...)
        parts = pdb_path.split('/')
        if len(parts) >= 3:
            category = parts[1] # F1_ASP_BB_N
        else:
            print(f"Path too short: {pdb_path}")
            category = 'misc'
            
        if category not in categories:
            categories[category] = []
            
        # Find best model for this PDB
        best_model = None
        best_score = 9999
        models = entry.get('models', [])
        if not models:
             # print(f"No models list for {pdb_path}")
             continue
             
        for model in models:
            if 'error' in model: continue
            if model.get('delta_total', 9999) < best_score:
                best_score = model['delta_total']
                best_model = model
        
        if best_model:
            categories[category].append({
                'pdb_path': pdb_path,
                'model': best_model,
                'score': best_score
            })
        else:
            pass
            # print(f"No valid models for {pdb_path}")
            
    # Select best from each category
    selected = []
    category_names = sorted(categories.keys())
    
    print(f"Found {len(category_names)} categories: {category_names}")
    
    # Round 1: Best from each category
    for cat in category_names:
        # Sort candidates in this category by score
        categories[cat].sort(key=lambda x: x['score'])
        if categories[cat]:
            selected.append(categories[cat][0])
            
    # Round 2: If we need more, pick next bests
    idx = 1
    while len(selected) < n:
        added_any = False
        for cat in category_names:
            if len(selected) >= n: break
            if idx < len(categories[cat]):
                selected.append(categories[cat][idx])
                added_any = True
        idx += 1
        if not added_any:
            break
            
    return selected[:n]

selected = select_diverse_candidates('merged_results.json', 10)

print(f"Selected {len(selected)} candidates:")
for item in selected:
    print(f"{item['pdb_path']} (Model {item['model']['model_num']}): {item['score']:.2f}")

# Save list for the next script
with open('selected_candidates.json', 'w') as f:
    json.dump(selected, f, indent=2)
