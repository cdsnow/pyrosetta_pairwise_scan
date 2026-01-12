import json
import os
import sys

def check_sterics(results_file):
    if not os.path.exists(results_file):
        print(f"Error: {results_file} not found.")
        return

    print(f"Loading {results_file}...")
    with open(results_file, 'r') as f:
        data = json.load(f)
    
    results = data.get('results', {})
    total_models = 0
    clashing_models = 0
    
    for entry in results.values():
        for model in entry.get('models', []):
            if 'error' in model:
                continue
            
            total_models += 1
            # Default to a high value if missing to count as a clash
            fa_rep = model.get('delta_fa_rep', 999.0)
            
            if fa_rep > 10.0:
                clashing_models += 1
                
    if total_models == 0:
        print("No models found in results.")
        return

    clash_pct = (clashing_models / total_models) * 100
    safe_models = total_models - clashing_models
    safe_pct = (safe_models / total_models) * 100

    print("\nSteric Clash Statistics (delta_fa_rep > 10.0):")
    print(f"Total models evaluated:  {total_models:,}")
    print(f"Clashing models:         {clashing_models:,} ({clash_pct:.1f}%)")
    print(f"Safe models:             {safe_models:,} ({safe_pct:.1f}%)")

if __name__ == "__main__":
    filename = sys.argv[1] if len(sys.argv) > 1 else 'merged_results.json'
    check_sterics(filename)
