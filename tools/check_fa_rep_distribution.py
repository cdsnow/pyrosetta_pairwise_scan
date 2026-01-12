import json
import sys

def main():
    with open('merged_results.json') as f:
        data = json.load(f)
    
    results = data.get('results', {})
    
    total_models = 0
    fa_rep_counts = {10: 0, 20: 0, 50: 0, 100: 0, 1000: 0, 10000: 0, 100000: 0}
    
    for entry in results.values():
        for m in entry.get('models', []):
            if 'error' in m: continue
            total_models += 1
            fa = m.get('delta_fa_rep', 999)
            
            for thresh in fa_rep_counts:
                if fa <= thresh:
                    fa_rep_counts[thresh] += 1
                    
    print(f"Total models: {total_models}")
    for thresh in sorted(fa_rep_counts.keys()):
        count = fa_rep_counts[thresh]
        pct = (count / total_models) * 100 if total_models else 0
        print(f"delta_fa_rep <= {thresh}: {count} ({pct:.1f}%)")

if __name__ == "__main__":
    main()
