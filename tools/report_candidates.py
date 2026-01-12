#!/usr/bin/env python
import os
import statistics
import sys

def count_models(pdb_path):
    count = 0
    try:
        with open(pdb_path, 'r') as f:
            for line in f:
                if line.startswith('MODEL'):
                    count += 1
    except Exception as e:
        return 0
    return count

def main():
    try:
        base_dir = sys.argv[1]
    except:
        print('Provide the directory to report on')
        sys.exit()

    if not os.path.isdir(base_dir):
        print(f"Directory '{base_dir}' not found.")
        return

    query_targets = sorted([d for d in os.listdir(base_dir) if os.path.isdir(os.path.join(base_dir, d))])
    
    print(f"Total number of query/targets: {len(query_targets)}")
    print("=" * 80)
    
    total_motifs = 0
    for qt in query_targets:
        qt_path = os.path.join(base_dir, qt)
        signatures = [f for f in os.listdir(qt_path) if f.endswith('.pdb')]
        
        print(f"Query/Target: {qt}")
        print(f"  Number of neighbor signatures: {len(signatures)}")
        
        if not signatures:
            print("    (No signatures found)")
            print("-" * 80)
            continue
            
        motif_counts = []
        for sig in signatures:
            sig_path = os.path.join(qt_path, sig)
            num_motifs = count_models(sig_path)
            sig_name = sig.replace('.pdb', '')
            motif_counts.append((sig_name, num_motifs))
            total_motifs += num_motifs
        
        # Extract just the counts for statistics
        counts = [c for _, c in motif_counts]
        
        if not counts:
             print("    (No models found in signatures)")
             print("-" * 80)
             continue

        min_motifs = min(counts)
        mean_motifs = statistics.mean(counts)
        median_motifs = statistics.median(counts)
        
        # Sort by count descending for top 10
        motif_counts.sort(key=lambda x: x[1], reverse=True)
        top_10 = motif_counts[:10]
        
        print(f"  Statistics (Motifs per Signature):")
        print(f"    Min:    {min_motifs}")
        print(f"    Mean:   {mean_motifs:.2f}")
        print(f"    Median: {median_motifs:.2f}")
        print(f"  Top 10 Signatures by Motif Count:")
        for name, count in top_10:
             print(f"    - {name:<40}: {count}")
             
        print("-" * 80)
    
    print(f"\nTOTAL MOTIFS ACROSS ALL TARGETS: {total_motifs}")

if __name__ == "__main__":
    main()
