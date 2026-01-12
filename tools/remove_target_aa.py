import os
import glob

def remove_chain_f(file_path):
    with open(file_path, 'r') as f:
        lines = f.readlines()
    
    new_lines = []
    for line in lines:
        if line.startswith('ATOM') or line.startswith('HETATM'):
            # Chain ID is at index 21
            if len(line) > 21 and line[21] == 'F':
                continue
        new_lines.append(line)
    
    with open(file_path, 'w') as f:
        f.writelines(new_lines)

def main():
    base_dir = 'candidates_v2'
    pdb_files = glob.glob(os.path.join(base_dir, '**', '*.pdb'), recursive=True)
    
    print(f"Processing {len(pdb_files)} files...")
    count = 0
    for pdb_file in pdb_files:
        remove_chain_f(pdb_file)
        count += 1
        if count % 100 == 0:
            print(f"Processed {count} files")
    
    print("Done.")

if __name__ == "__main__":
    main()
