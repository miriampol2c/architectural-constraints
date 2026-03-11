import os
import shutil
import json
import glob
from pathlib import Path
import csv
import fnmatch
import pandas as pd
import argparse
from Bio.PDB import MMCIFParser, PDBIO
import re

# T must be given with decimals!! eg 1.0
def parse_args():
    parser = argparse.ArgumentParser(description="Select the best AF3 design based on th pLDDT scores.")
    parser.add_argument("--pdb_folder", type=str, required=True, help="Absolute path to the folder containing .pdb files") # Tracking the filtered folder
    return parser.parse_args()

# FUNCTIONS

# Function to retrieve all the output folders
def find_output_folders(root_path):
    # Construct the search pattern
    pattern = os.path.join(root_path, '*MPNN_output', 'seqs', 'AF3')
    # Use glob to find all matching directories
    matching_folders = glob.glob(pattern)
    return matching_folders

# Function to count samples folded (even if folding has not been completed, the directory is created)
def count_final_folders(output_folder):
    # List all subdirectories in the given output folder
    subdirs = [d for d in os.listdir(output_folder) if os.path.isdir(os.path.join(output_folder, d))]
    return len(subdirs)

# Function to search for a file
def find_file(root_path, file_name):
    for root, dirs, files in os.walk(root_path):
        if file_name in files:
            return os.path.join(root, file_name)
    return None

# Function to search for several files
def find_files(directory, pattern):
    files_found = []
    for root, dirs, files in os.walk(directory):
        for file in files:
            if fnmatch.fnmatch(file, pattern):
                files_found.append(os.path.join(root, file))
    print(f"Files found: {files_found}")
    return files_found

# Function to copy files
def copy_file_to_folder(src_file, dest_folder):
    dest_folder = Path(dest_folder)
    if not dest_folder.exists():
        dest_folder.mkdir(parents=True)
    shutil.copy(src_file, dest_folder)

#  AF3 Function to merge csv files
def merge_csv_files_af3(input_dir, output_file):
    input_dir = Path(input_dir)
    csv_files = sorted(input_dir.glob("*_AF3_results.csv"))

    with open(output_file, "w", newline='') as fout:
        writer = None
        for i, csv_file in enumerate(csv_files):
            with open(csv_file, "r", newline='') as fin:
                reader = csv.reader(fin)
                header = next(reader)
                if writer is None:
                    writer = csv.writer(fout)
                    writer.writerow(header)  # Write header only once
                for row in reader:
                    writer.writerow(row)

# AF3 Function to select best design per target given a full csv file data
def select_design_af3(best_designs_path):
    best_designs_path = Path(best_designs_path)
    merged_csv = best_designs_path / "merged_best_plddt_af3.csv"
    output_csv = best_designs_path / "best_designs_selected.csv"

    # Check if the merged CSV exists and is not empty
    if not merged_csv.exists():
        print(f"Error: File does not exist: {merged_csv}")
        return
    if merged_csv.stat().st_size == 0:
        print(f"Error: File is empty: {merged_csv}")
        return

    # Load the CSV
    df = pd.read_csv(merged_csv)
    if df.empty:
        print(f"Error: CSV has no rows: {merged_csv}")
        return

    # Extract Target and Temperature
    # check t in lower case because is set by af3
    #df["Temperature"] = df["Sample_fullname"].str.extract(r'_t(\d{2})_')
    #df["Sample"] = df["Sample_fullname"].str.extract(r'_sample(\d+)$')
    # Extract temperature, allowing decimals (e.g., 0.1, 25.0, etc.)
    df["Temperature"] = df["Sample_fullname"].str.extract(r'_t([0-9.]+)_')
    # Extract sample number at the end (e.g., sample1, sample12
    df["Sample"] = df["Sample_fullname"].str.extract(r'_sample(\d+)$')

    # Select best design per (Target, Temperature)
    best_per_target_temp = df.loc[
        df.groupby(["Target", "Temperature"])["AF3_score"].idxmax()
    ].reset_index(drop=True)

    # Save filtered CSV
    best_per_target_temp.to_csv(output_csv, index=False)
    print(f"Saved best designs to: {output_csv}")

    # Directories
    per_target_dir = best_designs_path / "per_target"
    per_sample_dir = best_designs_path / "per_target_and_sample"
    per_sample_dir.mkdir(parents=True, exist_ok=True)

    # Extract sequence names
    seqs_names = []
    with open(output_csv, 'r', encoding='utf-8') as fi:
        reader_csv = csv.reader(fi)
        next(reader_csv)  # skip header
        for row in reader_csv:
            #if row:
            #    full_seq_name = row[0].strip()
            #    seq_name = full_seq_name.split(':')[0].strip()
            #    seqs_names.append(seq_name)
            seq_name = row[1]
            seqs_names.append(seq_name)

    print("Extracted sequence names:")
    for seq_name in seqs_names:
        print(f"- {seq_name}")

    print(f"Total sequences: {len(seqs_names)}")

    # Copy selected PDBs to final dir (per sample and target)
    for seq_name in seqs_names:
        source_pdb = per_target_dir / f"{seq_name}.pdb"
        dest_pdb = per_sample_dir / f"{seq_name}.pdb"

        if source_pdb.exists():
            shutil.copy2(source_pdb, dest_pdb)
            print(f"Copied: {source_pdb.name} → {per_sample_dir}")
        else:
            print(f"Missing: {source_pdb}")


def convert_cif_to_pdb(cif_file_path, output_dir):
    parser = MMCIFParser(QUIET=True)
    io = PDBIO()

    structure_id = os.path.basename(cif_file_path).replace('_model', '').replace ('.cif', '')
    structure = parser.get_structure(structure_id, cif_file_path)

    if output_dir is None:
        output_dir = os.path.dirname(cif_file_path)

    pdb_file_path = os.path.join(output_dir, f"{structure_id}.pdb")
    io.set_structure(structure)
    io.save(pdb_file_path)

    return pdb_file_path

def extract_target_name(full_path):
    # Traverse backwards through parts of the path to find the MPNN folder
    for part in reversed(Path(full_path).parts):
        match = re.match(r'^(.+?)_MPNN_output$', part)
        if match:
            return match.group(1)
    # Fallback: nothing matched
    return None

# MAIN
def main():
    args = parse_args()
    root_path = Path(args.pdb_folder)
    best_designs_path = root_path / "Best_designs_AF3"
    best_designs_path.mkdir(parents=True, exist_ok=True)

    # Retrieve output folders with AF3 predictions
    print("Checking targets with designed and folded sequences")
    output_folders = find_output_folders(root_path)

    # Count the number of output folders with AF2 preidctions
    number_of_output_folders = len(output_folders)
    print(f"Number of targets with AF3 folder: {number_of_output_folders}")

    # Iterate over each output folder and count its subdirectories (samples designed and totally or partially folded)
    for folder in output_folders:
        final_folder_count = count_final_folders(folder)
        
        if final_folder_count == 0:
            print(f"\nFolder '{folder}' contains {final_folder_count} subdirectories (designed samples)")
            continue
            
        if final_folder_count > 50:
            print(f"\nFolder '{folder}' contains {final_folder_count} subdirectories (designed samples)")
            print(f"50 samples per target is the maximum. Check the target folder, for now it has been discarded.")
            continue

        print(f"Folder '{folder}' contains {final_folder_count} subdirectories (designed samples)")
        
        # Attempt to extract target name calling a specific function
        target_name = extract_target_name(folder)

        if target_name:
            print(f"--> TARGET NAME: {target_name}")
        else:
            print("MPNN_output folder not found in path.")

        # Iterate through designed samples
        for subdir in os.listdir(folder):
            full_path = os.path.join(folder, subdir)
            if os.path.isdir(full_path):
                subdir_name = Path(full_path).name
                print(f"{'-' * 50}")
                print(f"Sample: {subdir_name}")
                print(f"{'-' * 50}")
                #print(f"Processing {target_name}, {subdir_name}")
                #print(f"{'-' * 50}")
                # AF3 already outputs a ransking scores, and the best strcuture according to it its already in the af3 root path
                output_csv = best_designs_path / f"{subdir_name}_AF3_results.csv" 
    
                # Step 1: Always create (or overwrite) the CSV file with header
                with open(output_csv, mode='w', newline='') as file:
                    writer = csv.writer(file)
                    writer.writerow(["Target", "Sample_fullname", "AF3_Model", "AF3_score"])

                # Step 2: Now process the ranking_debug.json and append to this CSV
                ranking_file = Path(full_path) / "ranking_scores.csv"

                if ranking_file.exists():
                    print(f"Ranking file found for {subdir_name}")
                    scores = {} 
                    with open(ranking_file, 'r') as file:
                        reader = csv.DictReader(file)
                        for row in reader:
                            n_model = int(row['sample'])
                            score = float(row['ranking_score'])
                            scores[n_model] = score
                        if scores:
                            with open(output_csv, mode='a', newline='') as file:
                                writer = csv.writer(file)
                                for model, score in scores.items():
                                    writer.writerow([target_name, subdir_name, model, score])

                            max_model = max(scores, key=scores.get)
                            max_value = scores[max_model]
                            print(f"Model with highest AF3 score: {max_model}, Value: {max_value}")

                            cif_file = find_file(full_path, f"{subdir_name}_model.cif")
                            if cif_file:
                                print(".cif  file found:", cif_file)
                                print(".cif file being converted to PDB")
                                out_path = best_designs_path / "per_target"
                                out_path.mkdir(parents=True, exist_ok=True)
                                convert_cif_to_pdb(cif_file, out_path) # directamente se guarda en out_path
                                #because The prediction with highest ranking is the one included in the root directory.
                                #est_file = out_path / f"{subdir_name}.pdb"
                                #shutil.copy(pdb_file, dest_file)
                                print(f"Converted and saved in: {out_path}")
                            else:
                                print("PDB file not found:", cif_file)
                        else:
                            print("No pLDDT data found in the JSON file")
                else:
                    print(f"\nRanking file not found for {subdir_name}")

    # Merge CSVs and select best design per target
    merge_csv_files_af3(best_designs_path, best_designs_path / "merged_best_plddt_af3.csv")
    select_design_af3(best_designs_path)
 
if __name__ == "__main__":
    main()
