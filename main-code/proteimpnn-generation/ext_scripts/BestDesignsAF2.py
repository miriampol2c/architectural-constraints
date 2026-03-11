import os
import shutil
import json
import glob
from pathlib import Path
import csv
import fnmatch
import pandas as pd
import argparse
import re

def parse_args():
    parser = argparse.ArgumentParser(description="Select the best AF2 design based on th pLDDT scores.")
    parser.add_argument("--pdb_folder", type=str, required=True, help="Absolute path to the folder containing .pdb files") # Tracking the filtered folder
    return parser.parse_args()

# FUNCTIONS

# Function to retrieve all the output folders
def find_output_folders(root_path):
    # Construct the search pattern
    pattern = os.path.join(root_path, '*MPNN_output', 'seqs', 'AF2')
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

#  AF2 Function to merge csv files
def merge_csv_files_af2(input_dir, output_file):
    input_dir = Path(input_dir)
    csv_files = sorted(input_dir.glob("*_AF2_results.csv"))

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

# AF2 Function to select best design per target given a full csv file data
def select_design_af2(best_designs_path):
    best_designs_path = Path(best_designs_path)
    merged_csv = best_designs_path / "merged_best_plddt_af2.csv"
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
    #df["Target_ID"] = df["Target"].str.extract(r'([^_]+)')  # before underscore
    #df["Temperature"] = df["Target"].str.extract(r'T=([\d.]+)').astype(float)
    # now names do not contain . or = (change for a complete new run)
    #df["Temperature"] = df["Sample"].str.extract(r'T([\d]+)').astype(float)
    #df["Temperature"] = df["Sample_fullname"].str.extract(r'_T(\d{2})_')
    #df["Sample"] = df["Sample_fullname"].str.extract(r'_sample(\d+)$')
    # Extract temperature, allowing decimals (e.g., 0.1, 25.0, etc.)
    df["Temperature"] = df["Sample_fullname"].str.extract(r'_T([0-9.]+)_')
    # Extract sample number at the end (e.g., sample1, sample12
    df["Sample"] = df["Sample_fullname"].str.extract(r'_sample(\d+)$')


    # Select best design per (Target, Temperature)
    best_per_target_temp = df.loc[
        df.groupby(["Target", "Temperature"])["pLDDT_value"].idxmax()
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

    # Copy selected PDBs
    for seq_name in seqs_names:
        source_pdb = per_target_dir / f"{seq_name}.pdb"
        dest_pdb = per_sample_dir / f"{seq_name}.pdb"

        if source_pdb.exists():
            shutil.copy2(source_pdb, dest_pdb)
            print(f"Copied: {source_pdb.name} → {per_sample_dir}")
        else:
            print(f"Missing: {source_pdb}")

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
    best_designs_path = root_path / "Best_designs_AF2"
    #per_target_path = best_designs_path / "per_target"
    #per_sample_path = best_designs_path / "per_target_and_sample"
    best_designs_path.mkdir(parents=True, exist_ok=True)
    #per_target_path.mkdir(parents=True, exist_ok=True)
    #per_sample_path.mkdir(parents=True, exist_ok=True)

    # Retrieve output folders with AF2 predictions
    print("Checking targets with designed and folded sequences")
    output_folders = find_output_folders(root_path)

    # Count the number of output folders with AF2 preidctions
    number_of_output_folders = len(output_folders)
    print(f"Number of targets with AF2 folder: {number_of_output_folders}")

    # Iterate over each output folder and count its subdirectories (samples designed and totally or partially folded)
    for folder in output_folders:
        final_folder_count = count_final_folders(folder)
        if final_folder_count > 50:
            print(f"\nFolder '{folder}' contains {final_folder_count} subdirectories (designed samples)")
            print(f"50 samples per target is the maximum. Check the target folder, for now it has been discarded.")
            continue

        print(f"\nFolder '{folder}' contains {final_folder_count} subdirectories (designed samples)")
        
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
                output_csv = best_designs_path / f"{subdir_name}_AF2_results.csv" 
    
                # Step 1: Always create (or overwrite) the CSV file with header
                with open(output_csv, mode='w', newline='') as file:
                    writer = csv.writer(file)
                    writer.writerow(["Target", "Sample_fullname", "AF2_Model", "pLDDT_value"])

                # Step 2: Now process the ranking_debug.json and append to this CSV
                ranking_file = Path(full_path) / "ranking_debug.json"
                if ranking_file.exists():
                    print(f"Ranking file found for {subdir_name}")
                    with open(ranking_file, 'r') as file:
                        data = json.load(file)
                        plddts = data.get("plddts", {})
                        if plddts:
                            with open(output_csv, mode='a', newline='') as file:
                                writer = csv.writer(file)
                                for model, plddt_value in plddts.items():
                                    writer.writerow([target_name, subdir_name, model, plddt_value])

                            max_model = max(plddts, key=plddts.get)
                            max_value = plddts[max_model]
                            print(f"Model with highest pLDDT value: {max_model}, Value: {max_value}")

                            pdb_file = find_file(full_path, f"relaxed_{max_model}.pdb")
                            if pdb_file:
                                print("PDB file found:", pdb_file)
                                out_path = best_designs_path / "per_target"
                                out_path.mkdir(parents=True, exist_ok=True)
                                dest_file = out_path / f"{subdir_name}.pdb"
                                shutil.copy(pdb_file, dest_file)
                                print(f"Copied to: {dest_file}")
                            else:
                                print("PDB file not found:", pdb_file)
                        else:
                            print("No pLDDT data found in the JSON file")
                else:
                    print(f"\nRanking file not found for {subdir_name}")
        
    # Merge CSVs and select best design per target
    merge_csv_files_af2(best_designs_path, best_designs_path / "merged_best_plddt_af2.csv")
    select_design_af2(best_designs_path)
 
if __name__ == "__main__":
    main()