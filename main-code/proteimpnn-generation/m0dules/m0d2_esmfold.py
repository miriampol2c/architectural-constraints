import os
import re
import glob
import subprocess
import argparse
import sys
import re
import json
import torch
import esm
import time
import biotite.structure.io as bsio # type: ignore
from Bio import SeqIO
from Bio.PDB import PDBParser, PDBIO, Select
from collections import defaultdict
import shutil 
from pathlib import Path
import csv
from Bio import AlignIO
from Bio.Align.Applications import MafftCommandline
from io import StringIO
from Bio.Align import MultipleSeqAlignment
from Bio.SeqIO.FastaIO import FastaWriter
from pathlib import Path
from Bio import PDB
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.PDB import PDBParser
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from pathlib import Path

def parse_args():
    parser = argparse.ArgumentParser(description="Run FrustMPNN pipeline on a folder of PDB files.")
    parser.add_argument("--pdb_folder", type=str, required=True, help="Absolute path to the folder containing .pdb files")
    parser.add_argument("--modelling", type=str, required=True, default="ESMFold", help="Structure prediction method ESMFold")

    return parser.parse_args()

###########################################
# FUNCTIONS TO PREDICT ESMFOLD STRUCTURES #
###########################################

os.environ["TORCH_HOME"] = os.path.expanduser("/gpfs/projects/bsc08/bsc008609/esm/")

def load_ESMFold():
    # Empty the cache before loading the model
    torch.cuda.empty_cache()
    print("Empty the cache")
    
    # Load model
    model = esm.pretrained.esmfold_v1()
   
    model = model.eval().cuda()
    print("Model loaded")
    #print("Model evaluation")

    # Optionally, uncomment to set a chunk size for axial attention. This can help reduce memory.
    # Lower sizes will have lower memory requirements at the cost of increased speed.
    # model.set_chunk_size(128)

    return model
    
def run_ESMFold(fasta_file, model):

    # Input file path
    input_file_path = os.path.abspath(fasta_file)
    # Directory where the input file is located
    directory_path = os.path.dirname(input_file_path)
    # Extract base filename (without extension)
    base_filename, _ = os.path.splitext(os.path.basename(input_file_path))
    # Define the output directory name
    output_dir_name = f"{base_filename}_ESMFold"
    # Full path to the output directory
    output_dir = os.path.join(directory_path, output_dir_name)

    # Create the directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)

    # Load sequences from multi-FASTA file
    sequences = []
    sequences_counter = 0
    with open(fasta_file, "r") as f:
        for record in SeqIO.parse(f, "fasta"):
            sequences.append((record.id, str(record.seq)))  # Store (sequence_id, sequence) tuple
            sequences_counter += 1  # Increment sequence counter

    # Predict structures for each sequence
    plddts = []
    for i, (sequence_id, sequence) in enumerate(sequences):
        with torch.no_grad():
            output = model.infer_pdb(sequence)
            
        if sequence_id.count("_") == 7 and ("relaxed" in sequence_id or "unrelaxed" in sequence_id):
            short_sequence_id = sequence_id
        elif '-' in sequence_id and '_' in sequence_id:
            short_sequence_id = '_'.join(sequence_id.split("_")[:5])
        elif sequence_id.count("_") >= 4:
            short_sequence_id = "_".join(sequence_id.split("_")[:4])
        else:
            # If sequence_id contains two or fewer "_", short_sequence_id is the entire sequence_id
            short_sequence_id = sequence_id                
        
        #output_filename = f"{output_dir}/{short_sequence_id}.pdb"
        output_dir = Path(output_dir)  # ensure it's a Path
        output_filename = output_dir / f"{short_sequence_id}.pdb"
        with open(output_filename, "w") as f:
            f.write(output)
            torch.cuda.empty_cache() 

        # Load structure and calculate pLDDT
        struct = bsio.load_structure(output_filename, extra_fields=["b_factor"])
        plddts.append((short_sequence_id, struct.b_factor.mean()))
        torch.cuda.empty_cache() 

    if plddts:
        # Write pLDDT values to output file
        #prefix = base_filename.split("_")[0] if "_" in base_filename else base_filename.split("-")[0]
        prefix = base_filename
        with open(f"{output_dir}/{prefix}_plddts", "w") as f:
            for short_sequence_id, plddt in plddts:
                f.write(f"{short_sequence_id}: {plddt}\n")

        # Write sorted pLDDT values to output file
        with open(f"{output_dir}/{prefix}_plddts_sorted", "w") as f:
            for short_sequence_id, plddt in sorted(plddts, key=lambda x: x[1], reverse=True):
                f.write(f"{short_sequence_id}: {plddt}\n")
                
        # Write best design in an individual file if possible
        best_design = sorted(plddts, key=lambda x: x[1], reverse=True)[0]
        if best_design:
            with open(f"{output_dir}/{prefix}_best_design", "w") as f:
                short_sequence_id, plddt = best_design
                f.write(f"{short_sequence_id}: {plddt}\n")
    else:
        print(f"Warning: No designs found for {base_filename}, skipping best_design file.")

def extract_matching_sequences(input_file, mfa_file, output_file, line_width=60):
    """
    Extract matching sequences from a multi-FASTA file based on a list of names.

    Parameters:
        input_file (str or Path): File containing the selected design names (e.g., best_design).
        mfa_file (str or Path): Path to the multi-FASTA file.
        output_file (str or Path): File to save matching sequences.
        line_width (int): Number of characters per FASTA line (default 60).
    """
    input_file = Path(input_file)
    mfa_file = Path(mfa_file)
    output_file = Path(output_file)

    if not input_file.exists() or not mfa_file.exists():
        print(f"Error: One or more input files do not exist.")
        return

    # Read selected names
    with input_file.open('r', encoding='utf-8') as fi:
        names_list = [line.split(":")[0].strip() for line in fi if line.strip()]

    # Read multi-FASTA content
    with mfa_file.open("r", encoding='utf-8') as fa:
        content = fa.read()
        sequences = content.split(">")[1:]  # skip leading empty string

    found = 0
    with output_file.open("a", encoding='utf-8') as fo:
        for seq in sequences:
            header, *body = seq.split("\n", 1)
            body = "".join(body).replace("\n", "")
            if header in names_list:
                #print(f"Found: {header}")
                fo.write(f">{header}\n")
                # Wrap sequence lines
                #chunks = [body[i:i+line_width] for i in range(0, len(body), line_width)]
                #fo.write("\n".join(chunks) + "\n")
                # Write full sequence in a line
                fo.write(body + "\n")  # Write the entire sequence in a single line
                found += 1

    #print(f"Done. {found} matching sequences written to: {output_file}")


def copy_best_design_pdbs(input_csv, source_dir, destination_dir, verbose=True):
    """
    Copy PDB files corresponding to the best designs listed in a CSV file.

    Parameters:
        input_csv (str or Path): Path to the CSV file with sequence names.
        source_dir (str or Path): Directory containing all PDB files.
        destination_dir (str or Path): Directory to copy the selected PDBs to.
        verbose (bool): Whether to print detailed logs (default: True).
    """
    input_csv = Path(input_csv)
    source_dir = Path(source_dir)
    destination_dir = Path(destination_dir)

    if not input_csv.exists():
        print(f"Error: Input CSV file not found: {input_csv}")
        return
    if not source_dir.exists():
        print(f"Error: Source directory not found: {source_dir}")
        return

    # Create destination directory if it doesn't exist
    destination_dir.mkdir(parents=True, exist_ok=True)

    # Extract sequence names from the CSV
    sequence_names = []
    with input_csv.open("r", encoding="utf-8") as f:
        reader = csv.reader(f)
        for row in reader:
            if row:
                name = row[0].strip().split(":")[0]
                sequence_names.append(name)

    #if verbose: # just for debug
        #print(f"Total sequences to copy: {len(sequence_names)}")

    # Copy files
    copied = 0
    missing = 0
    for name in sequence_names:
        source_file = source_dir / f"{name}.pdb"
        dest_file = destination_dir / f"{name}.pdb"
        if source_file.exists():
            shutil.copy2(source_file, dest_file)
            copied += 1
            #if verbose: # just for debug
                #print(f"Copied: {name}.pdb")
        else:
            missing += 1
            if verbose:
                print(f"Missing: {name}.pdb in {source_dir}")

    #print(f"Done. {copied} PDB files copied, {missing} missing.")

# MAIN

def main(pdb_folder, modelling):
    
    new_fasta_files = sorted(glob.glob(os.path.join(pdb_folder, "*_output", "seqs", "*.fasta")))
    #print("new_fasta_files")
    #print(new_fasta_files)
    # ESMFold modelling
    if modelling == "ESMFold":

        # Load the model calling the function
        print("\n*** ESMFold structure prediction ***")
        model = load_ESMFold()
        # Run the model for each fasta
        for fasta_file in new_fasta_files:
            run_ESMFold(fasta_file, model)
            # output from here should not be a dependency list, as esmfold is faster than af
        print("DONE")


        print("\n*** ESMFold design selection (based on higher pLDDT) ***")
        # Design selection per target
        # Track best designs files, sorted to avoid matching problems
        best_esm_models = sorted(glob.glob(os.path.join(pdb_folder, "*_output", "seqs", "*_ESMFold", "*_best_design")))
        #print("best_esm_models")
        #print(best_esm_models)
        all_esm_models = sorted(glob.glob(os.path.join(pdb_folder, "*_output", "seqs", "*_ESMFold")))
        #print("all_esm_models")
        #print(all_esm_models)

        out_dir = Path(pdb_folder) / "Best_designs_ESMFold"
        out_dir.mkdir(parents=True, exist_ok=True)
        #print(f"{out_dir} created")
        out_file = out_dir/ "best_designs_ESMFold.fasta"
        # Check if exists, if yes remove and create it again to avoid overwriting
        if out_file.exists():
            out_file.unlink()  # Remove the file
            #Create an empty file
            out_file.touch()
            print(f"{out_file} created")
        
        # Iterate in pairs to extract seq and pdb matchings
        # Convert to dicts keyed by basename without extension
        # Create lookup: map 4-char prefix → file
        #fasta_dict = {Path(f).name[:4]: f for f in new_fasta_files}
        #best_design_dict = {Path(f).name[:4]: f for f in best_esm_models}
        #folder_designs_dict = {Path(f).name[:4]: f for f in all_esm_models}
        #standard for any pdb name (pdb code, af2, etc)
        #print("fasta_dict")
        fasta_dict = {Path(f).name.split(".")[0]: f for f in new_fasta_files}
        #print(fasta_dict)
        #print("best_design_dict")
        best_design_dict = {Path(f).name.split("_best")[0]: f for f in best_esm_models}
        #print(best_design_dict)
        #print("folder_designs_dict")
        folder_designs_dict = {Path(f).name.split("_ESM")[0]: f for f in all_esm_models}
        #print(folder_designs_dict)
        

        # Get common prefixes
        common_keys_seq = sorted(set(best_design_dict.keys()) & set(fasta_dict.keys()))
        common_keys_pdb = sorted(set(best_design_dict.keys()) & set(folder_designs_dict.keys()))

        # Fasta matching
        if not common_keys_seq:
            print("No matching fasta keys found!")

        # Process seq matching pairs
        for prefix in common_keys_seq:
            best_design = best_design_dict[prefix]
            new_fasta_file = fasta_dict[prefix]
            #print(f"Processing matching prefix for '{prefix}': {best_design}, {new_fasta_file}")
            print(f"Processing matching prefix for '{prefix}' to get sequence")
            extract_matching_sequences(best_design, new_fasta_file, out_file)

        # pdb matching
        if not common_keys_pdb:
            print("No matching pdb keys found!")

        # Process pdb matching pairs
        for prefix in common_keys_pdb:
            best_design = best_design_dict[prefix]
            folder_designs = folder_designs_dict[prefix]
            #print(f"Processing matching prefix '{prefix}': {best_design}, {folder_designs}")
            print(f"Processing matching prefix for '{prefix}' to get structure")
            copy_best_design_pdbs(best_design, folder_designs, out_dir)

        print("DONE")

        return out_file, out_dir

    else:
        exit

if __name__ == "__main__":
    
    # Retrieve args
    args = parse_args()

    # Name args
    pdb_folder = args.pdb_folder
    modelling = str(args.modelling)

    main(pdb_folder, modelling)
