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
    parser.add_argument("--jobid", type=str, required=True, help="Name to track the FrustraEvo results")
    

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
    
def run_ESMFold(fasta_file, model, jobid):

    # Input file path
    input_file_path = os.path.abspath(fasta_file)
    # Directory where the input file is located
    directory_path = os.path.dirname(input_file_path)
    # Extract base filename (without extension)
    #base_filename, _ = os.path.splitext(os.path.basename(input_file_path))
    # Define the output directory name
    output_dir_name = f"{jobid}_ESMFold"
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
        prefix = jobid
        with open(f"{output_dir}/{prefix}_plddts", "w") as f:
            for short_sequence_id, plddt in plddts:
                f.write(f"{short_sequence_id}: {plddt}\n")

        # Write sorted pLDDT values to output file
        with open(f"{output_dir}/{prefix}_plddts_sorted", "w") as f:
            for short_sequence_id, plddt in sorted(plddts, key=lambda x: x[1], reverse=True):
                f.write(f"{short_sequence_id}: {plddt}\n")
                
        # Write best designs (one per target) in an individual file if possible
        best_per_basename = {}
        for short_sequence_id, plddt in plddts:
            basename = short_sequence_id.split("_")[0]

            if basename not in best_per_basename:
                best_per_basename[basename] = (short_sequence_id, plddt)
            else:
                _, best_plddt = best_per_basename[basename]
                if plddt > best_plddt:
                    best_per_basename[basename] = (short_sequence_id, plddt)

        # Write best designs file
        with open(f"{output_dir}/{prefix}_best_designs_ESMFold", "w") as f:
            for basename, (short_sequence_id, plddt) in sorted(best_per_basename.items()):
                f.write(f"{basename}\t{short_sequence_id}: {plddt}\n")

        # derive multifasta including best designs
        out_dir = Path(directory_path) / "Best_designs_ESMFold"
        out_dir.mkdir(parents=True, exist_ok=True)
        print(f"{out_dir} created")
        out_file = out_dir/ "best_designs_ESMFold.fasta"
        # Check if exists, if yes remove and create it again to avoid overwriting
        if out_file.exists():
            out_file.unlink()  # Remove the file
        else:
            out_file.touch()
            print(f"{out_file} created")
        
        # Collect all best short IDs
        best_ids = {short_id for (_, (short_id, _)) in best_per_basename.items()}

        # Open input FASTA and output FASTA
        with open(input_file_path, "r") as infile, open(out_file, "w") as outfile:
            for record in SeqIO.parse(infile, "fasta"):

                # If your short_sequence_id is EXACTLY the same as record.id:
                if record.id in best_ids:
                    # write fasta
                    SeqIO.write(record, outfile, "fasta")
                    # copy corresponding pdb file
                    pdb_source = Path(output_dir) / f"{record.id}.pdb"
                    pdb_destination = out_dir / f"{record.id}.pdb"

                    if pdb_source.exists():
                        shutil.copy2(pdb_source, pdb_destination)
                    else:
                        print(f"WARNING: {pdb_source} not found")

                # If short_sequence_id is only PART of the header, use this instead:
                # if any(short_id in record.id for short_id in best_ids):
                #     SeqIO.write(record, outfile, "fasta")

        print(f"Best designs FASTA written to {out_file}")
        

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
                print(f"Found: {header}")
                fo.write(f">{header}\n")
                # Wrap sequence lines
                #chunks = [body[i:i+line_width] for i in range(0, len(body), line_width)]
                #fo.write("\n".join(chunks) + "\n")
                # Write full sequence in a line
                fo.write(body + "\n")  # Write the entire sequence in a single line
                found += 1

    print(f"Done. {found} matching sequences written to: {output_file}")


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

    if verbose:
        print(f"Total sequences to copy: {len(sequence_names)}")

    # Copy files
    copied = 0
    missing = 0
    for name in sequence_names:
        source_file = source_dir / f"{name}.pdb"
        dest_file = destination_dir / f"{name}.pdb"
        if source_file.exists():
            shutil.copy2(source_file, dest_file)
            copied += 1
            if verbose:
                print(f"Copied: {name}.pdb")
        else:
            missing += 1
            if verbose:
                print(f"Missing: {name}.pdb in {source_dir}")

    print(f"Done. {copied} PDB files copied, {missing} missing.")

# MAIN

def main(pdb_folder, modelling, jobid):
    
    new_fasta_files = sorted(glob.glob(os.path.join(pdb_folder, "*all_targets.fasta")))
    #print(new_fasta_files)
    # ESMFold modelling
    if modelling == "ESMFold":

        # Load the model calling the function
        #print("\n*** ESMFold structure prediction ***")#
        model = load_ESMFold()
        # Run the model for each fasta
        for fasta_file in new_fasta_files:
            run_ESMFold(fasta_file, model, jobid) # best design selection included
            # output from here should not be a dependency list, as esmfold is faster than af
    else:
        exit

if __name__ == "__main__":
    
    # Retrieve args
    args = parse_args()

    # Name args
    pdb_folder = args.pdb_folder
    modelling = str(args.modelling)
    jobid = args.jobid

    main(pdb_folder, modelling, jobid)
