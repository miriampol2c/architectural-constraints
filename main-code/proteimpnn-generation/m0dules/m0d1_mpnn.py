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
    parser.add_argument("--num_seq_per_target", type=int, required=False, default=30, help="Number of sequences to be designed per target")
    parser.add_argument("--sampling_temp", type=float, required=False, default=0.1, help="Sampling temperature")

    return parser.parse_args()


######################################
# FUNCTIONS TO DESIGN MPNN SEQUENCES #
######################################

def wait_for_sbatch(job_output):
    match = re.search(r"Submitted batch job (\d+)", job_output)
    if not match:
        raise RuntimeError("Could not extract job ID from sbatch output.")
    job_id = match.group(1)

    print(f"Waiting for SLURM job {job_id} to finish...")

    while True:
        result = subprocess.run(["squeue", "--job", job_id],
                                capture_output=True, text=True)
        if job_id not in result.stdout:
            print(f"Job {job_id} completed.")
            break
        time.sleep(10)  # Poll every 10 seconds

# Set here for the moment
#generation_script="/gpfs/scratch/bsc08/bsc008609/PROTEINMPNN/ProteinMPNN/protein_mpnn_run.py"
#include seed as an argument to be defined?

def mpnn_design(pdb_folder, num_seq_per_target, sampling_temp):
    #generation_script = os.path.join(os.getcwd(), 'ext_functions', 'protein_mpnn_run.py')
    # para q funcion tiene q tomarse de aqui
    generation_script="/gpfs/scratch/bsc08/bsc008609/PROTEINMPNN/ProteinMPNN/protein_mpnn_run.py"

    # Get list of PDB files
    pdb_files = glob.glob(os.path.join(pdb_folder, "*.pdb"))
    print(f"Found {len(pdb_files)} PDB files\n")
    
    if not pdb_files:
        print(f"No .pdb files found in {pdb_folder}")
        return

    #print("List of PDB files:")
    #for pdb_file in pdb_files:
        #print(pdb_file)

    # Start design
    for pdb_file in pdb_files:
        #print(f"\nProcessing file: {pdb_file}")
        pdb_name = os.path.splitext(os.path.basename(pdb_file))[0]
        print(f"Processing {pdb_name}")
        output_dir = os.path.join(pdb_folder, f"{pdb_name}_MPNN_output")
        os.makedirs(output_dir, exist_ok=True)

        cmd = [
            "python", str(generation_script),
            "--pdb_path", str(pdb_file),
            "--out_folder", str(output_dir),  
            "--num_seq_per_target", str(num_seq_per_target),
            "--sampling_temp", str(sampling_temp),
            "--save_score", "1",
            "--save_probs", "1",
            "--seed", "32",
            "--batch_size", "1"
        ]

        #print("Running command:", " ".join(cmd))  # debug print
        result = subprocess.run(cmd, check=True, capture_output=True, text=True)
        print(result.stdout)
        print(result.stderr)
        #subprocess.run(cmd, check=True)  # check=True raises error if fails


def curate_and_filter_sequences(input_file, output_file_path):
    # Open input file for reading
    with open(input_file, 'r') as file:
        sequences = []
        current_sequence = None
        first_sequence_header = None

        # Read each line in the input file
        for line in file:
            line = line.strip()

            # Check if line is a header
            if line.startswith('>'):
                # Check if there was a previous sequence, and if it passed the filter
                if current_sequence is not None:
                    #print(f'--------------------------------------------------------------------')
                    if any(char in current_sequence for char in 'UXBZJ'):
                        print(f'Sequence(s) to be removed because contain non-canonical amino acids:')
                        # Print the sequence if it contains the specified items
                        print(header, current_sequence)
                    else:
                        # Append the sequence to the list if it does not contain the specified items
                        sequences.append((header, current_sequence))
                    #print(f'--------------------------------------------------------------------')
                        
                # Save the new header
                header = line
                current_sequence = ''

                if first_sequence_header is None:
                    first_sequence_header = os.path.basename(header.split()[0])  # Get the base name of the header

                # Exclude writing the first sequence to the output file
                continue
                # deberia dejar la ref? o añadirla solo a la hora de hacer el aligment??

            # Append sequence lines to the current sequence
            current_sequence += line

        # Check the last sequence in the file
        if current_sequence is not None:
            if any(char in current_sequence for char in 'UXBZJ'): 
                # Print the sequence if it contains the specified items
                print(header, current_sequence)
            else:
                # Append the sequence to the list if it does not contain the specified items
                sequences.append((header, current_sequence))

    # Open output file for writing
    with open(output_file_path, 'w') as file:
        # Write the filtered sequences to the output file
        for i, (header, sequence) in enumerate(sequences):
            if i > 0:
                # Append the first part of the first sequence header to the header of subsequent sequences
                header = f"{first_sequence_header.rstrip(',')}_{header[1:].replace(', ', '_')}"
                # remove dots between digits
                # header = re.sub(r'(?<=\d)\.(?=\d)', '', header) # provisional to remove dots problems in FE
                header = re.sub(r'=', '', header)               # Remove all '=' characters
                file.write(header + '\n')

                # Write the sequence in one line
                file.write(sequence + '\n')

                # Write the sequence with 60 amino acids per line
                #for j in range(0, len(sequence), 60):
                #    file.write(sequence[j:j+60] + '\n')
   
def modify_header(header):
    if '_' in header:
        # Case 1: Header contains at least 7 underscores and includes "relaxed" or "unrelaxed"
        if header.count("_") >= 7 and ("relaxed" in header or "unrelaxed" in header):
            short_header = "_".join(header.split("_")[:8])
        # Case 2: Header contains exactly 7 underscores
        elif header.count("_") == 7:
            short_header = "_".join(header.split("_")[:3])
        # Case 3: Header contains at least 8 underscores and does not include '-'
        elif header.count("_") >= 8 and '-' not in header:
            short_header = "_".join(header.split("_")[:4])
        # Case 4: Header contains both '-' and '_'
        elif '-' in header and '_' in header:
            short_header = '_'.join(header.split("_")[:5])
        # Default case: Return the original header
        else:
            short_header = header
    else:
        # Case: Header does not contain '_'
        short_header = header

    return short_header

def final_multifasta(input_file, final_file_path):
    # Create a list to store modified records
    modified_records = []

    with open(input_file, "r") as handle:
        for record in SeqIO.parse(handle, "fasta"):
            # Modify the header of each record
            modified_header = modify_header(record.description)
            record.description = modified_header
            record.id = modified_header

            # Append the modified record to the list
            modified_records.append(record)

    # Bio.SeqIO.write() to write sequences in FASTA format, by default it wraps lines at 60 characters. 
    # To write each sequence on a single line, you'll need to override the default line-wrapping behavior.
    # Write the modified records to the final file
    # with open(final_file_path, "w") as out_handle:
    #    SeqIO.write(modified_records, out_handle, "fasta")
    with open(final_file_path, "w") as out_handle:
        for record in modified_records:
            out_handle.write(f">{record.id}\n{str(record.seq)}\n")

def rename_mpnn_fasta(pdb_folder):
    
    # Check designs files
    fasta_files = glob.glob(os.path.join(pdb_folder, "*_output", "seqs", "*.fa"))
    
    # Curate designs files
    for fasta_file in fasta_files:
        # Check input fasta file
        absolute_fasta_path = os.path.abspath(os.path.join(os.getcwd(), fasta_file))
        #print("Processing file:", absolute_fasta_path) # useful for debug
        input_fasta_directory = os.path.dirname(absolute_fasta_path)
        # Define output fasta file
        output_fasta_name = os.path.basename(absolute_fasta_path).replace(".fa", "_fullheaders.fa")
        output_fasta_path = os.path.join(input_fasta_directory, output_fasta_name)
        absolute_output_path = os.path.abspath(os.path.join(os.getcwd(), output_fasta_path))
        # Final file
        final_file_name = os.path.basename(absolute_output_path).replace("_fullheaders.fa", ".fasta")
        final_file_path = os.path.join(input_fasta_directory, final_file_name)
        
        # Apply auxiliar_functions
        curate_and_filter_sequences(absolute_fasta_path, output_fasta_path)
        final_multifasta(absolute_output_path, final_file_path)
        
def split_multi_fasta(input_file, output_folder):
    # Create output folder if it does not exist
    os.makedirs(output_folder, exist_ok=True)
    
    # Save individual files
    for record in SeqIO.parse(input_file, "fasta"):
        output_file = os.path.join(output_folder, f"{record.id}.fasta")
        SeqIO.write(record, output_file, "fasta")

    #print(f"Multifasta completely divided. Individual fasta files saved in {output_folder}.")
   

# MAIN

def main(pdb_folder, num_seq_per_target, sampling_temp):
    
    # Design sequences for all the targets filtered
    print("\n*** ProteinMPNN design ***")
    mpnn_design(pdb_folder, num_seq_per_target, sampling_temp)
    #print("DONE")

    # Curate mpnn designs output before predict the structures (renaming)
    print("*** Renaming multifasta files ***")
    rename_mpnn_fasta(pdb_folder)
    print("DONE")
    
    # Once curated, the extension of the fasta files has changed. Sorted to avoid matching problems, because glob.glob does not guarantee sorting
    new_fasta_files = sorted(glob.glob(os.path.join(pdb_folder, "*_output", "seqs", "*.fasta")))
    print("\n*** Splitting multifasta files ***")
    for new_fasta_file in new_fasta_files:
        # splitting is needed just for af2, but was included by default
        # First split multifasta to fold the sequences, as the inputs are single sequences (introduce greasy?)
        # Define where individual fastas will be saved: same directory as where multifasta is
        multifasta_dir = os.path.dirname(new_fasta_file)
        output_folder = os.path.join(multifasta_dir, "Individual_designs")
        # Call split function
        split_multi_fasta(new_fasta_file, output_folder)
        #print("DONE")
        # Get list of split fasta files
        #individual_fastas = glob.glob(os.path.join(output_folder, "*.fasta"))
    print("DONE")
	
if __name__ == "__main__":

    # Retrieve args
    args = parse_args()

    # Name args
    pdb_folder = args.pdb_folder
    num_seq_per_target = str(args.num_seq_per_target)
    sampling_temp = str(args.sampling_temp)

    main(pdb_folder, num_seq_per_target, sampling_temp)
