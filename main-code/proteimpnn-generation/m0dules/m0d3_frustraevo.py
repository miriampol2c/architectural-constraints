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
import pandas as pd

def parse_args():
    parser = argparse.ArgumentParser(description="Run ProteinMPNN design pipeline on a folder of PDB files.")
    parser.add_argument("--pdb_folder", type=str, required=True, help="Absolute path to the folder containing .pdb files")
    parser.add_argument("--pdbs", type=str, required=True, help="Absolute path to the folder containing .pdb files")
    parser.add_argument("--multifasta", type=str, required=False, help="Multifasta file")
    parser.add_argument("--jobid", type=str, required=True, help="Name to track the FrustraEvo results")
    parser.add_argument("--ref", type=str, required=False, help="Choose the protein reference for the FrustraEvo analysis") # if no native id protein provided, default would be taking the first designed protein of the alignment once alph sorted
    parser.add_argument("--modelling", type=str, required=True, help="Folding algorithm used to predict designs") 
    
    return parser.parse_args()

def process_ref(ref, new_pdb_folder, out_dir, out_file):
    if ref:
        # Paths
        ref_pdb = Path(new_pdb_folder) / f"{ref}.pdb"
        ref_fasta = Path(new_pdb_folder) / f"{ref}_MPNN_output" / "seqs" / f"{ref}.fa"

        if not Path(ref_fasta).is_file():
            raise FileNotFoundError(f'The selected "Ref": {ref} does not exist!')

        # Copy PDB to out_dir
        out_dir = Path(out_dir)
        out_dir.mkdir(parents=True, exist_ok=True)
        shutil.copy(ref_pdb, out_dir / ref_pdb.name)

        print(f"Copied PDB ref to {out_dir}")

        # Extract first sequence from multifasta
        records = list(SeqIO.parse(ref_fasta, "fasta"))
        if not records:
            raise ValueError(f"No sequences found in {ref_fasta}")

        # Take first sequence
        first_seq = records[0]
        seq_str = str(first_seq.seq)  # <-- plain string, only the sequence

        # Format FASTA entry with your ref variable as header
        fasta_entry = f">{ref}\n{seq_str}\n"

        # Prepend to existing file content
        if os.path.exists(out_file):
            with open(out_file, "r") as f:
                existing_content = f.read()
        else:
            existing_content = ""

        with open(out_file, "w") as f:
            f.write(fasta_entry)
            f.write(existing_content)

        print(f"Copied sequence ref to {out_file}")

def prepare_input_FE_mafft_msa(input_fasta, ref_id=None, output_fasta=None):
    """
    Align sequences from an input FASTA file using MAFFT, saving the aligned result to output_fasta.
    Also moves all .pdb files in the same directory as input_fasta to a new subdirectory 'pdbs'.

    Args:
        input_fasta (str or Path): Path to the input FASTA file containing unaligned sequences.
        output_fasta (str or Path, optional): Path to write the aligned sequences.
    """

    input_fasta = Path(input_fasta)
    input_dir = input_fasta.parent

    # Default output path in the same directory
    if output_fasta is None:
        output_fasta = input_dir / "MSA.fasta"
    else:
        output_fasta = Path(output_fasta)

    # Run MAFFT
    print(f"Running MAFFT on {input_fasta}...")
    mafft_cline = MafftCommandline(input=str(input_fasta), maxiterate=1000, auto=True, quiet=False)
    stdout, stderr = mafft_cline()
    alignment = AlignIO.read(StringIO(stdout), "fasta")
    #AlignIO.write(alignment, output_fasta, "fasta")
    #print(f"MSA written to {output_fasta}")

    # Separate and sort
    if ref_id:
        ref_seq = [rec for rec in alignment if rec.id == ref_id]
        other_seqs = [rec for rec in alignment if rec.id != ref_id]
        sorted_alignment = ref_seq + sorted(other_seqs, key=lambda record: record.id.lower())
    else:
        sorted_alignment = sorted(alignment, key=lambda record: record.id.lower())

    # Write alignment
    sorted_alignment = MultipleSeqAlignment(sorted_alignment)
    AlignIO.write(sorted_alignment, output_fasta, "fasta")
    print(f"MSA written to {output_fasta}")

    # Move .pdb files to 'pdbs/' subfolder
    pdb_dir = input_dir / "pdbs"
    pdb_dir.mkdir(exist_ok=True)

    pdb_files = list(input_dir.glob("*.pdb"))
    for pdb_file in pdb_files:
        dest = pdb_dir / pdb_file.name
        shutil.move(str(pdb_file), str(dest))
        #print(f"Moved {pdb_file.name} to {pdb_dir}")


def rebuild_msa(template_msa_file, replicates_file, out_file):
    def base_id(header):
        """Extract base identifier (remove _T suffix and spaces)."""
        return header.split()[0].split("_")[0] # adding T0 instead of just T because could be confusing with the TED domaind annotation

    # --- Read input files ---
    template_records = list(SeqIO.parse(template_msa_file, "fasta"))
    replicates = list(SeqIO.parse(replicates_file, "fasta"))

    new_records = []
    invalid = False  # flag to stop building if any issue appears
    problems = []  # collect issues to report

    for t_rec in template_records:
        b_id = base_id(t_rec.id)
        matched_replicates = [r for r in replicates if base_id(r.id) == b_id]

        # --- Case 1: no replicate found ---
        if not matched_replicates:
            problems.append(f" No replicate found for {t_rec.id}")
            invalid = True
            continue

        template_ungapped = str(t_rec.seq).replace("-", "")
        template_len = len(template_ungapped)

        for r in matched_replicates:
            ungapped_seq = str(r.seq).replace("-", "")
            seq_len = len(ungapped_seq)

            # --- Case 2: length mismatch ---
            if seq_len != template_len:
                problems.append(
                    f" Length mismatch for {r.id} "
                    f"(replicate {seq_len} vs template {template_len})"
                )
                invalid = True
                continue

            # --- Apply template gap pattern ---
            new_seq = []
            j = 0
            for c in str(t_rec.seq):
                if c == "-":
                    new_seq.append("-")
                else:
                    new_seq.append(ungapped_seq[j])
                    j += 1

            r.seq = Seq("".join(new_seq))
            new_records.append(r)

    # --- Abort if any template or replicate was invalid ---
    if invalid:
        print(f"Skipping MSA build for {template_msa_file} due to inconsistencies:")
        for p in problems:
            print(f"   {p}")
        return  # do not write MSAs

    # --- Build MSA ---
    SeqIO.write(new_records, out_file, "fasta")
    print(f" Successfully built MSA → {out_file} ({len(new_records)} sequences)")


def prepare_input_FE_template_based(input_fasta, refid=None, jobid=None, output_fasta=None):
    """
    Align sequences from an input FASTA file using as a template the one from frustraevo when evaluating same family but with native members, saving the aligned result to output_fasta.
    Also moves all .pdb files in the same directory as input_fasta to a new subdirectory 'pdbs'.

    Args:
        input_fasta (str or Path): Path to the input FASTA file containing unaligned sequences.
        output_fasta (str or Path, optional): Path to write the aligned sequences.
    """

    input_fasta = Path(input_fasta)
    input_dir = input_fasta.parent # this is Best_designs_Xmodel, fasta file

    # Default output path in the same directory
    if output_fasta is None:
        output_fasta = input_dir / "MSA.fasta"
    else:
        output_fasta = Path(output_fasta)

    # Build template-based msa
    # first, find the template and copy within input_dir
    # we must copy the MSA_Clean_aux from native runs, because in these msas, 
    # the positions of the MSA that are gaps in the reference have not been yet trimmed by frustraevo
    # from such copy, derive the MSA to be used in our frustarevo --< pending Miguel
    #clean_aux_dir = "/gpfs/scratch/bsc08/bsc008609/CSA_enzymes/FrustMPNN/funfams_feb2026/clean_aux_msas/"
    # the one provided by miguel with msas from all families
    clean_aux_dir = "/gpfs/scratch/bsc08/bsc008092/CATH/enzymes/MIRIAM/MSA_Clean_aux_files/"
    # in here there are multiples msas look for 1 in specific
    clean_aux_file_template = os.path.join(clean_aux_dir, f"MSA_Clean_aux_{jobid}.fasta") # jobid from arguments
    # If you want to check that it exists:
    if not os.path.exists(clean_aux_file_template):
        raise FileNotFoundError(f"{clean_aux_file_template} not found")
    # input_fasta is the multifasta contaiing best designs
    rebuild_msa(clean_aux_file_template, input_fasta, output_fasta)

    # Move .pdb files to 'pdbs/' subfolder
    pdb_dir = input_dir / "pdbs"
    pdb_dir.mkdir(exist_ok=True)

    pdb_files = list(input_dir.glob("*.pdb"))
    for pdb_file in pdb_files:
        dest = pdb_dir / pdb_file.name
        shutil.move(str(pdb_file), str(dest))
        #print(f"Moved {pdb_file.name} to {pdb_dir}")
        
def run_FE(out_dir, jobid, ref=None):
    out_dir = Path(out_dir)
    pdb_folder = out_dir / "pdbs" 
    #print(pdb_folder)

    # --- Check if PDB folder exists and contains at least 3 files ---
    if not pdb_folder.exists() or not pdb_folder.is_dir():
        raise FileNotFoundError("PDB folder not found.")
    pdb_files = list(pdb_folder.glob("*.pdb"))
    if len(pdb_files) < 3: #3 is the minimum to avoid an error when measuring conservation in frustraevo
        raise ValueError("PDB folder must contain at least 3 PDB files.")

    # --- Check if MSA FASTA file exists and contains at least 3 sequences ---
    out_file_list = list(out_dir.glob("MSA.fasta"))
    if not out_file_list:
        raise FileNotFoundError("No MSA fasta file found in the output directory.")
    
    out_file = out_file_list[0]
    #print(out_file)
    msa_seqs = list(SeqIO.parse(out_file, "fasta"))
    if len(msa_seqs) < 3:
        raise ValueError("MSA file must contain at least 3 sequences.")

    # -- Set ref from first sequence of the MSA  if not provided
    if not ref:
        alignment = AlignIO.read(out_file, "fasta")
        ref = alignment[0].id
        print(ref)

    # Script paths
    fe_script = "/gpfs/scratch/bsc08/bsc008092/FrustraEvo/FrustraEvo-updates-testing/run_logo.py"
    r_script = "/gpfs/scratch/bsc08/bsc008092/FrustraEvo/FrustraEvo-updates-testing/Scripts"

    # Slurm command
    logdir = os.getcwd()
    slurm_command = f"""#!/usr/bin/env bash
#SBATCH --account=bsc08
#SBATCH --job-name=FE
#SBATCH --output={logdir}/debug/FE/FE_%j.out
#SBATCH --error={logdir}/debug/FE/FE_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --qos=gp_bscls


# to use gpp modules from acc
source /apps/modules/gpp.sh

#module purge
#export GREASY_LOGFILE=./debug/{ref}.log
module load intel
#module load greasy
module load anaconda/2023.07
module load R/4.3.2
source activate frustratometer-reqs
export PATH=/apps/GPP/ANACONDA/2023.07/envs/frustratometer-reqs/bin:$PATH

echo "Start at $(date)"
echo {ref}

python3 {fe_script} --JobId {jobid} --fasta {out_file} --RPath {r_script} --ref {ref} --pdb_db {pdb_folder} --SeqDist 3 --cmaps yes --Graphics TRUE

#as output will be saved in the current directory, move to the best designs folder 
mv FrustraEvo_{jobid} {out_dir}

echo "End at $(date)"
"""

    slurm_file = f"submit_FE_{ref}.sh"
    with open(slurm_file, "w") as f:
        f.write(slurm_command)
    
    print("jobid sent to the GPP queue")
    subprocess.run(["sbatch", slurm_file], text=True)
    print("WAIT FOR IT!")



# MAIN

def main(pdb_folder, out_dir, out_file, jobid, ref, modelling):

    # Include native ref within the FE input folder if provided 
    if ref:
        print("\n*** Including native reference ***")
        process_ref(ref, pdb_folder, out_dir, out_file)
        print("DONE")
            
        # Once folding and designs selection done, prepare FrustraEvo input and run
        if modelling == "ESMFold":
            fasta_best_designs = (os.path.join(pdb_folder, "Best_designs_ESMFold", "best_designs_ESMFold.fasta"))
        if modelling == "AF2":
            fasta_best_designs = (os.path.join(pdb_folder, "Best_designs_AF2", "per_target_and_sample", "Best_designs.fasta"))
        if modelling == "AF3":
            fasta_best_designs = (os.path.join(pdb_folder, "Best_designs_AF3", "per_target_and_sample", "Best_designs.fasta"))

        print("\n*** Preparing FrustraEvo input ***")
        # this is what we used to do
        #prepare_input_FE_mafft_msa(fasta_best_designs, ref)
        # currently using template-based msas (to follow the format of the native msa output from frustraevo)
        print(" < Currently trying to build a template-based MSA on top of the native pattern >")
        prepare_input_FE_template_based(fasta_best_designs, ref, jobid)

        msa_path = Path(fasta_best_designs)
        msa_for_fe = msa_path.parent / "MSA.fasta"
        if msa_for_fe.exists():
            print("DONE")
            print("\n*** Running FrustraEvo ***")
            run_FE(out_dir, jobid, ref) # use best designs directory as input
            #print("DONE") # termina en gpp
            # Once finished, remove sh and pdf --> why its outside?
            [os.remove(f) for f in glob.glob(f"submit_FE_{ref}.sh")]
            [os.remove(f) for f in glob.glob("Rplot*.*")] # no funciona esto 
        else:
            raise ValueError("It is not possible to run FrustraEvo because there is no MSA available")


    else:
        raise ValueError("IMPORTANT: You must provide a reference to run FrustraEvo")
        

if __name__ == "__main__":

    # Retrieve args
    args = parse_args()

    # Name args
    modelling = args.modelling
    pdb_folder = args.pdb_folder
    pdbs = args.pdbs
    multifasta = args.multifasta
    jobid = str(args.jobid)
    ref = args.ref

    main(pdb_folder, pdbs, multifasta, jobid, ref, modelling)