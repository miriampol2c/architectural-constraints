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
# FUNCTIONS TO PREDICT AF2/AF3 STRUCTURES #
###########################################

PDB_RESIDUE_TO_AA = {
    'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C',
    'GLU': 'E', 'GLN': 'Q', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I',
    'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F', 'PRO': 'P',
    'SER': 'S', 'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V'
}

def extract_sequences_from_pdb(pdb_file):
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("structure", pdb_file)
    chain_sequences = {}

    for model in structure:
        for chain in model:
            chain_id = chain.id
            seq_list = []
            for residue in chain:
                resname = residue.get_resname()
                if resname in PDB_RESIDUE_TO_AA:
                    seq_list.append(PDB_RESIDUE_TO_AA[resname])
            chain_sequences[chain_id] = "".join(seq_list)
        break  # Only first model

    return chain_sequences

def extract_pdb_sequences(input_dir):
    input_dir = Path(input_dir)
    output_file = input_dir / "Best_designs.fasta"

    pdb_files = list(input_dir.glob("*.pdb"))
    if not pdb_files:
        print(f"No PDB files found in {input_dir}")
        return

    seq_records = []
    for pdb_file in pdb_files:
        sequences = extract_sequences_from_pdb(pdb_file)
        if not sequences:
            print(f"No sequences found in {pdb_file.name}")
            continue
        for chain_id, seq in sequences.items():
            seq_id = f"{pdb_file.stem}"
            seq_records.append(SeqRecord(Seq(seq), id=seq_id, description=""))

    if seq_records:
        SeqIO.write(seq_records, output_file, "fasta")
        print(f"Saved {len(seq_records)} sequences to {output_file}")
    else:
        print("No sequences extracted from any PDB files.")

# function to parse individual fastas as json, using a template
def fasta2json(folder):
    individual_fastas = glob.glob(os.path.join(folder, "*.fasta"))

    # Create output directory next to the input folder
    folder_path = Path(folder)
    output_dir = folder_path.parent / "Individual_jsons"
    output_dir.mkdir(parents=True, exist_ok=True)

    for fasta_path in individual_fastas:  # This is already a full path
        fasta_path = Path(fasta_path)

        # Read the header and sequence from the FASTA file
        with open(fasta_path, "r") as f:
            lines = f.readlines()
            header = lines[0].strip().lstrip(">")
            sequence = "".join(line.strip() for line in lines[1:])  # Join all lines after header

        # JSON template
        template = {
            "name": "",
            "sequences": [
                {
                    "protein": {
                        "id": ["A"],
                        "sequence": ""
                    }
                }
            ],
            "modelSeeds": [1],
            "dialect": "alphafold3",
            "version": 1
        }

        # Fill in the template
        json_obj = template.copy()
        json_obj["name"] = header
        json_obj["sequences"][0]["protein"]["sequence"] = sequence

        # Output JSON file path
        output_file = output_dir / f"{fasta_path.stem}.json"
        with open(output_file, "w") as out_f:
            json.dump(json_obj, out_f, indent=4)
    
    return output_dir

def wait_for_sbatch(job_output, poll_interval=10):
    """
    Wait until a submitted SLURM job finishes.
    Works with both default and --parsable sbatch output.
    """

    # Clean output (remove whitespace, keep only last line)
    job_output = job_output.strip().splitlines()[-1]

    # Extract job ID
    match = re.search(r"(\d+)", job_output)
    if not match:
        raise RuntimeError(f"Could not extract job ID from sbatch output: {job_output}")
    job_id = match.group(1)

    print(f"Waiting for SLURM job {job_id} to finish...")

    while True:
        result = subprocess.run(["squeue", "--job", job_id],
                                capture_output=True, text=True)
        if job_id not in result.stdout:
            print(f"Job {job_id} completed.")
            break
        time.sleep(poll_interval)

def folding_AF2(folder):
    # Activate temporarily the specific AF2 module, must be included all in the same shell to preserve the module loaded
    #subprocess.run("module purge && module load singularity alphafold/2.3.2 && module load openmpi/4.1.5-gcc", shell=True, executable="/bin/bash")

    jobids = []

    # fastas already splitted by default
    individual_fastas = sorted(glob.glob(os.path.join(folder, "*_output", "seqs", "Individual_designs", "*.fasta")))
    print(f'Predicting AlphaFold2 structures for designs contained in {folder}')
    for fasta_file in individual_fastas:
        # define out dir (pred_dir)
        fasta_dir = os.path.dirname(fasta_file)
        fasta_dir_partial = os.path.dirname(fasta_dir)  # one level up, which is seqs 
        pred_dir = os.path.join(fasta_dir_partial, "AF2")

        filename = os.path.splitext(os.path.basename(fasta_file))[0]
        targetname = "_".join(os.path.splitext(os.path.basename(fasta_file))[0].split("_")[:1])
        
        # Define the command exactly as you would run it in Bash
        logdir = os.getcwd()  # get current absolute path
        slurm_command = f"""#!/bin/sh
#SBATCH -A
#SBATCH -q
#SBATCH --account=bsc08
#SBATCH --job-name=af2
#SBATCH --output={logdir}/debug/AF2/{filename}/af2_%j.out
#SBATCH --error={logdir}/debug/AF2/{filename}/af2_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --qos=acc_bscls

module purge
module load singularity
module load alphafold/2.3.2
module load openmpi/4.1.5-gcc

echo "Start at $(date)"

bsc_alphafold \\
--fasta_paths={fasta_file} \\
--output_dir={pred_dir} \\
--use_gpu_relax=False \\
--max_template_date=2023-12-31 \\
--model_preset=monomer

echo "End at $(date)"
"""
            
        with open(f"submit_af2_{filename}.sh", "w") as f:
            f.write(slurm_command)

        # Submit the job and capture output
        result = subprocess.run(["sbatch", f"submit_af2_{filename}.sh"], capture_output=True, text=True)
        # como monitorizo cuando acaba esto para seleccionar los diseños, lo mejor seguir usando dependencies?
        # Parse and store the job ID
        if result.returncode == 0:
            output = result.stdout.strip()
            print(f"Submission output: {output}")
            if output.startswith("Submitted batch job"):
                jobid = output.split()[-1]
                jobids.append(jobid)
            else:
                print(f"Unexpected sbatch output: {output}")
        else:
            print(f"Error submitting job: {result.stderr}")
    
    # Return the dependency string to use in follow-up sbatch
    if jobids:
        return f"--dependency=afterok:{':'.join(jobids)}"
    else:
        return None  # Or raise an exception depending on how critical the failure isdependency = f"--dependency=afterok:{':'.join(jobids)}"


def folding_AF3(folder):

    jobids = []

    output_folders = sorted(glob.glob(os.path.join(folder, "*_output", "seqs", "Individual_designs")))
    
    for output_folder in output_folders:
        # Go two levels up
        out_path = os.path.dirname(os.path.dirname(output_folder))  # until _output...
        target_name = os.path.basename(out_path).replace("_output", "")
        print(out_path)
        print(target_name)
        json_folder = fasta2json(output_folder) # the output is another output folder containing all the json files (input af3)
        
        print(f'Predicting AlphaFold3 structures for sequences contained in {output_folder}')
        if json_folder:
            pred_dir = os.path.join(folder, f"{target_name}_output", "seqs", "AF3")
            os.makedirs(pred_dir, exist_ok=True)
            logdir = os.getcwd()  # get current absolute path
            slurm_command = f"""#!/bin/sh
#SBATCH -A
#SBATCH -q
#SBATCH --account=bsc08
#SBATCH --job-name=af3
#SBATCH --output={logdir}/debug/AF3/{target_name}/af3_%j.out
#SBATCH --error={logdir}/debug/AF3/{target_name}/af3_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=40
#SBATCH --gres=gpu:1
#SBATCH --qos=acc_bscls

module purge
module load cuda/12.6
module load singularity alphafold/3.0.0

fasta={json_folder} # as input, a directory containing jsons (full path)
output={pred_dir}
weights="/gpfs/projects/bsc08/shared_tools/alphafold3_weights/"

echo "Start at $(date)"
bsc_alphafold $fasta $output $weights

echo "End at $(date)"
"""                         
            with open(f"submit_af3_{target_name}.sh", "w") as f:
                f.write(slurm_command)

            # Submit the job and capture output
            result = subprocess.run(["sbatch", f"submit_af3_{target_name}.sh"], capture_output=True, text=True)
            # como monitorizo cuando acaba esto para seleccionar los diseños, lo mejor seguir usando dependencies?
            # Parse and store the job ID --> one job per target (all tha samples together)
            if result.returncode == 0:
                output = result.stdout.strip()
                print(f"Submission output: {output}")
                if output.startswith("Submitted batch job"):
                    jobid = output.split()[-1]
                    jobids.append(jobid)
                else:
                    print(f"Unexpected sbatch output: {output}")
            else:
                print(f"Error submitting job: {result.stderr}")
    
    # Return the dependency string to use in follow-up sbatch
    if jobids:
        return f"--dependency=afterok:{':'.join(jobids)}"
    else:
        return None  # Or raise an exception depending on how critical the failure is


# MAIN

def main(pdb_folder, modelling):

    #Once curated, the extension of the fasta files has changed. Sorted to avoid matching problems, because glob.glob does not guarantee sorting
    #new_fasta_files = sorted(glob.glob(os.path.join(pdb_folder, "*_output", "seqs", "*.fasta")))
    
    # AF2 modelling (default modelling) 
    if modelling == "AF2":
        #new_fasta_files = sorted(glob.glob(os.path.join(pdb_folder, "*_output", "seqs", "Individual_designs", "*.fasta")))
        dependency_list_af2 = folding_AF2(pdb_folder)
        if dependency_list_af2:
            print(f"Will use dependency: {dependency_list_af2}")
            # Complete folding is the requiriment to select designs according to the pLDDT
            # use new_pdb_folder as input, the folder containing filtered pdbs 
            #print(new_pdb_folder)
            # make sure that the next sh is executable from here (until include code in here)
            selection_script = os.path.join(os.getcwd(), 'ext_scripts', 'BestDesignsAF2.sh')
            subprocess.run(["chmod", "+x", selection_script])
            result = subprocess.run([
            "sbatch",
            "--parsable", dependency_list_af2, 
            f"--export=ALL,new_pdb_folder={pdb_folder}",
            selection_script], capture_output=True, text=True)
             
            # Wait for job to finish (the alternative way is using dependencies, as we do for the folding)
            # but this time is just one job
            print("Raw sbatch output:", result.stdout)
            wait_for_sbatch(result.stdout)

            # Need to derive the multifasta from the best designs and built the msa
            folder_best_designs = (os.path.join(pdb_folder, "Best_designs_AF2", "per_target_and_sample"))
            extract_pdb_sequences(folder_best_designs)
            fasta_best_designs = (os.path.join(folder_best_designs, "Best_designs.fasta"))

            [os.remove(f) for f in glob.glob("submit_af2*.sh")]

            return fasta_best_designs, folder_best_designs
        else:
            print("AF2 jobs did not submit correctly. Cannot schedule dependent job.")

    # AF3 modelling 
    if modelling == "AF3":
        
        # DESDE AQUI CAMBIAR, EL input es el folder, hago glob glob en la funcion
        dependency_list_af3 = folding_AF3(pdb_folder)
        if dependency_list_af3:
            print(f"Will use dependency: {dependency_list_af3}")
            
            # Complete folding is the requiriment to select designs according to the pLDDT
            selection_script = os.path.join(os.getcwd(), 'ext_scripts', 'BestDesignsAF3.sh')
            subprocess.run(["chmod", "+x", selection_script])
            result = subprocess.run([
            "sbatch", "--parsable", dependency_list_af3, 
            f"--export=ALL,new_pdb_folder={pdb_folder}",
            selection_script], capture_output=True, text=True)
            
            # wait for folding to finish
            print("Raw sbatch output:", result.stdout)
            wait_for_sbatch(result.stdout)

            # Need to derive the multifasta from the best designs and built the msa
            folder_best_designs = (os.path.join(pdb_folder, "Best_designs_AF3", "per_target_and_sample"))
            extract_pdb_sequences(folder_best_designs)
            fasta_best_designs = (os.path.join(folder_best_designs, "Best_designs.fasta"))

            # remove extra sh 
            [os.remove(f) for f in glob.glob("submit_af3*.sh")]
            
            return fasta_best_designs, folder_best_designs

        else:
            print("AF3 jobs did not submit correctly. Cannot schedule dependent job.")
        

if __name__ == "__main__":
    
    # Retrieve args
    args = parse_args()

    # Name args
    pdb_folder = args.pdb_folder
    modelling = str(args.modelling)

    main(pdb_folder, modelling)