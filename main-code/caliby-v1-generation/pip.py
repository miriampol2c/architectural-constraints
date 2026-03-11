import argparse
import os
from os.path import basename
import time
from glob import glob
import pandas as pd
from pathlib import Path
import subprocess
import os
import shutil
import re
from Bio import SeqIO
from Bio.Seq import Seq
from collections import defaultdict

def parse_args():
    parser = argparse.ArgumentParser(description="Run Caliby on a PDBs directory.")
    parser.add_argument("--pdb_folder", type=str, required=True, help="Absolute path to the folder containing .pdb files")
    parser.add_argument("--num_ensembles_per_target", type=int, required=False, default=32, help="Number of ensembles to be generated per target")
    parser.add_argument("--num_seqs_per_target", type=int, required=False, default=30, help="Number of sequences to be designed per target, conditioned on ensembles")
    #parser.add_argument("--sampling_temp", type=float, required=False, default=0.1, help="Sampling temperature")
    #parser.add_argument("--modelling", type=str, required=False, default="AF2", help="Structure prediction method: AF2, AF3, or ESMFold")
    parser.add_argument("--jobid", type=str, required=True, help="Name to track the FrustraEvo results")
    parser.add_argument("--ref", type=str, required=True, help="Choose the protein reference for the FrustraEvo analysis") # if no native id protein provided, default would be taking the first designed protein of the alignment once alph sorted
    return parser.parse_args()

def fix_pdb_charge(pdb_path):
    fixed_lines = []

    with open(pdb_path, "r") as f:
        for line in f:
            if line.startswith(("ATOM", "HETATM")):
                charge = line[78:80]
                if "?" in charge:
                    line = line[:78] + " 0" + line[80:]
            fixed_lines.append(line)

    with open(pdb_path, "w") as f:
        f.writelines(fixed_lines)


def fix_all_pdbs_in_dir(input_dir):
    input_dir = Path(input_dir)

    pdb_files = list(input_dir.glob("*.pdb"))

    print(f"Found {len(pdb_files)} PDB files in {input_dir}")

    for pdb_file in pdb_files:
        fix_pdb_charge(pdb_file)

    print("All PDB files cleaned.")

# Function to clean input pdbs using Caliby's original code
# example how to run: cleaned_dir_pdbs = clean_input_pdbs("/.../new/native_ras")

def clean_input_pdbs(input_dir):
    
    # in here the input is the original folder containing structures in .pdb format
    input_dir = Path(input_dir) # to be sure we have the path

    # Build output directory name: <basename>_clean
    dir_basename = basename(str(input_dir))
    output_dir = input_dir.parent / f"{dir_basename}_clean" / "cifs"  # by default generates cif files, save in pdb just in case:

    output_dir.mkdir(parents=True, exist_ok=True)

    # caliby function to clean pdb cannot convert ' ?' to an integer (? present in some af models) --> so crash
    # add an additional function prior to that one within clean_input_pdbs function
    fix_all_pdbs_in_dir(input_dir)

    # including .resolve for path (to get absolute path in all cases) removes ambiguity working with subprocess and hydra
    cmd1 = [
        env_python,
        f"{caliby_dir}/caliby/data/preprocessing/atomworks/clean_pdbs.py",
        f"input_cfg.pdb_dir={input_dir.resolve()}",
        "num_workers=4",
        f"out_dir={output_dir.resolve()}",
        "hydra.run.dir=."] # not sure if this for the log
    
    #print("Running:", " ".join(map(str, cmd1)))
    result = subprocess.run(cmd1, check=True, capture_output=True, text=True)
    print(result.stdout)
    print(result.stderr)

    out_dir_pdbs = input_dir.parent / f"{dir_basename}_clean" / "pdbs"
    out_dir_pdbs.mkdir(parents=True, exist_ok=True)

    # maybe put it as a normal function
    cmd2 = [
        env_python,
        "./aux/cif2pdb.py",
        output_dir,
        out_dir_pdbs]
    
    #print("Running:", " ".join(map(str, cmd2)))
    result = subprocess.run(cmd2, check=True, capture_output=True, text=True)
    print(result.stdout)
    print(result.stderr)
    
    return out_dir_pdbs

def generate_protpardelle_ensembles(input_dir, n_ensembles):
    
    input_dir = Path(input_dir).resolve()
    output_dir = (input_dir / "protpardelle_ensembles").resolve()

    output_dir.mkdir(parents=True, exist_ok=True)

    pdb_list = (output_dir / "pdb_list").resolve()

    pdb_files = sorted(input_dir.glob("*.pdb"))
    with open(pdb_list, "w") as f:
        for p in pdb_files:
            f.write(p.name + "\n")

    cmd3 = [
        env_python,
        f"{caliby_dir}/caliby/eval/sampling/generate_ensembles.py",
        "--config-dir", f"{caliby_dir}/caliby/eval/sampling",
        "--config-name", "generate_ensembles",
        f"+model_params_path={model_params}",
        f"+input_cfg.pdb_dir={input_dir}",
        f"+input_cfg.pdb_name_list={pdb_list}",
        f"num_samples_per_pdb={n_ensembles}",
        f"+out_dir={output_dir}",
        "seed=0",
        "batch_size=1",
        f"++sampling_yaml_path={caliby_dir}/caliby/configs/protpardelle-1c/multichain_backbone_partial_diffusion.yaml",
        "hydra.run.dir=."]
    
    result = subprocess.run(cmd3, check=True, capture_output=True, text=True)
    print(result.stdout)
    print(result.stderr)

    return output_dir 

def generate_sequences_caliby(clean_dir_pdbs, ensembles_dir, n_seqs):

    # setting
    os.environ["HYDRA_FULL_ERROR"] = "1"

    CKPT_PATH = f"{model_params}/caliby/caliby.ckpt" # default and normal one,  there are 2 adittional version for solublee caliby still not tested
    SAMPLING_YAML = f"{caliby_dir}/caliby/configs/seq_des/atom_mpnn_inference.yaml"

    # Mirrors (only needs to be done once?)
    pdb_mirror = "./mirrors/pdb"
    ccd_mirror = "./mirrors/ccd"

    Path(pdb_mirror).mkdir(parents=True, exist_ok=True)
    Path(ccd_mirror).mkdir(parents=True, exist_ok=True)

    # Export mirrors 
    os.environ["PDB_MIRROR_PATH"] = pdb_mirror
    os.environ["CCD_MIRROR_PATH"] = ccd_mirror

    # here input_dir is "protpardelle ensembles"
    # Find the main ensemble folder (the one with cc95-...-rewind150)
    roots = glob(f"{ensembles_dir}/cc95-*rewind150")
    ensembles_cc95= roots[0] if roots else None

    seq_dir = f"{ensembles_dir}/seq_des_from_ensemble"
    list_dir = f"{seq_dir}/pdb_lists"
    Path(seq_dir).mkdir(parents=True, exist_ok=True)
    Path(list_dir).mkdir(parents=True, exist_ok=True)

    for pdb in clean_dir_pdbs.glob("*.pdb"):
      pdb_name = pdb.stem

      # we conditioned on a single ensemble (1 ensemble = 32 conformers), generated on top of the same target
      ensemble_path = Path(ensembles_cc95) / pdb_name

      if not ensemble_path.exists():
          print(f"[WARN] No ensemble found for {pdb_name}, skipping.")
          continue

      # Create individual list files to track ensembles used for designing on top of each pdb
      list_file = f"{list_dir}/{pdb_name}_info"
      with open(list_file, "w") as f:
          f.write(pdb_name + "\n")
          print(f"Single ensemble to be used to sample sequences on top of {pdb_name}: {ensemble_path}\n")
      
      # Output directory for designed sequences
      out_designs = f"{seq_dir}/{pdb_name}"
      Path(out_designs).mkdir(exist_ok=True)

      cmd4 = [
          env_python,
          f"{caliby_dir}/caliby/eval/sampling/seq_des_multi_ensemble.py",
          f"ckpt_path={CKPT_PATH}",
          f"input_cfg.conformer_dir={ensembles_cc95}", # where all the ensembles from the family members are saved
          f"input_cfg.pdb_name_list={list_file}", # selected ensembles to be used for this specific target
          f"seq_des_cfg.atom_mpnn.ckpt_path={CKPT_PATH}",
          f"seq_des_cfg.atom_mpnn.sampling_cfg={SAMPLING_YAML}",
          f"sampling_cfg_overrides.num_seqs_per_pdb={n_seqs}",
          "seed=0",
          f"out_dir={out_designs}",
          "hydra.run.dir=."
      ]

      result = subprocess.run(cmd4, check=True, capture_output=True, text=True)
      print(result.stdout)
      print(result.stderr)

    return seq_dir

def merge_csv(seqs_dir):

    # Path to save the merged CSV
    merged_csv = os.path.join(seqs_dir, "seq_des_outputs_from_all_targets.csv")

    # Merge logic
    all_dfs = []
    for root, dirs, files in os.walk(seqs_dir):
        for file in files:
            if file == "seq_des_outputs.csv":
                full_path = os.path.join(root, file)
                try:
                    df = pd.read_csv(full_path)
                    # Add PDB name (subfolder) as a column
                    pdb_name = os.path.basename(os.path.dirname(full_path))
                    df.insert(0, "PDB_ID", pdb_name)
                    all_dfs.append(df)
                    #print(f"[INFO] Loaded {file} from {pdb_name} ({len(df)} rows)")
                except Exception as e:
                    print(f"[WARNING] Could not read {full_path}: {e}")

    # Combine all DataFrames
    if all_dfs:
        merged_df = pd.concat(all_dfs, ignore_index=True)
        merged_df.to_csv(merged_csv, index=False)
        print(f"\nMerged {len(all_dfs)} CSV files -> {merged_csv}")
        print(f"   Total rows: {len(merged_df)}")
    else:
        print("No seq_des_outputs.csv files found.")
    
    return merged_csv


def csv2fasta(seqs_dir):

    # Get merged csv
    csv_file = merge_csv(seqs_dir)    
    # Path to save the FASTA
    fasta_path_full = os.path.join(seqs_dir,"seq_des_outputs_from_all_targets.fasta")

    # Read the CSV
    df = pd.read_csv(csv_file)

    # Write to multi-FASTA
    with open(fasta_path_full, "w") as fasta_file:
        for _, row in df.iterrows():
            out_pdb = row["out_pdb"]  # column with PDB/sequence identifier
            seq = row["seq"]          # column with sequence
            header = Path(out_pdb).stem  # basename without extension
            fasta_file.write(f">{header}\n{seq}\n")

    print(f"Multi-fasta (complete) saved to: {fasta_path_full}")


    # save best U (one per target)
    # Path to save the FASTA
    fasta_path_bestU = os.path.join(seqs_dir,"seq_des_best_U_per_target.fasta")


    # Select the best sequence per example_id based on column 'U'
    # Assuming lower 'U' is better
    best_df = df.loc[df.groupby("example_id")["U"].idxmin()]
    with open(fasta_path_bestU, "w") as fasta_file:
        for _, row in best_df.iterrows():
            out_pdb = row["out_pdb"]  # column with PDB/sequence identifier
            seq = row["seq"]          # column with sequence
            header = Path(out_pdb).stem  # basename without extension
            fasta_file.write(f">{header}\n{seq}\n")

    print(f"Multi-fasta (best U) saved to: {fasta_path_bestU}")   

    return fasta_path_full, fasta_path_bestU

def extract_pdbs(pdb_folder, fasta_file, output_folder):

    fasta_ids = {record.id for record in SeqIO.parse(fasta_file, "fasta")}
 
    for pdb_file in os.listdir(pdb_folder):
        if not pdb_file.endswith(".pdb"):
            continue

        # check if the PDB filename (without extension) matches any FASTA ID
        pdb_id = os.path.splitext(pdb_file)[0]
        if pdb_id in fasta_ids:
            src = os.path.join(pdb_folder, pdb_file)
            dst = os.path.join(output_folder, pdb_file)
            shutil.copy2(src, dst)
            #print(f"Copied {pdb_file} → {output_folder}")

def check_len(clean_dir_pdbs, ensembles_dir, all_models_dir): # corregir
    # logic: number of CA atoms = sequence length

    for pdb in Path(clean_dir_pdbs).glob("*.pdb"):  # for each family member
        results = []

        pdb_name = pdb.stem
        print(f"Checking {pdb_name}")
        #ensemble_pattern = f"{ensembles_dir}/cc95-*rewind150/{pdb_name}"
        #ensemble_pattern = (Path(ensembles_dir)/ "cc95-*rewind150"/ pdb_name).resolve()
        #ensemble_pattern = str(ensemble_pattern)
        
        ensembles_dir = Path(ensembles_dir).resolve()
        
        for pdb_path in ensembles_dir.glob(f"cc95-*rewind150/{pdb_name}/*.pdb"):
            length = 0

        #for pdb_path in Path().glob(f"{ensemble_pattern}/*.pdb"):  # check conformers len
            #length = 0

            with open(pdb_path, "r") as f:
                for line in f:
                    if line.startswith("ATOM") and line[12:16].strip() == "CA":
                        length += 1

            pdb_id = pdb_path.stem
            results.append((pdb_id, length))
        
        all_models_dir = Path(all_models_dir).resolve()
        for pdb_path in all_models_dir.glob(f"*{pdb_name}*.pdb"):

        #for pdb_path in Path().glob(f"{all_models_dir}/*{pdb_name}*.pdb"):  # check designs len (esmfold models)
            length = 0

            with open(pdb_path, "r") as f:
                for line in f:
                    if line.startswith("ATOM") and line[12:16].strip() == "CA":
                        length += 1

            pdb_id = pdb_path.stem
            results.append((pdb_id, length))

        # print(results) # for debug
        # check equal len for all id contained in results, for each member
        length_groups = defaultdict(list)

        for pdb_id, length in results:
            length_groups[length].append(pdb_id)

        if len(length_groups) == 1: # ==1 because we want a single unique lengths
            length = next(iter(length_groups))
            #print(f"All structures have length {length}")
        else:
            print("Length mismatch detected:")
            for length, ids in length_groups.items():
                print(f"Length {length}:")
                for i in ids:
                    print(f"  - {i}")     

# from here revisar
# 3-letter → 1-letter mapping
three_to_one = {
    "ALA":"A","ARG":"R","ASN":"N","ASP":"D","CYS":"C",
    "GLN":"Q","GLU":"E","GLY":"G","HIS":"H","ILE":"I",
    "LEU":"L","LYS":"K","MET":"M","PHE":"F","PRO":"P",
    "SER":"S","THR":"T","TRP":"W","TYR":"Y","VAL":"V"
}

# Common residue renames (modified residues → canonical) #include here as many changes as needed
residue_renames = {
    "MSE": "MET", "HIE": "HIS", "HSD": "HIS",
    "CYX": "CYS", "CY1": "CYS", "KCX": "LYS"
}

def get_pdb_sequence(pdb_path, chain=None):
    """
    Extract concatenated amino acid sequence from a PDB file,
    including renaming of common modified residues.

    Args:
        pdb_path (str or Path): Path to PDB file
        chain (str, optional): Only extract this chain; default None = all chains

    Returns:
        str: concatenated sequence
    """
    seq_dict = {}
    last_residues = {}
    concat_seq = ""

    with open(pdb_path, "r") as f:
        for line in f:
            if line.startswith("ATOM"):
                resname = line[17:20].strip()
                # apply renaming
                resname = residue_renames.get(resname, resname)
                
                chain_id = line[21].strip()
                resnum = line[22:26].strip()

                # skip if specific chain is requested
                if chain is not None and chain_id != chain:
                    continue

                # initialize tracking
                if chain_id not in last_residues:
                    last_residues[chain_id] = None

                # avoid double counting atoms of same residue
                if last_residues[chain_id] == resnum:
                    continue
                last_residues[chain_id] = resnum

                # convert 3-letter code → 1-letter, unknowns → 'X'
                aa = three_to_one.get(resname, "X")

                # append to chain sequence
                if chain_id not in seq_dict:
                    seq_dict[chain_id] = []
                seq_dict[chain_id].append(aa)

    # concatenate sequences of all chains in alphabetical order
    for chain_id in sorted(seq_dict.keys()):
        concat_seq += "".join(seq_dict[chain_id])

    return concat_seq

def process_ref(ref, new_pdb_folder, out_dir, out_file): # need to change paths
    if ref:
        # Paths
        ref_pdb = Path(new_pdb_folder) / "pdbs"/ f"{ref}.pdb" # clean pdb of the ref
        #ref_fasta = Path(new_pdb_folder) / f"{ref}_MPNN_output" / "seqs" / f"{ref}.fa" # en este caso extraerla del cristal
        # we get direcly the sequence
        ref_fasta_str = get_pdb_sequence(ref_pdb, chain=None)
        if not ref_fasta_str:
            raise FileNotFoundError(f'The selected "Ref": {ref} could not be extracted')

        # Copy PDB to out_dir
        out_dir = Path(out_dir)
        #out_dir.mkdir(parents=True, exist_ok=True)
        shutil.copy(ref_pdb, out_dir / ref_pdb.name)

        print(f"Copied PDB ref to {out_dir}")

        # Format FASTA entry with your ref variable as header
        fasta_entry = f">{ref}\n{ref_fasta_str}\n"

        # Prepend to existing file content
        if os.path.exists(out_file): # aqui cambiar por el best US
            with open(out_file, "r") as f:
                existing_content = f.read()
        else:
            existing_content = ""

        with open(out_file, "w") as f:
            f.write(fasta_entry)
            f.write(existing_content)

        print(f"Copied sequence ref to {out_file}") 

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
        
# run esmfold using different module 
def run_esmfold(dir, jobid): # better do not use

    # even if we build the module, in here we need the script directly, change path if needed
    esmfold = "/gpfs/scratch/bsc08/bsc008609/CALIBY/funfams/aux/m0d_esmfold.py"

    #to_run = m0d_esmfold.main(dir, "ESMFold", jobid)

    # important to set bash instead of just sh in the shebang
    slurm_command = f"""#!/bin/bash
#SBATCH -A
#SBATCH -q
#SBATCH --account=bsc08
#SBATCH --job-name=esmfold
#SBATCH --output=./debug/esmfold_{jobid}_%j.out
#SBATCH --error=./debug/esmfold_{jobid}_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=40
#SBATCH --gres=gpu:1
#SBATCH --qos=acc_bscls

module purge
module load anaconda cuda/11.8
source $(conda info --base)/etc/profile.d/conda.sh
source activate esmfold
    
python3 {esmfold} --pdb_folder {dir} --modelling ESMFold --jobid {jobid}
   """
    jobs = []
    slurm_file = f"submit_esmfold_{jobid}.sh"
    with open(slurm_file, "w") as f:
        f.write(slurm_command)
    
    #print("jobid sent to the ACC queue")
    result = subprocess.run(["sbatch", slurm_file], capture_output=True, text=True, check=True)
    # as this time is just one job, we do not need several dependencies
    print("Raw sbatch output:", result.stdout)
    wait_for_sbatch(result.stdout)
    
    # maybe borrar el sh de esmfold once finished? 
    # no return in this case
    #all_models_dir, out_file_esmfold, out_dir_esmfold
    #return all_models_dir, out_file_esmfold, out_dir_esmfold 


def run_esmfold_inline(dir, jobid):
    esmfold = "/gpfs/scratch/bsc08/bsc008609/CALIBY/funfams/aux/m0d_esmfold.py"
    python_bin = "/apps/ACC/ANACONDA/2023.07/envs/esmfold/bin/python"  # <- exact path from `which python` in esm env

    cmd = f"""
module unload miniforge
module load anaconda/2023.07 cuda/11.8

{python_bin} {esmfold} --pdb_folder {dir} --modelling ESMFold --jobid {jobid}
"""

    result = subprocess.run(cmd, shell=True, executable="/bin/bash",
                            capture_output=True, text=True)

    #print("===== STDOUT =====")
    #print(result.stdout)
    #print("===== STDERR =====")
    #print(result.stderr)

    if result.returncode != 0:
        raise RuntimeError(f"ESMFold failed for job {jobid}. See STDERR above.")
    

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
    
    # the one provided by miguel with msas from all families
    clean_aux_dir = "/gpfs/scratch/bsc08/bsc008092/CATH/enzymes/MIRIAM/MSA_Clean_aux_files/"
    # in here there are multiples msas look for 1 in specific
    clean_aux_file_template = os.path.join(clean_aux_dir, f"MSA_Clean_aux_{jobid}.fasta") # jobid from arguments
    if not os.path.exists(clean_aux_file_template):
        raise FileNotFoundError(f"{clean_aux_file_template} not found")
    # input_fasta is the multifasta containing best designs
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

############################

# Fix  global model variables
model_params="/gpfs/scratch/bsc08/bsc008609/CALIBY/v2/caliby-main/model_params"
caliby_dir="/gpfs/scratch/bsc08/bsc008609/CALIBY/v2/caliby-main"

# Set python version form the proper environment to avoid issues (also global variable)
env_python="/gpfs/apps/MN5/ACC/MINIFORGE/24.3.0-0/envs/caliby_fix/bin/python"


def main():
    
    args = parse_args() # Retrieve args

    print("### IMPORTANT CONSIDERATIONS ###\n" \
    "Please be sure that your input PDB folder:\n" \
    "1)only contain PDB files, not any other folder or file\n" \
    "2)contains as many PDB as sequences included in you auxiliar MSA (same headers and sequences)\n"
    "It is strongly recommended that you cancel this task and run it again once you are sure!")
    # may pending: include path to aux msa by command line
    

    if all([args.pdb_folder, args.num_ensembles_per_target, args.num_seqs_per_target, args.jobid]):

        dir_basename = basename(str(args.pdb_folder))
        print(f"TARGET FAMILY - {dir_basename}") # Retrieve name of the protein family (for funfams must coincide with args.jobid)

        print(" \n<<< SEQUENCE DESIGN BY CALIBY >>>\n")
        
        # caliby function to clean pdb cannot convert ' ?' to an integer (? present in some af models) --> so crash
        # add an additional function prior to that one within clean_input_pdbs function
        print(f"[INFO] Cleaning input structures...") # using script provided by caliby authors
        clean_dir_pdbs = clean_input_pdbs(args.pdb_folder)

        print(f"[INFO] Generating ensembles...")
        ensembles_dir = generate_protpardelle_ensembles(clean_dir_pdbs, args.num_ensembles_per_target)

        print(f"[INFO] Sampling sequences ensemble-conditioned...")
        sequences_dir = generate_sequences_caliby(clean_dir_pdbs, ensembles_dir, args.num_seqs_per_target)

        print(f"[INFO] Generating multi-fasta for all designs...") # merging csv (caliby sequences output format)
        fasta_path_full, fasta_path_bestU = csv2fasta(sequences_dir) # fasta_path_full=designs for all family members, fasta_path_bestU= best design per family member
        
        print("\n<<< STRUCTURE PREDICTION BY ESMFOLD >>>\n") # Run ESMFold in its proper env

        print(f"[INFO] Predicting designs structure with ESMFold...")
        #run_esmfold(sequences_dir, args.jobid) # maybe longer waiting time because its launching a new job
        # we include a function to run in the same job just loading new modules/env(new function implemented feb26)
        run_esmfold_inline(sequences_dir, args.jobid) # best plddt models selection included even if not nedeed
        
        print(f"[INFO] Grouping best designs per target according to U scores...")
        all_models_dir = Path(sequences_dir) / f"{args.jobid}_ESMFold"
        best_U_models_dir = Path(sequences_dir) / "Best_designs_U"
        best_U_models_dir.mkdir(parents=True, exist_ok=True)
        extract_pdbs(all_models_dir, fasta_path_bestU, best_U_models_dir)
        shutil.move(fasta_path_bestU, best_U_models_dir)

        print(f"[INFO] Checking output lengths...") # Check if designs are that longer as its native target
        check_len(clean_dir_pdbs, ensembles_dir, all_models_dir)  # (compare native.len, ensmeble.len and design.len), for each family member

    else:
        raise ValueError("IMPORTANT: You must provide a directory path conatining PDB files to run the pipeline")
    
    print("\n<<< FRUSTRAEVO >>>\n") # Run FrustraEvo in its proper env, and in GPP

    if args.ref:
    
        print(f"[INFO] Including native reference...")
        # new path for gasta because we move it previously
        fasta_path_bestU = Path(best_U_models_dir) / "seq_des_best_U_per_target.fasta"
        process_ref(args.ref, clean_dir_pdbs, best_U_models_dir, fasta_path_bestU)

        print(f"[INFO] Preparing FrustraEvo input...")
        prepare_input_FE_template_based(fasta_path_bestU, args.ref, args.jobid)
        print(" < Currently using native template-based MSA >")
        
        print(f"[INFO] Running FrustraEvo...")
        run_FE(best_U_models_dir, args.jobid, args.ref) # use designs lower U
        # Once finished, remove sh and pdf --> why its outside?
        [os.remove(f) for f in glob.glob(f"submit_FE_{ref}.sh")]
        [os.remove(f) for f in glob.glob("Rplot*.*")] # no funciona esto 
    else:
        raise ValueError("IMPORTANT: You must provide a reference to run FrustraEvo")

if __name__ == "__main__":
    start = time.perf_counter()
    main()
    end = time.perf_counter()
    print(f"Elapsed time (not counting FrustraEvo): {end - start:.6f} seconds")
