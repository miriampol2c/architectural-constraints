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
    parser = argparse.ArgumentParser(description="Preprocess input PDBs before design pipeline.")
    parser.add_argument("--pdb_folder", type=str, required=True, help="Absolute path to the folder containing .pdb files")

    return parser.parse_args()

####################################################
# FUNCTIONS TO FILTER AND CURATE PROTEIN BACKBONES #
####################################################

def check_backbone_complete(pdb_input):
    """
    Check if the backbone atoms (N, CA, C, O, CB) are complete.
    Accepts either a file path or a BioPython Structure object.
    """
    if isinstance(pdb_input, str):
        parser = PDBParser(QUIET=True)
        structure = parser.get_structure("pdb", pdb_input)
    else:
        structure = pdb_input  # Assume it's already a Structure object

    required_atoms = {"N", "CA", "O", "CB"}
    residue_atoms = defaultdict(set)
    glycine_residues = set()

    for model in structure:
        for chain in model:
            for residue in chain:
                if not residue.id[0] == " ":
                    continue  # Skip heteroatoms and water

                res_id = (chain.id, residue.id[1], residue.id[2])
                resname = residue.get_resname()
                if resname == "GLY":
                    glycine_residues.add(res_id)

                for atom in residue:
                    name = atom.get_name()
                    if name in required_atoms or name == "C":
                        residue_atoms[res_id].add(name)

    for res_id, atoms in residue_atoms.items():
        if not {"N", "CA", "O", "C"}.issubset(atoms):
            return False
        if res_id not in glycine_residues and "CB" not in atoms:
            return False

    return True

def complete_backbone_with_modeller(pdb_path):
    #modeller_script_path = '/gpfs/scratch/bsc08/bsc008609/CSA_enzymes/MissingAtoms.py'
    
    modeller_script_path = os.path.join(os.getcwd(), 'ext_scripts', 'MissingAtoms.py')
    # Use the direct path to the python binary inside your environment
    modeller_python = '/apps/GPP/ANACONDA/2023.07/envs/frustratometer-reqs/bin/python'

    try:
        print("Running Modeller with direct Python call")
        subprocess.run([
            modeller_python, modeller_script_path, pdb_path
        ], check=True)
        #return pdb_path

    except subprocess.CalledProcessError as e:
        print(f"Modeller failed: {e}")
        #return None

def filter_pdb_structure(pdb_path):
    """
    Applies residue name and HETATM filtering.
    Overwrites the input PDB only if modifications were made.

    Returns:
        True if structure was modified and saved,
        False if no modifications were needed,
        None if parsing failed.
    """
    residue_renames = {
        "MSE": "MET", "HIE": "HIS", "HSD": "HIS",
        "CYX": "CYS", "CY1": "CYS", "KCX": "LYS"
    }

    standard_residues = {
        "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY",
        "HIS", "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER",
        "THR", "TRP", "TYR", "VAL"
    }

    class ProteinResidueSelector(Select):
        def accept_residue(self, residue):
            return residue.get_resname() in standard_residues

    parser = PDBParser(QUIET=True)
    io = PDBIO()

    try:
        structure = parser.get_structure("model", pdb_path)
        changes_made = False

        for model in structure:
            for chain in model:
                if not chain.id.strip():
                    chain.id = "A"
                for residue in chain:
                    resname = residue.get_resname()
                    if resname in residue_renames:
                        residue.resname = residue_renames[resname]
                        changes_made = True

                    # Normalize HETATM to ATOM
                    hetflag, resseq, icode = residue.id
                    if hetflag != " " and resname in standard_residues:
                        residue.id = (" ", resseq, icode)
                        changes_made = True

        if changes_made:
            io.set_structure(structure)
            io.save(pdb_path, select=ProteinResidueSelector())
            print(f"Filtered and saved: {os.path.basename(pdb_path)}")
            return True
        else:
            print(f"No filtering needed: {os.path.basename(pdb_path)}")
            return False

    except Exception as e:
        print(f"Failed to process {os.path.basename(pdb_path)}: {e}")
        return None

def clean_pdbs(input_folder, output_path):
        
    if not os.path.exists(output_path):
        os.makedirs(output_path)
    
    pdb_files = [f for f in os.listdir(input_folder) if f.lower().endswith(".pdb")]

    for fname in pdb_files:
        input_pdb = os.path.join(input_folder, fname)
        output_pdb = os.path.join(output_path, fname)
        # print(output_pdb) is a path

        # Step 1: Copy original to output folder
        shutil.copy2(input_pdb, output_pdb)
        print(f"Copied {fname} to {output_path}")
        
        # Step 2: Filter residues and fix atoms
        filter_pdb_structure(output_pdb)

        # Step 3: Complete missing backbone atoms with Modeller if needed (only atoms, no full residues)
        if check_backbone_complete(output_pdb): # false by default i think?
            print(f"{fname} backbone is incomplete (atoms missing)")
            print(f"Call Modeller to complete ...")
            complete_backbone_with_modeller(output_pdb) # gte model object, but the fuction save the pdb
            print(f"Done")
        else:
            print(f"{fname} backbone is completed (all atoms already included)")
    
    return output_path


# MAIN

def main(pdb_folder):

    # Filter and clean input structures, define new input folder to save them
    # Always tracking absolute path
    pdb_folder = Path(pdb_folder).resolve()
    new_pdb_folder = pdb_folder.parent / (pdb_folder.name + "_filtered")
    os.makedirs(new_pdb_folder, exist_ok=True)
    os.chmod(new_pdb_folder, 0o777)  # Temporal: Everyone: read, write, execute (not recommended unless necessary)
    # Note:
    # Filtered pdbs will be add to the new input folder (working directory)
    # Pdbs with no filtering needed will be copy directly
    new_dir = clean_pdbs(pdb_folder, new_pdb_folder)

    return new_dir
    
if __name__ == "__main__":

    # Retrieve args
    args = parse_args()

    # Name args
    pdb_folder = args.pdb_folder
    
    main(pdb_folder)
