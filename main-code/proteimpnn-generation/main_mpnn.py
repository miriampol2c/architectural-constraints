import argparse
import os
import time
from m0dules import m0d0_preprocessing, m0d1_mpnn, m0d2_esmfold, m0d2_af, m0d3_frustraevo

def parse_args():
    parser = argparse.ArgumentParser(description="Run FrustMPNN pipeline on a folder of PDB files.")
    parser.add_argument("--pdb_folder", type=str, required=True, help="Absolute path to the folder containing .pdb files")
    parser.add_argument("--num_seq_per_target", type=int, required=False, default=30, help="Number of sequences to be designed per target")
    parser.add_argument("--sampling_temp", type=float, required=False, default=0.1, help="Sampling temperature")
    parser.add_argument("--modelling", type=str, required=False, default="AF2", help="Structure prediction method: AF2, AF3, or ESMFold")
    parser.add_argument("--jobid", type=str, required=True, help="Name to track the FrustraEvo results")
    parser.add_argument("--ref", type=str, required=False, help="Choose the protein reference for the FrustraEvo analysis") # if no native id protein provided, default would be taking the first designed protein of the alignment once alph sorted
    return parser.parse_args()

def main():

    # Retrieve args
    args = parse_args()
    
    # New dir contains filtered pdbs
    #new_dir = m0d0_preprocessing.main(args.pdb_folder) # not needed if data is curated
    # if preprocessing, then we will use new_dir
    
    # Start design using the filtered pdbs
    # be sure you include as many family members as included in the native msa
    m0d1_mpnn.main(args.pdb_folder, args.num_seq_per_target, args.sampling_temp)
    
    # Predict the structures and select best designs
    if args.modelling == "ESMFold":
    	out_file, out_dir = m0d2_esmfold.main(args.pdb_folder, args.modelling)
    	#print(out_file), #print(out_dir)
    	# in here return out_file and out_dir, where the final designs are saved, needd for mod3
    
    if args.modelling == "AF2":
    	out_file, out_dir = m0d2_af.main(args.pdb_folder, args.modelling)
    
    if args.modelling == "AF3":
    	out_file, out_dir = m0d2_af.main(args.pdb_folder, args.modelling)
    
    # Prepare final input and run Frustraevo
    m0d3_frustraevo.main(args.pdb_folder, out_dir, out_file, args.jobid, args.ref, args.modelling)
    
if __name__ == "__main__":
    start = time.perf_counter()
    main()
    end = time.perf_counter()
    print(f"\nElapsed time (not counting FrustraEvo): {end - start:.6f} seconds")

