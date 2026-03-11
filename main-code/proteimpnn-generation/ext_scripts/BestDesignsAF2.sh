#!/bin/bash
#SBATCH -A
#SBATCH -q
#SBATCH --account=bsc08
#SBATCH --job-name=BestDesignsAF2
#SBATCH --output=BestDesignsAF2_%j.out
#SBATCH --error=BestDesignsAF2_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --qos=acc_bscls

# Relocate logs
LOGDIR="$PWD/debug/AF2"
mkdir -p "$LOGDIR"

mv -f BestDesignsAF2_${SLURM_JOB_ID}.out "$LOGDIR/BestDesignsAF2_${SLURM_JOB_ID}.out"
mv -f BestDesignsAF2_${SLURM_JOB_ID}.err "$LOGDIR/BestDesignsAF2_${SLURM_JOB_ID}.err"

module purge
module load anaconda cuda/11.8 && source activate esmfold

# Given a filtered pdb folder (main .py script)
#new_pdb_folder="$1"

# new pdb folder has been defined and exported in the original py script

echo "Running BestDesignsAF2.py with PDB folder: ${new_pdb_folder}"

python3 /gpfs/scratch/bsc08/bsc008609/CSA_enzymes/FrustMPNN/ext_scripts/BestDesignsAF2.py --pdb_folder "${new_pdb_folder}"


