#!/bin/bash
#SBATCH --job-name=cstopo
#SBATCH --array=0-31
#SBATCH --time=04:00:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=1
#SBATCH --output=logs/%A_%a.out

module load python
python run_trousers.py --N 40000 --seeds 128 \
    --task $SLURM_ARRAY_TASK_ID --ntasks 32 --out results/

# then, once the array completes:
#   python run_trousers.py --merge results/
