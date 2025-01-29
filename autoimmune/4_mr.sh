#!/bin/bash
#SBATCH --job-name=4_mr
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --time=3:00:00
#SBATCH --mem=5GB
#SBATCH --account=sscm013522

# module add languages/R/4.3.3
R CMD BATCH scripts/autoimmune_mab_pregnancy/autoimmune/4_mr.R scripts/autoimmune_mab_pregnancy/autoimmune/4_mr.Rout