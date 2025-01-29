#!/bin/bash
#SBATCH --job-name=proxies
#SBATCH --partition=mrcieu
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --time=12:00:00
#SBATCH --mem=15GB
#SBATCH --account=sscm013522

module add languages/r/4.3.3
R CMD BATCH scripts/autoimmune_mab_pregnancy/autoimmune/2_a_sel_proxies.R scripts/autoimmune_mab_pregnancy/autoimmune/2_a_sel_proxies.Rout