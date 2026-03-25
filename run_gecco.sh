#!/bin/bash
#SBATCH -p amd_512
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 32

# module load miniforge/24.11
# mamba create -n gecco_env python=3.9 -c conda-forge -c bioconda gecco=0.9.8 biopython=1.79 -y

module load miniforge/24.11
source activate /public3/home/scg4618/.local/share/mamba/envs/gecco_env

gecco run \
    -g Cs_JDA12.genome.fa \
    --proteins C_siamense_JDA_12.fa \
    -o Cs_JDA12.gecco_output \
    -m 0.5 \
    -j 16 \
    --antismash-sideload ./Cs_JDA12.antismash8.output/Cs_JDA12.genome.json

python3 generate_gecco_report.py \
  --gff ./C_siamense_JDA_12.gff3 \
  --gff-feature gene \
  --gff-id-attr ID \
  --gff-min-overlap 0.5
