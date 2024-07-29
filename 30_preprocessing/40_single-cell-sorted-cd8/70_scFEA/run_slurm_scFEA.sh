#!/bin/bash
#SBATCH --job-name=scFEA
#SBATCH --output=scFEA.out
#SBATCH --error=%x_%j.err
#SBATCH --partition=standard
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=02:00:00
module purge
eval "$(conda shell.bash hook)"
conda activate scFEA_2024
# Define directories and files
DATA_DIR="../../../../scFEA/data/"
INPUT_DIR="/data/projects/2021/MicrobialMetabolites/single-cell-sorted-cd8/results/scFEA/tmp"
RES_DIR="/data/projects/2021/MicrobialMetabolites/single-cell-sorted-cd8/results/scFEA/tmp"
TEST_FILE="matrix.csv"
MODULE_GENE_FILE="module_gene_complete_mouse_m168.csv"

# Run Python script
srun python ../../../../scFEA/src/scFEA.py  --input_dir $INPUT_DIR --res_dir $RES_DIR --test_file $TEST_FILE --moduleGene_file $MODULE_GENE_FILE
