#!/bin/bash
#SBATCH -o job_%j_%N.out
#SBATCH --partition=gpu
#SBATCH -J hw3
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --cpus-per-task=1
#SBATCH --gres=gpu:1
#SBATCH -t 01:30:00

export OMP_NUM_THREADS=4
./HeatDifu
