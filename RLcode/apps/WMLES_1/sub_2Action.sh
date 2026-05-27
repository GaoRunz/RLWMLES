#!/bin/bash
#SBATCH -p amd_512
#SBATCH -N 1
#SBATCH -n 2
#SBATCH --mem=24G

source /public3/soft/modules/module.sh
module load cmake/3.24.1 anaconda/2023.07-2-hxl gcc/8.1.0-wjl mpich/3.1.4-gcc8.1.0
source ~/.bash_profile
smarties.py -r zzz_train_aNew_LES1_2w_02_bad_action2_fre5

