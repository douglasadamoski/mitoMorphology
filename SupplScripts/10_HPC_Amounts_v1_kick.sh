#!/usr/bin/bash 
#
#
# example: nohup ./10_HPC_Amounts_v1_kick.sh "A99.28_ectopica_parkin_controle_minMitoLabel_300.RData" &
# kill -9 PID
#
echo "Filename now: $1";

# activate environment
source ~/anaconda3/etc/profile.d/conda.sh
conda activate R413

# make the run
/mnt/nfs/home/douglas/anaconda3/envs/R413/bin/Rscript --vanilla 10_HPC_Amounts_v1.R $1 > logAmounts.$1.txt

# Fecha o conda
conda deactivate

# Sair
exit 0
