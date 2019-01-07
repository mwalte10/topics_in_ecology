#!/bin/csh

#$ -M mwalte10@nd.edu
#$ -m abe
#$ -N dec_23_
#$ -t 390-400

module load bio/R/3.4.0

Rscript model_dec_23.R $SGE_TASK_ID