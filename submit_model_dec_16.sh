#!/bin/csh

#$ -N sp9_
#$ -t 1-1000

module load bio/R/3.3.1-gcc

Rscript model_dec_16.R $SGE_TASK_ID