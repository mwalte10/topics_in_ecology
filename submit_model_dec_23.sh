#!/bin/csh

#$ -N dec_23_
#$ -t 1-20

module load bio/R/3.4.0

Rscript model_dec_23.R $SGE_TASK_ID