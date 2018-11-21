#!/bin/csh

#$ -N sp9_
#$ -t 1-400

module load bio/R/3.3.1-gcc

Rscript model_validation.R $SGE_TASK_ID