#!/bin/csh

#$ -N sp9_
#$ -t 1-8000

module load bio/R/3.3.1-gcc

Rscript model_native_parm.R $SGE_TASK_ID