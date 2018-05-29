#!/bin/csh

#$ -N model_sp9_
#$ -t 1-400

module load bio/R/3.3.1-gcc

Rscript model_fixed.R $SGE_TASK_ID
