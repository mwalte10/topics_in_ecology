#!/bin/csh

#$ -N test_sp9_
#$ -t 1-400

module load bio/R/3.3.1-gcc

Rscript model_nov_6.R $SGE_TASK_ID