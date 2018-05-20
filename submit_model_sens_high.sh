#!/bin/csh

#$ -N high_pop_
#$ -t 1-400

module load bio/R/3.3.1-gcc

Rscript model_sens_high.R $SGE_TASK_ID
