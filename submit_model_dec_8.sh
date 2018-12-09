#!/bin/csh

#$ -N output.new_
#$ -t 1-5

module load bio/R/3.3.1-gcc

Rscript model_dec_8.R $SGE_TASK_ID