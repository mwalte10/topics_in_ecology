#!/bin/csh

#$ -N beta.inf_
#$ -t 1

module load bio/R/3.3.1-gcc

Rscript testing.R $SGE_TASK_ID