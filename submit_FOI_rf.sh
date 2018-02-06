#!/bin/csh

#$ -N FOI.rf_
#$ -t 1-101

module load bio/R/3.3.1-gcc

Rscript FOI_random_forest.R $SGE_TASK_ID