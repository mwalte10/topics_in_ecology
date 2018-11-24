#!/bin/csh

#$ -N test_
#$ -t 1-2020

module load bio/R/3.3.1-gcc

Rscript model_nov_24.R $SGE_TASK_ID