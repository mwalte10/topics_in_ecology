#!/bin/csh

#$ -N test_
#$ -t 1-25

module load bio/R/3.3.1-gcc

Rscript test.R $SGE_TASK_ID
