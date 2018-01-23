#!/bin/csh

#$ -N test_vac
#$ -t 1-10

module load bio/R/3.3.1-gcc

Rscript test_vac.R $SGE_TASK_ID
