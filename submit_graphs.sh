#!/bin/csh

#$ -N graphs_
#$ -t 1-6

module load bio/R/3.3.1-gcc

Rscript test_vac_parm.R $SGE_TASK_ID