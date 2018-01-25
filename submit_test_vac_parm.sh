#!/bin/csh

#$ -N vac_parm_test
#$ -t 1-4

module load bio/R/3.3.1-gcc

Rscript test_vac_parm.R $SGE_TASK_ID