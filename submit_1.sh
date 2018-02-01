#!/bin/csh

#$ -N vac_
#$ -t 4

module load bio/R/3.3.1-gcc

Rscript test_vac_parm.R $SGE_TASK_ID