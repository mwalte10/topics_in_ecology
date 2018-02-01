#!/bin/csh

#$ -N vac_
#$ -t 4

module load bio/R/3.3.1-gcc

Rscript big_vac.R $SGE_TASK_ID