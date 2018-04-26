#!/bin/csh

#$ -N beta.inf_
#$ -t 1

module load bio/R/3.3.1-gcc

Rscript diff_vac.R $SGE_TASK_ID