#!/bin/csh

#$ -N low_beta_vac
#$ -t 1-100

module load bio/R/3.3.1-gcc

Rscript low_beta_vac.R $SGE_TASK_ID
