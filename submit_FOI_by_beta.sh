#!/bin/csh

#$ -N FOI_by_beta_
#$ -t 1-20

module load bio/R/3.4.0

Rscript FOI_by_beta.R $SGE_TASK_ID