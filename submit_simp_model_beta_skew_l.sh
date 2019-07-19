#!/bin/csh

#$ -N beta_
#$ -t 1-100

module load bio/R/3.4.0

Rscript simp_model_beta_skew_l.R $SGE_TASK_ID