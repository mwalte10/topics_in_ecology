#!/bin/csh

#$ -N simp_
#$ -t 1-120

module load bio/R/3.4.0

Rscript simp_model_7_7_spec_analysis.R $SGE_TASK_ID