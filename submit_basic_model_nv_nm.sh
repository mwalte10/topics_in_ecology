#!/bin/csh

#$ -N beta.rf_
#$ -t 1-51

module load bio/R/3.3.1-gcc

Rscript basic_model_nv_nm.R $SGE_TASK_ID