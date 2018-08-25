#!/bin/csh

#$ -N test_
#$ -t 1-400

module load bio/R/3.3.1-gcc

Rscript model_age_dist.R $SGE_TASK_ID