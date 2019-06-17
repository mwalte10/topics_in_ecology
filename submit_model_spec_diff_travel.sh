#!/bin/csh

#$ -N spec_
#$ -t 1-400

module load bio/R/3.4.0

Rscript model_spec_diff_travel.R $SGE_TASK_ID