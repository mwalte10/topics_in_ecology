#!/bin/csh

#$ -N june_20_
#$ -t 1-20

module load bio/R/3.4.0

Rscript model_june_20_travel.R $SGE_TASK_ID