#!/bin/csh

#$ -N feb_18_
#$ -t 1

module load bio/R/3.4.0

Rscript model_feb_18.R $SGE_TASK_ID