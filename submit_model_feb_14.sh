#!/bin/csh

#$ -N feb_14_
#$ -t 1-20

module load bio/R/3.4.0

Rscript model_feb_14.R $SGE_TASK_ID