#!/bin/csh

#$ -N jan_22_
#$ -t 1-20

module load bio/R/3.4.0

Rscript model_jan_22.R $SGE_TASK_ID