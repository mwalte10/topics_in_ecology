#!/bin/csh

#$ -N mar_20_notravel_
#$ -t 1-20

module load bio/R/3.4.0

Rscript model_mar_20_notravel.R $SGE_TASK_ID