#!/bin/csh

#$ -N mar_20_travel_
#$ -t 1-20

module load bio/R/3.4.0

Rscript model_mar_20_travel.R $SGE_TASK_ID