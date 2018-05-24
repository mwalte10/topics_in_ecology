#!/bin/csh

#$ -N model_travel_null_
#$ -t 1-400

module load bio/R/3.3.1-gcc

Rscript model_travel_null.R $SGE_TASK_ID
