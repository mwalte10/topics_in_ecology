#!/bin/csh

#$ -N model_vac_
#$ -t 1-25

module load bio/R/3.3.1-gcc

Rscript model_travel_null.R $SGE_TASK_ID
