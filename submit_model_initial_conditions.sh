#!/bin/csh

#$ -N initial_conditions_
#$ -t 1-20

module load bio/R/3.3.1-gcc

Rscript model_inital_conditions.R $SGE_TASK_ID