#!/bin/csh

#$ -N dengvaxia_rainbow
#$ -t 1-3

module load bio/R/3.3.1-gcc

Rscript basic_model_movement.R $SGE_TASK_ID