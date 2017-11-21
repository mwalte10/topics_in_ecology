#!/bin/csh

#$ -N dengvaxia_basic_ses
#$ -t 1-100

module load bio/R/3.3.1-gcc

Rscript basic_vaccination_movement.R $SGE_TASK_ID