#!/bin/csh

#$ -N time_dependent_vac
#$ -t 1-4

module load bio/R/3.3.1-gcc

Rscript time_dependent_vac.R $SGE_TASK_ID