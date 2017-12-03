#!/bin/csh

#$ -N dengvaxia_vac
#$ -t 1-50

module load bio/R/3.3.1-gcc

Rscript basic_vac.R $SGE_TASK_ID