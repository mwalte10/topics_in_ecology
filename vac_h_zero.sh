#!/bin/csh

#$ -N dengvaxia_final
#$ -t 1-50

module load bio/R/3.3.1-gcc

Rscript vac_mov_real.R $SGE_TASK_ID