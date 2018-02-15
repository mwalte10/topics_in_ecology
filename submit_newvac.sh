#!/bin/csh

#$ -N newvac_
#$ -t 1-6

module load bio/R/3.3.1-gcc

Rscript new.vaccine.R $SGE_TASK_ID