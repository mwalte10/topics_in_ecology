#!/bin/csh

#$ -N movement_sp9_
#$ -t 1

module load bio/R/3.3.1-gcc

Rscript movement.R $SGE_TASK_ID
