#!/bin/csh

#$ -N movement_
#$ -t 1-400

module load bio/R/3.3.1-gcc

Rscript movement.R $SGE_TASK_ID
