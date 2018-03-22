#!/bin/csh

#$ -N beta.inf_
#$ -t 1-6

module load bio/R/3.3.1-gcc

Rscript track_i.R $SGE_TASK_ID