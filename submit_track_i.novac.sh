#!/bin/csh

#$ -N beta.inf_
#$ -t 1-11

module load bio/R/3.3.1-gcc

Rscript track_i.novac.R $SGE_TASK_ID