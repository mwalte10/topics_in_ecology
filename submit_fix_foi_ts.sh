#!/bin/csh

#$ -N feb_18_
#$ -t 1-20

module load bio/R/3.4.0

Rscript fix_foi_ts.R $SGE_TASK_ID