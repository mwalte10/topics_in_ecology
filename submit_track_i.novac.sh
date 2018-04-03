#!/bin/csh

#$ -N beta.nv.inf_
#$ -t 1

module load bio/R/3.3.1-gcc

Rscript track_i.novac.R $SGE_TASK_ID