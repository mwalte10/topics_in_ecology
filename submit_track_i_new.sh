#!/bin/csh

#$ -N beta.inf_
#$ -t 1-400

module load bio/R/3.3.1-gcc

Rscript track_i_new.R $SGE_TASK_ID
