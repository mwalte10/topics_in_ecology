#!/bin/csh

#$ -N jan_30_change_
#$ -t 1-20

module load bio/R/3.4.0

Rscript model_jan_30_trav.R $SGE_TASK_ID