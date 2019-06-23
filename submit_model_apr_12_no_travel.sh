#!/bin/csh

#$ -N apr_12_no_
#$ -t 1

module load bio/R/3.4.0

Rscript model_apr_12_no_travel.R $SGE_TASK_ID