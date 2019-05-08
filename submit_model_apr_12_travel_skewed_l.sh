#!/bin/csh

#$ -N apr_12_
#$ -t 1-20

module load bio/R/3.4.0

Rscript model_apr_12_travel_skewed_l.R $SGE_TASK_ID