#!/bin/csh

#$ -N apr_10_FOI_
#$ -t 1-20

module load bio/R/3.4.0

Rscript model_apr_10_FOI.R $SGE_TASK_ID