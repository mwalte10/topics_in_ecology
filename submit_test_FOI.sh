#!/bin/csh

#$ -N test_FOI_
#$ -t 1-20

module load bio/R/3.4.0

Rscript test_FOI.R $SGE_TASK_ID