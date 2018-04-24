#!/bin/csh

#$ -N low_beta.novac
#$ -t 1-100

module load bio/R/3.3.1-gcc

Rscript low_beta.novac.R $SGE_TASK_ID
