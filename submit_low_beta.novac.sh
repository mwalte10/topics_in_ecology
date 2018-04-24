#!/bin/csh

#$ -N low_beta.novac.rf
#$ -t 1-101

module load bio/R/3.3.1-gcc

Rscript low_beta.novac.R $SGE_TASK_ID
