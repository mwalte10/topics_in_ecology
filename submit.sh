#!/bin/csh

#$ -N dengvaxia_ses
#$ -t 1-400

module load bio/R/3.3.1-gcc

Rscript script.R $SGE_TASK_ID