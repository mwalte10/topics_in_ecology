#!/bin/csh

#$ -N model_sens_
#$ -t 1-400

module load bio/R/3.3.1-gcc

Rscript model_sens.R $SGE_TASK_ID
