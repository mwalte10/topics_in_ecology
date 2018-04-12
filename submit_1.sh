#!/bin/csh

#$ -N vac_
#$ -t 4

module load bio/R/3.3.1-gcc

<<<<<<< HEAD
Rscript test_vac_parm.R $SGE_TASK_ID
=======
Rscript big_vac.R $SGE_TASK_ID
>>>>>>> 2e0b4bf207a64e845b3982dc3803073089900541
