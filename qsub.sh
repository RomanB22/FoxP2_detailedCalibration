#!/bin/bash
#$ -cwd
#$ -N BluePyOpt
#$ -q cpu.q
#$ -pe 1
#$ -l h_vmem=10G
#$ -l h_rt=3:00:00
#$ -o /ddn/rbarav/FoxP2_detailedCalibration/BBP_1.run
#$ -e /ddn/rbarav/FoxP2_detailedCalibration/BBP_1.err

python -u mainDetailed.py 
