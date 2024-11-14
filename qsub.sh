#!/bin/bash
#$ -cwd
#$ -N BluePyOpt
#$ -q cpu.q
#$ -l h_vmem=10G
#$ -l h_rt=10:00:00
#$ -o /ddn/rbarav/FoxP2_detailedCalibration/BBP_1.run
#$ -e /ddn/rbarav/FoxP2_detailedCalibration/BBP_1.err

python -u mainDetailed.py 
