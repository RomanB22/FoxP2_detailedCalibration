#!/bin/bash
#$ -cwd
#$ -N BluePyOpt
#$ -q cpu.q
#$ -l h_vmem=5G
#$ -l h_rt=30:00:00
#$ -o /ddn/rbarav/FoxP2_detailedCalibration/BBP_0.run
#$ -e /ddn/rbarav/FoxP2_detailedCalibration/BBP_0.err

python -u mainDetailed.py 
