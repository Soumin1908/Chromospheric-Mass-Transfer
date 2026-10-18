#!/usr/bin/env python
#SBATCH --time=20:00:00 # walltime
#SBATCH --nodes=2 # number of nodes
#SBATCH --exclusive
#SBATCH -J "rp2000" # job name
#SBATCH --output=output.txt
print('hello world!')
import os 
path = '/resnick/groups/elbadry_lab/sbhatta2/mesa-r24.03.1/binary/chromosphere_paper/chromosphere_accretion_ritter_p2000/'
os.chdir(path)
os.system('./rn')