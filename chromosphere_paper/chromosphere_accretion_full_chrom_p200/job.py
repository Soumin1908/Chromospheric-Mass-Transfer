#!/usr/bin/env python
#SBATCH --time=20:00:00 # walltime
#SBATCH --nodes=2 # number of nodes
#SBATCH --exclusive
#SBATCH -J "fcp200" # job name
#SBATCH --output=output.txt
print('hello world!')
import os 
path = '/resnick/groups/elbadry_lab/sbhatta2/mesa-r24.03.1/binary/chromosphere_paper/chromosphere_accretion_full_chrom_p200/'
os.chdir(path)
os.system('./rn')