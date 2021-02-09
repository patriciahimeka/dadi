import sys
import os
import numpy as np
import dadi
import pylab
import Models_2D
import Optimize_Functions
from datetime import datetime

'''
Run second round of parameter optimizations for populations 'Ogasawara' (C) and Ryukyu (O).

usage: python 02_second_optimizations-C-D.py
Requires the Models_2D.py and Optimize_Functions.py scripts to be in same working directory. 

Adapted from Dan Protik: https://github.com/dportik/dadi_pipeline

Patricia Wepfer
August 2018
'''
#keep track of start time
t_begin = datetime.now()

#===========================================================================
#**************
#pop_ids is a list which should match the populations headers of your SNPs file columns
pop_ids=["Ogasawara", "Ryukyu"]
#projection sizes, in ALLELES not individuals
proj_1 = [26,26]

fs_1 = dadi.Spectrum.from_file('C-O.sfs')

#grid choice
pts = [30,40,50]
#prefix for output file naming
outfile = "Ch-Ok"
#**************


print '\n', '\n', "Data for spectrum:"
print pop_ids , "projection", proj_1
print "sample sizes", fs_1.sample_sizes
print "Segregating sites",fs_1.S(), '\n', '\n'

#======================================================================================
# Give fitted parameters from Round 1:

#"no_mig"
# 3 values
no_mig_params = [0.1528,	0.6922,	0.065]

#"sym_mig"
#4 Values
sym_mig_params = [0.2522,	0.8359,	1.7027,	4.9596]

#"asym_mig"
#5 Values
asym_mig_params = [0.262,	1.2659,	2.0665,	0.6688,	2.9319]


#======================================================================================
# specify search 
fs = fs_1
#integer to control number of replicates per model
reps = int(30)
#max number of iterations per optimization step (though see dadi user group for explanation)
maxiter = int(30)


#======================================================================================
# Call functions

# Divergence with no migration
Optimize_Functions.Optimize_Round2(pts, fs, outfile, reps, maxiter, "no_mig", no_mig_params)

# Split into two populations, with continuous symmetric migration.
Optimize_Functions.Optimize_Round2(pts, fs, outfile, reps, maxiter, "sym_mig", sym_mig_params)

# Split into two populations, with continuous asymmetric migration.
Optimize_Functions.Optimize_Round2(pts, fs, outfile, reps, maxiter, "asym_mig", asym_mig_params)

#===========================================================================
#clock the amount of time to complete the script
t_finish = datetime.now()
elapsed = t_finish - t_begin
print '\n', '\n', "-----------------------------------------------------------------------------------------------------"
print "Finished all analyses!"
print "Total time: {0} (H:M:S)".format(elapsed)
print "-----------------------------------------------------------------------------------------------------", '\n', '\n'
#===========================================================================


