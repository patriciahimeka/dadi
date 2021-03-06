import sys
import os
import numpy as np
import dadi
import pylab
import Models_2D
import Optimize_Functions
from datetime import datetime

'''
Run first round of parameter optimizations for populations 'Daito' (D) and 'Ryukyu' (O).

usage: python 01_first_optimizations-D-O.py
Requires the Models_2D.py and Optimize_Functions.py scripts to be in same working directory. 

Adapted from Dan Protik: https://github.com/dportik/dadi_pipeline

Patricia Wepfer
August 2018

'''
#keep track of start time
t_begin = datetime.now()

#===========================================================================
# Load data and display basic infos

#**************
#pop_ids is a list which should match the populations headers of your SNPs file columns
pop_ids=["D", "O"]
#projection sizes, in ALLELES not individuals
proj_1 = [16,26]

fs_1 = dadi.Spectrum.from_file('D-O.sfs')

#grid choice
pts = [30,40,50]
#prefix for output file naming
outfile = "D-O"
#**************

print '\n', '\n', "Data for spectrum:"
print pop_ids , "projection", proj_1
print "sample sizes", fs_1.sample_sizes
print "Segregating sites",fs_1.S(), '\n', '\n'

#======================================================================================
# Prepare settings to run the optimization round 2 function defined in the 
# script 'Optimize_Functions.py'.

#**************
#spectrum object name (as defined above)
fs = fs_1
#integer to control number of replicates per model
reps = int(30)
#max number of iterations per optimization step (though see dadi user group for explanation)
maxiter = int(20)

#======================================================================================
# Call the function with the relevant arguments.

# Divergence with no migration
Optimize_Functions.Optimize_Round1(pts, fs, outfile, reps, maxiter, "no_mig", no_mig_params)

# Split into two populations, with continuous symmetric migration.
Optimize_Functions.Optimize_Round1(pts, fs, outfile, reps, maxiter, "sym_mig", sym_mig_params)

# Split into two populations, with continuous asymmetric migration.
Optimize_Functions.Optimize_Round1(pts, fs, outfile, reps, maxiter, "asym_mig", asym_mig_params)

#===========================================================================
#clock the amount of time to complete the script
t_finish = datetime.now()
elapsed = t_finish - t_begin
print '\n', '\n', "-----------------------------------------------------------------------------------------------------"
print "Finished all analyses!"
print "Total time: {0} (H:M:S)".format(elapsed)
print "-----------------------------------------------------------------------------------------------------", '\n', '\n'
#===========================================================================


