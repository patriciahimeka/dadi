import sys
import os
import numpy as np
import dadi
import pylab
import Models_2D
import Optimize_Functions
from datetime import datetime

'''
Third round of optimaztion for model fitting
Populations Daito (D) and Ryukyu (O)

usage: python 2D_03_third_optimizations_D-O.py
Requires the Models_2D.py and Optimize_Functions.py scripts to be in same working directory. 

Scripts adapted from Dan Portik: https://github.com/dportik/dadi_pipeline

Patricia Wepfer
August 2018
'''
#keep track of start time
t_begin = datetime.now()


#**************
#pop_ids is a list which should match the populations headers of your SNPs file columns
pop_ids=["Da", "Ok"]
#projection sizes, in ALLELES not individuals
proj_1 = [16,26]

#Convert this dictionary into folded AFS object
#[polarized = False] creates folded spectrum object
#fs_1 = dadi.Spectrum.from_data_dict(dd1, pop_ids=pop_ids, projections = proj_1, polarized = False)
fs_1 = dadi.Spectrum.from_file('D-O.sfs')


#print relevant info to screen
print '\n', '\n', "Data for spectrum:"
print "projection", proj_1
print "sample sizes", fs_1.sample_sizes
print "Segregating sites",fs_1.S(), '\n', '\n'

#======================================================================================

#"no_divergence"
#leave blank, no parameters
no_divergence_params = []

#"no_mig"
#3 Values
no_mig_params = [0.1431,	0.339,	0.0155]

#"sym_mig"
#4 Values
sym_mig_params = [1.0059,	1.8537,	2.2456,	2.9731]

#"asym_mig"
#5 Values
asym_mig_params = [0.5353,	1.5052,	4.985,	2.3124,	1.8503]


#======================================================================================
#Input some of the basic reusable arguments here specific to your data set

#These are specific to your data set:
#**************
#grid choice
pts = [30,40,50]
#prefix for output file naming
outfile = "Da-Oki"

#These can be left alone, unless you want more searches:
#spectrum object name (we defined this above)
fs = fs_1
#integer to control number of replicates per model
reps = int(50)
#max number of iterations per optimization step (though see dadi user group for explanation)
maxiter = int(50)


#======================================================================================
# Now call the function with the relevant arguments.

# Standard neutral model, populations never diverge
Optimize_Functions.Optimize_Round3(pts, fs, outfile, reps, maxiter, "no_divergence", no_divergence_params)

# Split into two populations, no migration.
Optimize_Functions.Optimize_Round3(pts, fs, outfile, reps, maxiter, "no_mig", no_mig_params)

# Split into two populations, with continuous symmetric migration.
Optimize_Functions.Optimize_Round3(pts, fs, outfile, reps, maxiter, "sym_mig", sym_mig_params)

# Split into two populations, with continuous asymmetric migration.
Optimize_Functions.Optimize_Round3(pts, fs, outfile, reps, maxiter, "asym_mig", asym_mig_params)


#===========================================================================
#clock the amount of time to complete the script
t_finish = datetime.now()
elapsed = t_finish - t_begin
print '\n', '\n', "-----------------------------------------------------------------------------------------------------"
print "Finished all analyses!"
print "Total time: {0} (H:M:S)".format(elapsed)
print "-----------------------------------------------------------------------------------------------------", '\n', '\n'
#===========================================================================



