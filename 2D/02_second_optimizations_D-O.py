import sys
import os
import numpy as np
import dadi
import pylab
import Models_2D
import Optimize_Functions
from datetime import datetime

'''
Run second round of parameter optimizations for populations 'Daito' (D) and 'Ryukyu' (O).

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
pop_ids=["Daito", "Ryukyu"]
#projection sizes, in ALLELES not individuals
proj_1 = [16,26]

fs_1 = dadi.Spectrum.from_file('D-O.sfs')

#grid choice
pts = [30,40,50]
#prefix for output file naming
outfile = "Da-Ok"
#**************


print '\n', '\n', "Data for spectrum:"
print pop_ids , "projection", proj_1
print "sample sizes", fs_1.sample_sizes
print "Segregating sites",fs_1.S(), '\n', '\n'

#======================================================================================
# Give parameters from round 1

#"no_mig"
#3 Values
no_mig_params = [0.4783,	0.675,	0.042]

#"sym_mig"
#4 Values
sym_mig_params = [0.4354,	0.7863,	5.2871,	0.6828]

#"asym_mig"
#5 Values
asym_mig_params = [0.4587,	1.3089,	4.5929,	4.1985,	1.0622]

#"anc_sym_mig"
#5 Values
anc_sym_mig_params = [1.5358,	1.5817,	2.185,	1.7732,	0.0377]


#======================================================================================
#These can be left alone, unless you want more searches:
#spectrum object name (we defined this above)
fs = fs_1
#integer to control number of replicates per model
reps = int(30)
#max number of iterations per optimization step (though see dadi user group for explanation)
maxiter = int(30)


#======================================================================================
# Now call the function with the relevant arguments.

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


