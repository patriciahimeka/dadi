import sys
import os
import numpy as np
import dadi
import pylab
import Models_2D
import Optimize_Functions
from datetime import datetime

'''
usage: python dadi_2D_02_second_optimizations.py

Requires the Models_2D.py and Optimize_Functions.py scripts to be in same working directory. 
This is where all the population models and functions are stored for this script. 

Script will perform optimizations from multiple starting points using a
2-fold perturbed set of USER SELECTED starting values for parameters. 

Dan Portik
daniel.portik@uta.edu -> danielportik@email.arizona.edu
October 2017
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
# Now prepare to run the optimization round 2 function, which is defined in the 
# script 'Optimize_Functions.py'.

#"no_mig"
# 3 values
no_mig_params = [0.1528,	0.6922,	0.065]

#"sym_mig"
#4 Values
sym_mig_params = [0.2522,	0.8359,	1.7027,	4.9596]

#"asym_mig"
#5 Values
asym_mig_params = [0.262,	1.2659,	2.0665,	0.6688,	2.9319]

#"anc_sym_mig
# 5 values
anc_sym_mig_params= [0.4173,	0.7981, 4.1216,	8.3792,	0.0666]

# anc_asym_mig
# 6 values
anc_asym_mig_params= [0.7397,	1.778,	1.416,	0.2986,	2.5663,	0.0727]

#"sec_cont_sym_mig"
#5 values
sec_contact_sym_mig_params = [0.12942371,	0.43741949,	3.38312095,	0.51726987,	3.64316629]

#sec_cont_asym_mig
#6 values
sec_contact_asym_mig_params = [0.5035,	1.9039,	0.7924,	0.7909,	0.3214,	1.9369]

#sym_mig_two_epoch
#6 values
sym_mig_twoepoch_params = [0.4165,	1.6277,	0.376,	1.1284,	0.2854,	3.4936]

#asym_mig_two_epoch
#8 values
asym_mig_twoepoch_params = [0.161,	1.247,	0.4348,	1.6131,	3.436,	0.657,	3.8778,	0.1429]

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
#Optimize_Functions.Optimize_Round2(pts, fs, outfile, reps, maxiter, "no_mig", no_mig_params)

# Split into two populations, with continuous symmetric migration.
#Optimize_Functions.Optimize_Round2(pts, fs, outfile, reps, maxiter, "sym_mig", sym_mig_params)

# Split into two populations, with continuous asymmetric migration.
#Optimize_Functions.Optimize_Round2(pts, fs, outfile, reps, maxiter, "asym_mig", asym_mig_params)

# Ancien symmetric migration
#Optimize_Functions.Optimize_Round2(pts, fs, outfile, reps, maxiter, "anc_sym_mig", anc_sym_mig_params)

# Anceitn asym migration
#Optimize_Functions.Optimize_Round2(pts, fs, outfile, reps, maxiter, "anc_asym_mig", anc_asym_mig_params)

# Secondary contact symm. migration
Optimize_Functions.Optimize_Round2(pts, fs, outfile, reps, maxiter, "sec_contact_sym_mig", sec_contact_sym_mig_params)

# Secondary contact asym migration
#Optimize_Functions.Optimize_Round2(pts, fs, outfile, reps, maxiter, "sec_contact_asym_mig", sec_contact_asym_mig_params)

# Two epoch model sym mig
#Optimize_Functions.Optimize_Round2(pts, fs, outfile, reps, maxiter, "sym_mig_twoepoch", sym_mig_twoepoch_params)

# Two epoch asym mig
#Optimize_Functions.Optimize_Round2(pts, fs, outfile, reps, maxiter, "asym_mig_twoepoch", asym_mig_twoepoch_params)

#===========================================================================
#clock the amount of time to complete the script
t_finish = datetime.now()
elapsed = t_finish - t_begin
print '\n', '\n', "-----------------------------------------------------------------------------------------------------"
print "Finished all analyses!"
print "Total time: {0} (H:M:S)".format(elapsed)
print "-----------------------------------------------------------------------------------------------------", '\n', '\n'
#===========================================================================


