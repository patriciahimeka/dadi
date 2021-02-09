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
# Now prepare to run the optimization round 2 function, which is defined in the 
# script 'Optimize_Functions.py'.

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

#"anc_asym_mig"
#6 Values
anc_asym_mig_params = [1.2483,	3.9629,	4.518,	0.5015,	4.1635,	0.0882]

#"sec_contact_sym_mig"
#5 Values
sec_contact_sym_mig_params = [2.6136,	4.1536,	1.0742,	1.2434,	7.1761]

#"sec_contact_asym_mig"
#6 Values
sec_contact_asym_mig_params = [0.3176,	1.7005,	10.9231,	1.353,	1.7466,	0.3526]

#sym_mig_two_epoch
#6	values
sym_mig_twoepoch_params =	[2.5245,	4.4661,	1.7164,	0.8577,	0.1647,	8.2925]

#asym_mig_two_epoch
#8	values
asym_mig_twoepoch_params =[1.1702,	2.2964,	3.933,	3.1242,	2.5145,	1.5807,	2.1291,	2.7903]

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

# Ancient symmetric migration
Optimize_Functions.Optimize_Round2(pts, fs, outfile, reps, maxiter, "anc_sym_mig", anc_sym_mig_params)

#Ancient asymetric migration
Optimize_Functions.Optimize_Round2(pts, fs, outfile, reps, maxiter, "anc_asym_mig", anc_asym_mig_params)

# Secondary contact symm. migration
Optimize_Functions.Optimize_Round2(pts, fs, outfile, reps, maxiter, "sec_contact_sym_mig", sec_contact_sym_mig_params)

# Secondary contact asym migration
Optimize_Functions.Optimize_Round2(pts, fs, outfile, reps, maxiter, "sec_contact_asym_mig", sec_contact_asym_mig_params)

# Two epoch model sym mig
Optimize_Functions.Optimize_Round2(pts, fs, outfile, reps, maxiter, "sym_mig_twoepoch", sym_mig_twoepoch_params)

# Two epoch asym mig
Optimize_Functions.Optimize_Round2(pts, fs, outfile, reps, maxiter, "asym_mig_twoepoch", asym_mig_twoepoch_params)

#===========================================================================
#clock the amount of time to complete the script
t_finish = datetime.now()
elapsed = t_finish - t_begin
print '\n', '\n', "-----------------------------------------------------------------------------------------------------"
print "Finished all analyses!"
print "Total time: {0} (H:M:S)".format(elapsed)
print "-----------------------------------------------------------------------------------------------------", '\n', '\n'
#===========================================================================


