import sys
import os
import numpy as np
import dadi
import pylab
import Models_2D
import Optimize_Functions
from datetime import datetime

'''
usage: python dadi_2D_01_first_optimizations.py

Requires the Models_2D.py and Optimize_Functions.py scripts to be in same working directory. 
This is where all the population models and functions are stored for this script. 

Dan Portik
daniel.portik@uta.edu -> danielportik@email.arizona.edu
October 2017
'''
#keep track of start time
t_begin = datetime.now()


#===========================================================================

#**************
#pop_ids is a list which should match the populations headers of your SNPs file columns
pop_ids=["C", "D"]
#projection sizes, in ALLELES not individuals
proj_1 = [26,16]

fs_1=dadi.Spectrum.from_file('C-D.sfs')

#grid choice
pts = [30,40,50]
#prefix for output file naming
outfile = "C-D"

#These can be left alone, unless you want more searches:
#spectrum object name (we defined this above)
fs = fs_1
#integer to control number of replicates per model
reps = int(30)
#max number of iterations per optimization step (though see dadi user group for explanation)
maxiter = int(20)


print '\n', '\n',pop_ids, "Data for spectrum:"
print "projection", proj_1
print "grid pts", pts
print "sample sizes", fs_1.sample_sizes
print "Segregating sites",fs_1.S(), '\n', '\n'
#======================================================================================
# Now call the function with the relevant arguments.

# There are many models to test here. A brief definition is given for each, but the actual
# models are defined in the Models_2D.py script. The first 15 were implemented in Portik 
# et al. 2016 (doi: 10.1111/mec.14266), the following 9 are newer for various projects.

# Here it is set up to call each model one by one sequentially, which could finish relatively quickly.
# If it takes too long, create multiple verisions of this script, block out some models (use hashes or delete),
# and execute one version for every core you have available. It will greatly speed up these steps,
# and sometimes if extrapolations fail the script will crash too and this could prevent it from
# happening too many times.


# 1 Standard neutral model, populations never diverge
#Optimize_Functions.Optimize_Round1(pts, fs, outfile, reps, maxiter, "no_divergence")

# 2 Split into two populations, no migration.
#Optimize_Functions.Optimize_Round1(pts, fs, outfile, reps, maxiter, "no_mig")

# 3 Split into two populations, with continuous symmetric migration.
Optimize_Functions.Optimize_Round1(pts, fs, outfile, reps, maxiter, "sym_mig")

# 4 Split into two populations, with continuous asymmetric migration.
#Optimize_Functions.Optimize_Round1(pts, fs, outfile, reps, maxiter, "asym_mig")

# 5 Split with continuous symmetric migration, followed by isolation.
#Optimize_Functions.Optimize_Round1(pts, fs, outfile, reps, maxiter, "anc_sym_mig")

# 6 Split with continuous asymmetric migration, followed by isolation.
#Optimize_Functions.Optimize_Round1(pts, fs, outfile, reps, maxiter, "anc_asym_mig")

# 7 Split with no gene flow, followed by period of continuous symmetrical gene flow.
#Optimize_Functions.Optimize_Round1(pts, fs, outfile, reps, maxiter, "sec_contact_sym_mig")

# 8 Split with no gene flow, followed by period of continuous asymmetrical gene flow.
#Optimize_Functions.Optimize_Round1(pts, fs, outfile, reps, maxiter, "sec_contact_asym_mig")

# Newer Models.
# Split into two populations, with continuous symmetric migration, rate varying across two epochs.
#Optimize_Functions.Optimize_Round1(pts, fs, outfile, reps, maxiter, "sym_mig_twoepoch")

# Split into two populations, with continuous asymmetric migration, rate varying across two epochs.
#Optimize_Functions.Optimize_Round1(pts, fs, outfile, reps, maxiter, "asym_mig_twoepoch")



#===========================================================================
#clock the amount of time to complete the script
t_finish = datetime.now()
elapsed = t_finish - t_begin
print '\n', '\n', "-----------------------------------------------------------------------------------------------------"
print "Finished all analyses!"
print "Total time: {0} (H:M:S)".format(elapsed)
print "-----------------------------------------------------------------------------------------------------", '\n', '\n'
#===========================================================================


