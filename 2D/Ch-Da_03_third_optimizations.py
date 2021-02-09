import sys
import os
import numpy as np
import dadi
import pylab
import Models_2D
import Optimize_Functions
from datetime import datetime

'''
usage: python dadi_2D_03_third_optimizations.py

Requires the Models_2D.py and Optimize_Functions.py scripts to be in same working directory. 
This is where all the population models and functions are stored for this script. 

Dan Portik
daniel.portik@uta.edu -> danielportik@email.arizona.edu
October 2017
'''
#keep track of start time
t_begin = datetime.now()


#**************
# PARAMETERS TO ENTER PAIRWISE
#pop_ids is a list which should match the populations headers of your SNPs file columns
pop_ids=["Ch", "Da"]
#projection sizes, in ALLELES not individuals
proj_1 = [26,16]

#grid choice
pts = [30,40,50]
#prefix for output file naming
outfile = "C-D"
3
fs_1 = dadi.Spectrum.from_file('C-D.sfs')

#"no_mig"
#3 Values
no_mig_params = [0.3183,0.6859,0.1296]

#"sym_mig"
#4 Values
sym_mig_params =[0.6443,1.294,0.6798,1.4109]

#"asym_mig"
#5 Values
asym_mig_params = [1.476,3.9291,0.3783,0.1625,5.9615]

#"anc_sym_mig"
#7 Values
anc_sym_mig_params =[1.6899,3.3513,0.3238,6.4628,0.0858]

#"anc_asym_mig"
#8 Values
anc_asym_mig_params = [0.2206,0.342,1.8993,3.9617,4.9511,0.0118]

#"sec_contact_sym_mig"
#5 Values
sec_contact_sym_mig_params = [0.4517,0.9016,0.8876,0.0038,0.9984]

#"sec_contact_asym_mig"
#6 Values
sec_contact_asym_mig_params = [0.8097,1.8909,0.6445,0.3532,0.0087,2.0606]

#sym_mig_two_epoch
#6	values
sym_mig_twoepoch_params = [1.2964,2.51,0.4982,0.334,0.0717,4.4047]

#asym_mig_two_epoch
#8	values
asym_mig_twoepoch_params = [0.3133,0.4143,17.3388,2.8563,1.0669,2.1757,1.369,0.3208]

#**************

#print relevant info to screen
print '\n', '\n',pop_ids,  "Data for spectrum:"
print "projection", proj_1
print "sample sizes", fs_1.sample_sizes
print "Segregating sites",fs_1.S(), '\n', '\n'

#======================================================================================

#These can be left alone, unless you want more searches:
#integer to control number of replicates per model
reps = int(50)
#max number of iterations per optimization step (though see dadi user group for explanation)
maxiter = int(50)

#======================================================================================
# Now call the function with the relevant arguments.

fs= fs_1

# Split into two populations, no migration.
Optimize_Functions.Optimize_Round3(pts, fs, outfile, reps, maxiter, "no_mig", no_mig_params)

# Split into two populations, with continuous symmetric migration.
Optimize_Functions.Optimize_Round3(pts, fs, outfile, reps, maxiter, "sym_mig", sym_mig_params)

# Split into two populations, with continuous asymmetric migration.
Optimize_Functions.Optimize_Round3(pts, fs, outfile, reps, maxiter, "asym_mig", asym_mig_params)

# Split with continuous symmetric migration, followed by isolation.
Optimize_Functions.Optimize_Round3(pts, fs, outfile, reps, maxiter, "anc_sym_mig", anc_sym_mig_params)

# Split with continuous asymmetric migration, followed by isolation.
Optimize_Functions.Optimize_Round3(pts, fs, outfile, reps, maxiter, "anc_asym_mig", anc_asym_mig_params)

# Split with no gene flow, followed by period of continuous symmetrical gene flow.
Optimize_Functions.Optimize_Round3(pts, fs, outfile, reps, maxiter, "sec_contact_sym_mig", sec_contact_sym_mig_params)

# Split with no gene flow, followed by period of continuous asymmetrical gene flow.
Optimize_Functions.Optimize_Round3(pts, fs, outfile, reps, maxiter, "sec_contact_asym_mig", sec_contact_asym_mig_params)

#	Two	epoch	model	sym	mig	
Optimize_Functions.Optimize_Round3(pts,	fs,outfile, reps, maxiter,"sym_mig_twoepoch",sym_mig_twoepoch_params)
						
#	Two	epoch	asym	mig		
Optimize_Functions.Optimize_Round3(pts,	fs, outfile, reps, maxiter, "asym_mig_twoepoch", asym_mig_twoepoch_params)

#===========================================================================
#clock the amount of time to complete the script
t_finish = datetime.now()
elapsed = t_finish - t_begin
print '\n', '\n', "-----------------------------------------------------------------------------------------------------"
print "Finished all analyses!"
print "Total time: {0} (H:M:S)".format(elapsed)
print "-----------------------------------------------------------------------------------------------------", '\n', '\n'
#===========================================================================


