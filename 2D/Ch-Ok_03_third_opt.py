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

daniel.portik@uta.edu -> danielportik@email.arizona.edu
October 2017
'''
#keep track of start time
t_begin = datetime.now()


#===========================================================================
#get snps file and convert into allele frequency spectrum object in dadi

#**************
#snps1 = "/FULL PATH TO /dadi_2pops_North_South_snps.txt"

#Create python dictionary from snps file
#dd1 = dadi.Misc.make_data_dict(snps1)

#**************
#pop_ids is a list which should match the populations headers of your SNPs file columns
pop_ids=["Ch", "Ok"]
#projection sizes, in ALLELES not individuals
proj_1 = [26,26]

#Convert this dictionary into folded AFS object
#[polarized = False] creates folded spectrum object
#fs_1 = dadi.Spectrum.from_data_dict(dd1, pop_ids=pop_ids, projections = proj_1, polarized = False)
fs_1 = dadi.Spectrum.from_file('C-O.sfs')


#print relevant info to screen
print '\n', '\n', "Data for spectrum:"
print "projection", proj_1
print "sample sizes", fs_1.sample_sizes
print "Segregating sites",fs_1.S(), '\n', '\n'

#======================================================================================

#"no_mig"
#3 Values
no_mig_params = [0.1721,0.6979,0.0727]

#"sym_mig"
#4 Values
sym_mig_params = [0.514,1.707,0.8254,2.1003]

#"asym_mig"
#5 Values
asym_mig_params = [0.3736,1.7884,1.5563,0.441,1.5978]

#"anc_sym_mig"
#5 Values
anc_sym_mig_params = [0.4122,1.181,1.2668,2.2156,0.014]

#"anc_asym_mig"
#6 Values
anc_asym_mig_params = [0.3201,1.143,1.4887,1.1573,1.2949,0.0045]

#"sec_contact_sym_mig"
#5 Values
sec_contact_sym_mig_params = [0.4018,1.2587,1.1656,0.2423,1.8125]

#"sec_contact_asym_mig"
#6 Values
sec_contact_asym_mig_params = [0.4627,2.7375,1.3249,0.2427,0.4278,2.7157]

#sym_mig_two_epoch
#6	values
sym_mig_twoepoch_params =	[0.4958,1.4542,0.1684,0.8967,0.4218,1.048 ]

#asym_mig_two_epoch
#8	values
asym_mig_twoepoch_params =[0.4341,1.6836,0.3663,1.1174,1.284,0.212,3.1821,0.2906]


#======================================================================================
#Input some of the basic reusable arguments here specific to your data set

#These are specific to your data set:
#**************
#grid choice
pts = [30,40,50]
#prefix for output file naming
outfile = "Ch_Ok"

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
#Optimize_Functions.Optimize_Round3(pts, fs, outfile, reps, maxiter, "no_divergence", no_divergence_params)

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
Optimize_Functions.Optimize_Round3(pts,	fs, outfile, reps, maxiter, "sym_mig_twoepoch",sym_mig_twoepoch_params)
							
#	Two	epoch	asym	mig			
Optimize_Functions.Optimize_Round3(pts,	fs, outfile,reps,maxiter, "asym_mig_twoepoch", asym_mig_twoepoch_params)	

#===========================================================================
#clock the amount of time to complete the script
t_finish = datetime.now()
elapsed = t_finish - t_begin
print '\n', '\n', "-----------------------------------------------------------------------------------------------------"
print "Finished all analyses!"
print "Total time: {0} (H:M:S)".format(elapsed)
print "-----------------------------------------------------------------------------------------------------", '\n', '\n'
#===========================================================================



