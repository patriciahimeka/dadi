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

#"anc_sym_mig"
#5 Values
anc_sym_mig_params = [0.7858,	1.2965,	3.8012,	2.6331,	0.0154]

#"anc_asym_mig"
#6 Values
anc_asym_mig_params = [2.0371,	6.303,	2.3277,	0.6339,	9.9115,	0.0768]

#"sec_contact_sym_mig"
#5 Values
sec_contact_sym_mig_params = [1.4501,	2.7931,	1.5826,	0.3157,	5.3692]

#"sec_contact_asym_mig"
#6 Values
sec_contact_asym_mig_params = [0.3108,	2.5649,	10.76,	0.7624,	5.7009,	0.5631]

#sym_mig_two_epoch
#6	values
sym_mig_twoepoch_params =	[2.0384,	3.8796,	0.6922,	1.0689,	0.2146,	7.4406]

#asym_mig_two_epoch
#8	values
asym_mig_twoepoch_params =[0.4944,	1.6558,	13.9356,	3.5927,	5.1233,	2.0408,	0.9365,	0.9516]



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


# Split with continuous symmetric migration, followed by isolation.
Optimize_Functions.Optimize_Round3(pts, fs, outfile, reps, maxiter, "anc_sym_mig", anc_sym_mig_params)

# Split with continuous asymmetric migration, followed by isolation.
Optimize_Functions.Optimize_Round3(pts, fs, outfile, reps, maxiter, "anc_asym_mig", anc_asym_mig_params)


# Split with no gene flow, followed by period of continuous symmetrical gene flow.
Optimize_Functions.Optimize_Round3(pts, fs, outfile, reps, maxiter, "sec_contact_sym_mig", sec_contact_sym_mig_params)

# Split with no gene flow, followed by period of continuous asymmetrical gene flow.
Optimize_Functions.Optimize_Round3(pts, fs, outfile, reps, maxiter, "sec_contact_asym_mig", sec_contact_asym_mig_params)


# Split into two populations, with continuous symmetric migration, rate varying across two epochs.
Optimize_Functions.Optimize_Round3(pts, fs, outfile, reps, maxiter, "sym_mig_twoepoch", sym_mig_twoepoch_params)

# Split into two populations, with continuous asymmetric migration, rate varying across two epochs.
Optimize_Functions.Optimize_Round3(pts, fs, outfile, reps, maxiter, "asym_mig_twoepoch", asym_mig_twoepoch_params)

#===========================================================================
#clock the amount of time to complete the script
t_finish = datetime.now()
elapsed = t_finish - t_begin
print '\n', '\n', "-----------------------------------------------------------------------------------------------------"
print "Finished all analyses!"
print "Total time: {0} (H:M:S)".format(elapsed)
print "-----------------------------------------------------------------------------------------------------", '\n', '\n'
#===========================================================================



