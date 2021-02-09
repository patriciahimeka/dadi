import sys
import os
import numpy as np
import dadi
import matplotlib
from datetime import datetime
import Models_3D
'''
Optimization functions used in Wepfer et al., "The oceanographic isolation of the Ogasawara Islands and genetic divergence in a reef-building coral".

Script taken from Dan Portik's pipeline: https://github.com/dportik/dadi_pipeline/blob/master/README.md#V  (version 2 for python 2.7, 2018)
CITE: https://github.com/dportik/dadi_pipeline/blob/master/README.md#V

usage: python dadi_3D_02_second_optimizations_C-D-O.py
Requires the Models_3D.py script to be in same working directory.

############################################
Written for Python 2.7
Python modules required:
-Numpy
-Scipy
-Matplotlib
-dadi
############################################

Patricia Wepfer
August 2018
'''
t_begin = datetime.now()

#===========================================================================
#**************
#projection sizes, in ALLELES not individuals
proj_1 = [26,16,26]
#pop_ids is a list which should match the populations headers of your SNPs file columns
pop_ids=['C','D','O']

#Convert this dictionary into folded AFS object
#[polarized = False] creates folded spectrum object
#fs_1 = dadi.Spectrum.from_data_dict(dd1, pop_ids=pop_ids, projections = proj_1, polarized = False)
fs_1=dadi.Spectrum.from_file('C-D-O.sfs')

print '\n', '\n', "Data for spectrum:"
print "projection", proj_1
print "sample sizes", fs_1.sample_sizes
print "Segregating sites",fs_1.S(), '\n', '\n'


#======================================================================================
#create function to take in grid size and freq spectrum for a population, then run all 
#models 'x' times from optimized parameters from heavily perturbed basic starting values
#write output of parameter optimizations, theta, likelihood, and AIC for each replicate for 
#each model to file

def Three_Pop_Models(pts, fs, outfile, reps, y, model_name, params):
    print '\n',"============================================================================"
    print "Beginning analysis of {}".format(model_name)
    print "============================================================================"

    #create output file
    outname = "Round2_{0}_{1}_optimized.txt".format(outfile,model_name)
    fh_out = open(outname, 'a')
    fh_out.write("Model"+'\t'+"param_set"+'\t'+"Replicate"+'\t'+"log-likelihood"+'\t'+"theta"+'\t'+"AIC"+'\t'+"optimized_params"+'\n')
    fh_out.close()
    
    #variable to control number of loops per model (1 to x)
    x = int(reps) + int(1)
    
    ######################################################################################################################
    #Basic models
    ######################################################################################################################
    
    if model_name == "split_nomig":
        
        #####################################
        #Split with no migration
        #####################################
        print "---------------------------------------------------"
        print "Split with No Migration",'\n','\n'
        
        #first call a predefined model
        model_call = Models_3D.split_nomig

        #create an extrapolating function 
        func_exec = dadi.Numerics.make_extrap_log_func(model_call)

        #create parameter list for optimization, set bounds for search
        lower_bound = [0.01, 0.01, 0.01, 0.01, 0, 0]
        upper_bound = [30, 30, 30, 30, 10, 10]
        print "parameter set = [nu1, nuA, nu2, nu3, T1, T2]"

        for i in range(1,x):
            fh_out = open(outname, 'a')
            fh_out.write("Split with No Migration"+'\t')
            fh_out.write("parameter set = [nu1, nuA, nu2, nu3, T1, T2]"+'\t')
            fh_out.write("{}\t".format(i))
            print '\n', "Replicate {}:".format(i)
            print "starting parameters = ", params

            #perturb initial guesses
            params_perturbed = dadi.Misc.perturb_params(params, fold=2, upper_bound=upper_bound, lower_bound=lower_bound)

            #run optimization 
            params_opt = dadi.Inference.optimize_log_fmin(params_perturbed, fs, func_exec, pts,lower_bound=lower_bound, upper_bound=upper_bound, verbose=1, maxiter=y)
            print '\n',"optimized parameters = ", params_opt
            
            #simulate the model with the optimized parameters
            sim_model = func_exec(params_opt, fs.sample_sizes, pts)

            #calculate likelihood
            ll = dadi.Inference.ll_multinom(sim_model, fs)
            ll = np.around(ll, 2)
            print "likelihood = ", ll
            fh_out.write("{}\t".format(ll))

            #calculate theta
            theta = dadi.Inference.optimal_sfs_scaling(sim_model, fs)
            theta = np.around(theta, 2)
            print "Theta = ", theta
            fh_out.write("{}\t".format(theta))

            #calculate AIC 
            aic = ( -2*( float(ll))) + (2*6)
            print "AIC = ", aic, '\n', '\n'
            fh_out.write("{}\t".format(aic))

            for p in params_opt:
                p = np.around(p, 4)
                fh_out.write("{}\t".format(p))
            fh_out.write('\n')
            fh_out.close()

        print "---------------------------------------------------", '\n'


    elif model_name == "split_asymmig_all":

        #####################################
        #Split with asymmetric migration
        #####################################
        print "---------------------------------------------------"
        print "Split with asymmetric Migration",'\n','\n'

        #first call a predefined model
        model_call = Models_3D.split_asymmig_all

        #create an extrapolating function
        func_exec = dadi.Numerics.make_extrap_log_func(model_call)

        #create parameter list for optimization, set bounds for search
        lower_bound = [0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0, 0]
        upper_bound = [30, 30, 30, 30, 20, 20, 20, 20, 20, 20, 20, 10, 10]
        print "parameter set = [nu1, nuA, nu2, nu3, mA1, mA2, m12, m21, m13, m31, m3, T1, T2]"

        for i in range(1,x):
            fh_out = open(outname, 'a')
            fh_out.write("Split with Asymmetric Migration"+'\t')
            fh_out.write("parameter set = [nu1, nuA, nu2, nu3, mA1, mA2, m12, m21, m13, m31, m3, T1, T2]"+'\t')
            fh_out.write("{}\t".format(i))
            print '\n', "Replicate {}:".format(i)
            print "starting parameters = ", params

            #perturb initial guesses
            params_perturbed = dadi.Misc.perturb_params(params, fold=2, upper_bound=upper_bound, lower_bound=lower_bound)

            #run optimization
            params_opt = dadi.Inference.optimize_log_fmin(params_perturbed, fs, func_exec, pts,lower_bound=lower_bound, upper_bound=upper_bound, verbose=1, maxiter=y)
            print '\n',"optimized parameters = ", params_opt

            #simulate the model with the optimized parameters
            sim_model = func_exec(params_opt, fs.sample_sizes, pts)

            #calculate likelihood
            ll = dadi.Inference.ll_multinom(sim_model, fs)
            ll = np.around(ll, 2)
            print "likelihood = ", ll
            fh_out.write("{}\t".format(ll))

            #calculate theta
            theta = dadi.Inference.optimal_sfs_scaling(sim_model, fs)
            theta = np.around(theta, 2)
            print "Theta = ", theta
            fh_out.write("{}\t".format(theta))

            #calculate AIC
            aic = ( -2*( float(ll))) + (2*10)
            print "AIC = ", aic, '\n', '\n'
            fh_out.write("{}\t".format(aic))

            for p in params_opt:
                p = np.around(p, 4)
                fh_out.write("{}\t".format(p))
            fh_out.write('\n')
            fh_out.close()

        print "---------------------------------------------------", '\n'




    elif model_name == "starsplit":

        #####################################
        #Split between 3 pops at same time
        #####################################
        print "---------------------------------------------------"
        print "Starsplit between 3 populations, 2nd Optimization",'\n','\n'

        #first call a predefined model
        model_call = Models_3D.starsplit

        #create an extrapolating function
        func_exec = dadi.Numerics.make_extrap_log_func(model_call)

        #create parameter list for optimization, set bounds for search
        lower_bound = [0.01, 0.01, 0.01, 0.01,  0.01, 0.01, 0.01, 0.01, 0]
        upper_bound = [30, 30, 30, 20, 20, 20, 20, 20, 10]
        params=[1,1,1,1,1,1,1,1,1]
        print "parameter set = [nu1, nu2, nu3, m12, m21, m13, m31, m3, T]"

        for i in range(1,x):
            fh_out = open(outname, 'a')
            fh_out.write("Starsplit"+'\t')
            fh_out.write("parameter set = [nu1, nu2, nu3, m12, m21, m13, m31, m3, T]"+'\t')
            fh_out.write("{}\t".format(i))
            print '\n', "Replicate {}:".format(i)

            #perturb initial guesses
            params_perturbed = dadi.Misc.perturb_params(params, fold=2, upper_bound=upper_bound, lower_bound=lower_bound)

            #run optimization
            params_opt = dadi.Inference.optimize_log_fmin(params_perturbed, fs, func_exec, pts,lower_bound=lower_bound, upper_bound=upper_bound, verbose=1, maxiter=y)
            print '\n',"optimized parameters = ", params_opt

            #simulate the model with the optimized parameters
            sim_model = func_exec(params_opt, fs.sample_sizes, pts)

            #calculate likelihood
            ll = dadi.Inference.ll_multinom(sim_model, fs)
            ll = np.around(ll, 2)
            print "likelihood = ", ll
            fh_out.write("{}\t".format(ll))

            #calculate theta
            theta = dadi.Inference.optimal_sfs_scaling(sim_model, fs)
            theta = np.around(theta, 2)
            print "Theta = ", theta
            fh_out.write("{}\t".format(theta))

            #calculate AIC
            aic = ( -2*( float(ll))) + (2*10)
            print "AIC = ", aic, '\n', '\n'
            fh_out.write("{}\t".format(aic))

            for p in params_opt:
                p = np.around(p, 4)
                fh_out.write("{}\t".format(p))
            fh_out.write('\n')
            fh_out.close()

        print "---------------------------------------------------", '\n'


    elif model_name == "split_symmig_all":
        
        #####################################
        #Split with symmetric migration
        #####################################
        print "---------------------------------------------------"
        print "Split with Symmetric Migration",'\n','\n'
        
        #first call a predefined model
        model_call = Models_3D.split_symmig_all

        #create an extrapolating function 
        func_exec = dadi.Numerics.make_extrap_log_func(model_call)

        #create parameter list for optimization, set bounds for search
        lower_bound = [0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0, 0]
        upper_bound = [30, 30, 30, 30, 20, 20, 20, 20, 10, 10]
        print "parameter set = [nu1, nuA, nu2, nu3, mA, m1, m2, m3, T1, T2]"
        
        for i in range(1,x):
            fh_out = open(outname, 'a')
            fh_out.write("Split with Symmetric Migration"+'\t')
            fh_out.write("parameter set = [nu1, nuA, nu2, nu3, mA, m1, m2, m3, T1, T2]"+'\t') 
            fh_out.write("{}\t".format(i))
            print '\n', "Replicate {}:".format(i)
            print "starting parameters = ", params

            #perturb initial guesses
            params_perturbed = dadi.Misc.perturb_params(params, fold=2, upper_bound=upper_bound, lower_bound=lower_bound)

            #run optimization
            params_opt = dadi.Inference.optimize_log_fmin(params_perturbed, fs, func_exec, pts,lower_bound=lower_bound, upper_bound=upper_bound, verbose=1, maxiter=y)
            print '\n',"optimized parameters = ", params_opt

            #simulate the model with the optimized parameters
            sim_model = func_exec(params_opt, fs.sample_sizes, pts)

            #calculate likelihood
            ll = dadi.Inference.ll_multinom(sim_model, fs)
            ll = np.around(ll, 2)
            print "likelihood = ", ll
            fh_out.write("{}\t".format(ll))

            #calculate theta
            theta = dadi.Inference.optimal_sfs_scaling(sim_model, fs)
            theta = np.around(theta, 2)
            print "Theta = ", theta
            fh_out.write("{}\t".format(theta))

            #calculate AIC 
            aic = ( -2*( float(ll))) + (2*10)
            print "AIC = ", aic, '\n', '\n'
            fh_out.write("{}\t".format(aic))

            for p in params_opt:
                p = np.around(p, 4)
                fh_out.write("{}\t".format(p))
            fh_out.write('\n')
            fh_out.close()

        print "---------------------------------------------------", '\n'


    elif model_name == "split_symmig_adjacent":
        
        #####################################
        #Split with adjacent symmetric migration 
        #####################################
        print "---------------------------------------------------"
        print "Split with Adjacent Symmetric Migration",'\n','\n'
        
        #first call a predefined model
        model_call = Models_3D.split_symmig_adjacent

        #create an extrapolating function 
        func_exec = dadi.Numerics.make_extrap_log_func(model_call)

        #create parameter list for optimization, set bounds for search
        lower_bound = [0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0, 0]
        upper_bound = [30, 30, 30, 30, 20, 20, 20, 10, 10]
        print "parameter set = [nu1, nuA, nu2, nu3, mA, m1, m2, T1, T2]"
        
        for i in range(1,x):
            fh_out = open(outname, 'a')
            fh_out.write("Split with Adjacent Symmetric Migration"+'\t')
            fh_out.write("parameter set = [nu1, nuA, nu2, nu3, mA, m1, m2, T1, T2]"+'\t')
            fh_out.write("{}\t".format(i))
            print '\n', "Replicate {}:".format(i)
            print "starting parameters = ", params

            #perturb initial guesses
            params_perturbed = dadi.Misc.perturb_params(params, fold=2, upper_bound=upper_bound, lower_bound=lower_bound)

            #run optimization 
            params_opt = dadi.Inference.optimize_log_fmin(params_perturbed, fs, func_exec, pts,lower_bound=lower_bound, upper_bound=upper_bound, verbose=1, maxiter=y)
            print '\n',"optimized parameters = ", params_opt

            #simulate the model with the optimized parameters
            sim_model = func_exec(params_opt, fs.sample_sizes, pts)

            #calculate likelihood
            ll = dadi.Inference.ll_multinom(sim_model, fs)
            ll = np.around(ll, 2)
            print "likelihood = ", ll
            fh_out.write("{}\t".format(ll))

            #calculate theta
            theta = dadi.Inference.optimal_sfs_scaling(sim_model, fs)
            theta = np.around(theta, 2)
            print "Theta = ", theta
            fh_out.write("{}\t".format(theta))

            #calculate AIC 
            aic = ( -2*( float(ll))) + (2*9)
            print "AIC = ", aic, '\n', '\n'
            fh_out.write("{}\t".format(aic))

            for p in params_opt:
                p = np.around(p, 4)
                fh_out.write("{}\t".format(p))
            fh_out.write('\n')
            fh_out.close()

        print "---------------------------------------------------", '\n'

    

        
#======================================================================================
# Finally, execute model with appropriate arguments
# Three_Pop_Models(pts, fs, outfile, reps, maxiter, model_name, model_params):
# pts = grid choice (list of three numbers, ex. [20,30,40]
# fs = spectrum object name
# outfile = prefix for output naming (will result in "[prefix]_[model_name]_optimized_round1.txt")
# reps = integer to control number of replicates, ex. 10
# maxiter = max number of iterations per optimization step (not intuitive! see dadi user group)
# model_name = from this list ["split_nomig", "split_symmig_all", "split_symmig_adjacent", "refugia_1",
#        "refugia_2", "refugia_3", "ancmig_3", "ancmig_2", "ancmig_1"]
# model_params = list of parameter values to perturb and start optimizations from
#		ex. some_params = [4.787,0.465,7.071,1.879,0.181,0.635]



#===========================================================================
# enter best param values for each model here, presumably you will get these
# from the outputs of the previous script, "dadi_3D_01_first_optimizations.py"

#************** "split_nomig"
# 6 Values
#split_nomig_params = [4.7878,0.4657,7.0718,1.8793,0.1819,0.6351]

#************** "split_symmig_all"
# 10 Values
#split_symmig_all_params = [0.3645,2.3541,2.8694,0.1192,2.9680,0.2465,3.4061,2.2545,5.2956,0.3193]


#************** "split_asymmig_all"
# 13 Values
split_asymmig_all_params = [0.3376,	0.2146,	3.7781,	3.5552,	3.2884,	6.4483,	3.8234,	0.2693,	0.6021,	0.622,	0.2677,	0.6463,	4.3261]


#************** "starsplit"
# 9 Values
starsplit_params = [0.1785,	0.852,	0.3932,	2.1366,	3.1953,	2.42,	0.2563,	5.5322,	0.2892 ]


#************** "split_symmig_adjacent"
# 9 Values
#split_symmig_adjacent_params = [3.5506,3.6095,3.7157,1.1518,3.04188,0.5378,0.1777,5.8779,0.6790]




#**************
#Input some of the basic reusable arguments here
pts = [30,40,50]
fs = fs_1
outfile = "C-D-O"
reps = int(50)
maxiter = int(20)

#**************
# Here it is set up to call each model one by one sequentially, but this will take a very long time.
# I recommend blocking out all models except one (use a hash or delete), make several
# copies of the script and execute one model version for every core you have available.
# It will greatly speed up these steps, and sometimes if extrapolations fail the
# script will crash too.
# There are 20 models to test here. 

#Models from the Models_3D.py script 
#Three_Pop_Models(pts, fs, outfile, reps, maxiter, "split_nomig", split_nomig_params)
#Three_Pop_Models(pts, fs, outfile, reps, maxiter, "split_symmig_all", split_symmig_all_params)
Three_Pop_Models(pts, fs, outfile, reps, maxiter, "split_asymmig_all", split_asymmig_all_params)
#Three_Pop_Models(pts, fs, outfile, reps, maxiter, "starsplit", starsplit_params)
#Three_Pop_Models(pts, fs, outfile, reps, maxiter, "split_symmig_adjacent", split_symmig_adjacent_params)



#===========================================================================
#clock it!

t_finish = datetime.now()
elapsed = t_finish - t_begin

print '\n', '\n', "-----------------------------------------------------------------------------------------------------"
print "Finished all analyses!"
print "Total time: {0} (H:M:S)".format(elapsed)
print "-----------------------------------------------------------------------------------------------------", '\n', '\n'
#===========================================================================

