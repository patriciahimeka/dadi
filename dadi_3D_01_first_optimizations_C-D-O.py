import sys
import os
import numpy as np
import dadi
import matplotlib
from datetime import datetime
import Models_3D
'''
usage: python dadi_3D_01_first_optimizations_C-D-O.py

Requires the Models_3D.py script to be in same working directory.
This is where all the population model functions are stored. This is written for the
models specifically found in that script. 

Script will perform optimizations from multiple starting points using a
3-fold perturbed set of random starting values for parameters. The output for
each model is a tab-delimited text file which can be opened and sorted to
find the best scoring replicate. The parameter values from this run should
be used as the starting values in the next script, "dadi_3D_02_second_optimizations.py".
The order of output parameters matches the input order, so they can be copied
and pasted into the next script within the appropriate list. 

Requires user to edit sections of code marked with #**************

You'll absolutely need to provide the path to your SNPs input file
along with your specific projections and population labels. 


############################################
Written for Python 2.7
Python modules required:
-Numpy
-Scipy
-Matplotlib
-dadi
############################################

Dan Portik
daniel.portik@uta.edu
April 2017
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

def Three_Pop_Models(pts, fs, outfile, reps, y, model_name):
    print '\n',"============================================================================"
    print "Beginning analysis of {}".format(model_name)
    print "============================================================================"

    #create output file
    outname = "Round1_{0}_{1}_optimized.txt".format(outfile,model_name)
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
        params=[1,1,1,1,1,1]
        print "parameter set = [nu1, nuA, nu2, nu3, T1, T2]"

        for i in range(1,x):
            fh_out = open(outname, 'a')
            fh_out.write("Split with No Migration"+'\t')
            fh_out.write("parameter set = [nu1, nuA, nu2, nu3, T1, T2]"+'\t')
            fh_out.write("{}\t".format(i))
            print '\n', "Replicate {}:".format(i)

            #perturb initial guesses
            params_perturbed = dadi.Misc.perturb_params(params, fold=3, upper_bound=upper_bound, lower_bound=lower_bound)

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



    elif model_name == "starsplit":

        #####################################
        #Split between 3 pops at same time
        #####################################
        print "---------------------------------------------------"
        print "Starsplit between 3 populations",'\n','\n'

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
            fh_out.write("parameter set = [nu1, nu2, nu3, m12, m21, m13 ,m31, m3, T]"+'\t')
            fh_out.write("{}\t".format(i))
            print '\n', "Replicate {}:".format(i)

            #perturb initial guesses
            params_perturbed = dadi.Misc.perturb_params(params, fold=3, upper_bound=upper_bound, lower_bound=lower_bound)

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





    elif model_name == "split_asymmig_all":
        
        #####################################
        #Split with asymmetric migration
        #####################################
        print "---------------------------------------------------"
        print "Split with Asymmetric Migration",'\n','\n'
        
        #first call a predefined model
        model_call = Models_3D.split_asymmig_all

        #create an extrapolating function 
        func_exec = dadi.Numerics.make_extrap_log_func(model_call)

        #create parameter list for optimization, set bounds for search
        lower_bound = [0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0, 0]
        upper_bound = [30, 30, 30, 30, 20, 20, 20, 20, 20, 20, 20, 10, 10]
        params=[1,1,1,1,1,1,1,1,1,1,1,1,1]
        print "parameter set = [nu1, nuA, nu2, nu3, mA1, mA2, m12, m21, m13, m31, m3, T1, T2]"
        
        for i in range(1,x):
            fh_out = open(outname, 'a')
            fh_out.write("Split with Asymmetric Migration"+'\t')
            fh_out.write("parameter set = [nu1, nuA, nu2, nu3, mA1, mA2, m12, m21, m13, m31, m3, T1, T2]"+'\t') 
            fh_out.write("{}\t".format(i))
            print '\n', "Replicate {}:".format(i)

            #perturb initial guesses
            params_perturbed = dadi.Misc.perturb_params(params, fold=3, upper_bound=upper_bound, lower_bound=lower_bound)

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
        params=[1,1,1,1,1,1,1,1,1,1]
        print "parameter set = [nu1, nuA, nu2, nu3, mA, m1, m2, m3, T1, T2]"
        
        for i in range(1,x):
            fh_out = open(outname, 'a')
            fh_out.write("Split with Symmetric Migration"+'\t')
            fh_out.write("parameter set = [nu1, nuA, nu2, nu3, mA, m1, m2, m3, T1, T2]"+'\t') 
            fh_out.write("{}\t".format(i))
            print '\n', "Replicate {}:".format(i)

            #perturb initial guesses
            params_perturbed = dadi.Misc.perturb_params(params, fold=3, upper_bound=upper_bound, lower_bound=lower_bound)

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
        params=[1,1,1,1,1,1,1,1,1]
        print "parameter set = [nu1, nuA, nu2, nu3, mA, m1, m2, T1, T2]"
        
        for i in range(1,x):
            fh_out = open(outname, 'a')
            fh_out.write("Split with Adjacent Symmetric Migration"+'\t')
            fh_out.write("parameter set = [nu1, nuA, nu2, nu3, mA, m1, m2, T1, T2]"+'\t')
            fh_out.write("{}\t".format(i))
            print '\n', "Replicate {}:".format(i)

            #perturb initial guesses
            params_perturbed = dadi.Misc.perturb_params(params, fold=3, upper_bound=upper_bound, lower_bound=lower_bound)

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

    
        
#====================================================================================

#**************
#Input some of the basic reusable arguments here
pts = [30,40,50]
fs = fs_1
outfile = "C-D-O"
reps = int(20)
maxiter = int(15)

#**************

#Models from the Models_3D.py script
#Three_Pop_Models(pts, fs, outfile, reps, maxiter, "split_nomig")
#Three_Pop_Models(pts, fs, outfile, reps, maxiter, "split_symmig_all")
Three_Pop_Models(pts, fs, outfile, reps, maxiter, "split_asymmig_all")
#Three_Pop_Models(pts, fs, outfile, reps, maxiter, "starsplit")

#===========================================================================
#clock it!

t_finish = datetime.now()
elapsed = t_finish - t_begin

print '\n', '\n', "-----------------------------------------------------------------------------------------------------"
print "Finished all analyses!"
print "Total time: {0} (H:M:S)".format(elapsed)
print "-----------------------------------------------------------------------------------------------------", '\n', '\n'
#===========================================================================


