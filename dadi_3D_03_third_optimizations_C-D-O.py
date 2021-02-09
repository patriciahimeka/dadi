import sys
import os
import numpy as np
import dadi
import matplotlib
from datetime import datetime
import Models_3D
'''
usage: python dadi_3D_03_third_optimizations.py

Requires the Models_3D.py script to be in same working directory.
This is where all the population model functions are stored. This is written for the
models specifically found in that script. 

Script will perform optimizations from multiple starting points using a
1-fold perturbed set of random starting values for parameters. The output for
each model is a tab-delimited text file which can be opened and sorted to
find the best scoring replicate. 

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

#load fs
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
    print "Beginning 3rd optimizations analysis of {}".format(model_name)
    print "============================================================================"

    #create output file
    outname = "Round3_{0}_{1}_optimized.txt".format(outfile,model_name)
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
            params_perturbed = dadi.Misc.perturb_params(params, fold=1, upper_bound=upper_bound, lower_bound=lower_bound)

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
            params_perturbed = dadi.Misc.perturb_params(params, fold=1, upper_bound=upper_bound, lower_bound=lower_bound)

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
            params_perturbed = dadi.Misc.perturb_params(params, fold=1, upper_bound=upper_bound, lower_bound=lower_bound)

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
            params_perturbed = dadi.Misc.perturb_params(params, fold=1, upper_bound=upper_bound, lower_bound=lower_bound)

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
            fh_out.write("parameter set = [nu1, nu2, nu3, m12, m21, m13, m31, m3, T]"+'\t')
            fh_out.write("{}\t".format(i))
            print '\n', "Replicate {}:".format(i)

            #perturb initial guesses
            params_perturbed = dadi.Misc.perturb_params(params, fold=1, upper_bound=upper_bound, lower_bound=lower_bound)

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





#======================================================================================
# Finally, execute model with appropriate arguments

# model_name = from this list ["split_nomig", "split_symmig_all", "split_symmig_adjacent", "refugia_1",
#        "refugia_2", "refugia_3", "ancmig_3", "ancmig_2", "ancmig_1"]


#===========================================================================
# enter best param values for each model here, presumably you will get these
# from the outputs of the previous script, "dadi_3D_01_first_optimizations.py"

#************** "split_nomig"
# 6 Values
split_nomig_params = [2.1127,0.9016,7.4879,0.8309,0.0737,0.5939]

#************** "split_symmig_all"
# 10 Values
split_symmig_all_params = [0.6047,1.8755,2.7982,0.1294,1.1115,0.4329,1.9536,0.5993,1.5957,0.3824]

#************** "split_symmig_adjacent"
# 9 Values
split_symmig_adjacent_params = [2.911,6.4277,4.4721,1.6473,3.1932,0.2043,0.2706,6.2816,1.7949]

#************** "split_asymmig_all"
# 13 Values
split_asymmig_all_params = [0.2849,	0.3374,	0.545,	0.8678,	1.9271,	2.0997,	1.2006,	0.4029,	1.0596,	0.7082,	4.6644,	1.5479,	0.8986]

#************** "starsplit"
# 9 Values
starsplit_params = [0.3041,	0.5704,	1.0414,	0.9484,	0.4721,	1.108,	0.6374,	3.7644,	0.5398]


#**************
#Input some of the basic reusable arguments here
pts = [30,40,50]
fs = fs_1
outfile = "C-D-O"
reps = int(50)
maxiter = int(100)

#**************
# Here it is set up to call each model one by one sequentially, but this will take a very long time.
# I recommend blocking out all models except one (use a hash or delete), make several
# copies of the script and execute one model version for every core you have available.
# It will greatly speed up these steps, and sometimes if extrapolations fail the
# script will crash too.
# There are 20 models to test here. 

#Models from the Models_3D.py script

#Three_Pop_Models(pts, fs, outfile, reps, maxiter, "starsplit", starsplit_params)
Three_Pop_Models(pts, fs, outfile, reps, maxiter, "split_asymmig_all", split_asymmig_all_params)

#Three_Pop_Models(pts, fs, outfile, reps, maxiter, "split_nomig", split_nomig_params)
#Three_Pop_Models(pts, fs, outfile, reps, maxiter, "split_symmig_all", split_symmig_all_params)
#Three_Pop_Models(pts, fs, outfile, reps, maxiter, "split_symmig_adjacent", split_symmig_adjacent_params)
#Three_Pop_Models(pts, fs, outfile, reps, maxiter, "refugia_1", refugia_1_params)
#Three_Pop_Models(pts, fs, outfile, reps, maxiter, "refugia_2", refugia_2_params)
#Three_Pop_Models(pts, fs, outfile, reps, maxiter, "refugia_3", refugia_3_params)
#Three_Pop_Models(pts, fs, outfile, reps, maxiter, "ancmig_3", ancmig_3_params)
#Three_Pop_Models(pts, fs, outfile, reps, maxiter, "ancmig_2", ancmig_2_params)
#Three_Pop_Models(pts, fs, outfile, reps, maxiter, "ancmig_1", ancmig_1_params)



#===========================================================================
#clock it!

t_finish = datetime.now()
elapsed = t_finish - t_begin

print '\n', '\n', "-----------------------------------------------------------------------------------------------------"
print "Finished all analyses!"
print "Total time: {0} (H:M:S)".format(elapsed)
print "-----------------------------------------------------------------------------------------------------", '\n', '\n'
#===========================================================================

