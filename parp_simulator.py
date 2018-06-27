#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  parp_simulator.py
#
#  Copyright 2018
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with this program; if not, write to the Free Software
#  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
#  MA 02110-1301, USA.
#
#
"""
From parpsimulator.m:

    %PARPSIMULATOR  Class to simulate free diffusion of PARP molecules
    %
    % parpsimulator Properties:
    %   pxSize - Physical units of the simulated grid in microns/pixel
    %   deltaT - Physical units of each time step in seconds/time step
    %   diffusionCoeff - Diffusion coefficient in microns^2/seconds
    %   mparpFrac - Percent fraction of mobile PARP particles
    %   unbleachedFrac - Percent fraction of photobleached particles
    %   numParticles - Number of particles to simulate
    %   numSteps - Number of time steps to simulate
    %   outputMovie - If true, a movie will be saved
    %
    % parpsimulator Methods:
    %   simulate - Run the simulation using specified parameters
    %   parameterSweep - Search for best fit diffusion coefficient and mobile PARP fraction
    %
    % Example:
    %   %Create an instance of the class
    %   sim = parpsimulator
    %
    %   %Run the simulation with default settings
    %
"""
import os.path
import numpy as np
import sims
import optimization  #exec(open('sims.py').read())
import matplotlib.pyplot as plt


def wrapper(data_file, roi_file, mask_file, bound_d, exp_time, percent_bleached, sigmaD, sigmaF, mcmc_temp, offset, dam_ind, mcmc_steps, file_name):
    #exp_time is the sim_len but dont want the user to have to do the calcualations
    #default sigmaD =2    parameters in the mcmc - int
    #default sigmaf =.05  parameters in the mcmc - float
    #recomend/default  mcmc_steps = 200 - int
    #default bound_d = 20  - int the upper bound for the for you think d could possibly be
    #default offset = 10.5 - float
    #default dam_ind = 6.0, int
    #default mcmc_temp = 1 - int  when to calc likelyhood ratio, devide be estmate of noise ***NOT TEMPURATURE OF EXPERIMENT
    #default percent_bleached = .54 that is when grean(gfp) is used, that is ammount that becomes bleached - float
    #mcmc_steps = 200 suggestion -int

    #parse all data files then create the roi and the nucleus:
    mask = sims.parse_mask(mask_file)
    roi_cords = sims.parse_roi(roi_file)
    roi = sims.Polygon([(roi_cords[1], roi_cords[0]), (roi_cords[1], roi_cords[0]+roi_cords[2]), (roi_cords[1]+roi_cords[3], roi_cords[0]+roi_cords[2]), (roi_cords[1]+roi_cords[3], roi_cords[0])])
    nuc = sims.Polygon(list(zip(mask[0,:], mask[1,:])))

    data_pre, data = sims.parse_data(data_file, offset, dam_ind) # offset origianlly 10.5
    data_norm = data[1,:] / np.mean(data_pre[1,:])

    sim_len = mcmc_steps #recomended
    x0, y0 = sims.init_sim(12000, nuc)
    s1 = .18
    s2 = .18
    N = 100
    L = 3

    #OP, Error, AP, bool_flag_1, bool_flag_2, Iterate_ended = optimization.MCMC(4, .18, percent_bleached, nuc, roi, mcmc_steps, mcmc_temp, sigmaD, sigmaF, 0, 1, 0, bound_d, sim_len, data, data_pre, data_norm, x0, y0) #.18 is timestep
    #bool_flag_1 and bool_flag_2 and interate ends are for debugging, if either is true then the simulation has gone wrong

    results = optimization.CF(percent_bleached, nuc, roi, 0, 1, 0, bound_d, s1, s2, N, L, x0, y0, sim_len, data, data_norm)
    print("done with cf")
    x = results[2,:].argmin()
    results[2,x]

    exp_time = int(exp_time/.18) #translate time to steps

    stuck_norm = sims.simulate(results[0,x], results[1,x], 0.5, nuc, roi, sim_len, x0, y0)
    stuck_time = np.arange(sim_len+1) * 0.18 #converts array indices into seconds
    ret_error = sims.compute_error(data, data_norm, stuck_time, stuck_norm)

    fig, ax = plt.subplots(2)
    ax[0].scatter(results[0,:], results[1,:], 10., results[2,:], ".")
    ax[1].plot(stuck_time, stuck_norm, ".", label = "Simulation")
    ax[1].plot(data[0,:], data_norm, ".", label = "Data")
    ax[1].legend()
    ax[1].set_xlabel("Time (s)")
    ax[1].set_ylabel("Fraction of Proteins Bound/Baseline")
    plt.savefig(file_name+"fig") #saves figure to directory
    savepath = 'C:/GitHub/McDiff/app'
    os.path.join(savepath+'/graph', file_name+"_Figure")
    print("figure saved")

    #CSV to return to the user, as in simulated data results saved in results
    save_path = 'C:/GitHub/app/static/'
    newfile = open(file_name+"_MCMC_results.csv", 'w')

    for i in range(len(ret_error)):
        text = "{0} {1} {2}\n".format(ret_error[i], results[0,i], results[1,i])
        newfile.write(text)
    #trying here to attach the data_file name onto return results name
    os.path.join(save_path+'/results', file_name+"_MCMC_results.csv")
    #print(data_file+"_MCMC_results.csv")


    resid_plot_name = file_name+'_resid_plot'
    fit_data_name = file_name+"_MCMC_results.csv"

    #results[0,0:N] = D
    #results[1,0:N] = F
    #results[2,0:N] = E
    #these are all arrays
    D_final = results[0,0:N]
    F_final = results[1,0:N]
    ret_error = results[2,0:N]

    return resid_plot_name, fit_data_name, ret_error, D_final, F_final



def main(args):
    mask_file = "./test_files/1.31.18_GFPP1_Hela_1min_002NuclMask.txt"
    roi_file = "./test_files/1.31.18_GFPP1_Hela_1min_002ROI.txt"
    data_file = "./test_files/1.31.18_GFPP1_Hela_1min_002.csv"
    bound_d = 12
    exp_time = 30
    sigmaD = 2
    sigmaF = .05
    mcmc_temp = 1
    offset = 10.5
    mcmc_steps = 250
    percent_bleached = .56
    file_name = "FILE"

    #data_file, roi_file, mask_file, bound_d, exp_time, sigmaD, sigmaF, mcmc_temp, offset, mcmc_steps = input('Please input in this order: \n data_file, roi_file, mask_file, bound_d, exp_time, sigmaD, sigmaF, mcmc_temp, offset, mcmc_steps' )
    print(wrapper(data_file, roi_file, mask_file, bound_d, exp_time, percent_bleached, sigmaD, sigmaF, mcmc_temp, offset, mcmc_steps, file_name))
    return 0

if __name__ == '__main__':
    import sys
    sys.exit(main(sys.argv))
