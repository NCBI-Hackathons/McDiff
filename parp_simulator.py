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

def wrapper(data_file, roi_file, mask_file, bound_d, exp_time, sigmaD, sigmaF, mcmc_temp, offset, mcmc_steps):
    #exp_time is the sim_len but dont want the user to have to do the calcualations
    #default sigmaD =2    parameters in the mcmc
    #default sigmaf =.05  parameters in the mcmc
    #recomend/default  mcmc_steps = 200
    #default bound_d = 20  the upper bound for the for you think d could possibly be
    #default offset = 10.5
    #default mcmc_temp = 1 when to calc likelyhood ratio, devide be estmate of noise ***NOT TEMPURATURE OF EXPERIMENT
    #default percent_bleached = .46 that is when grean(gfp) is used 

    exec(open("sims.py").read())

    #parse all data files then create the roi and the nucleus
    mask = parse_mask(mask_file)
    roi_cords = parse_roi(roi_file)
    roi = Polygon([(roi_cords[1], roi_cords[0]), (roi_cords[1], roi_cords[0]+roi_cords[2]), (roi_cords[1]+d[3], roi_cords[0]+roi_cords[2]), (roi_cords[1]+roi_cords[3], roi_cords[0])])
    nuc = Polygon(list(zip(mask[0,:], mask[1,:])))

    data_pre, data = parse_data(data_file, offset) # offset origianlly 10.5
    data_norm = data[1,:] / np.mean(data_pre[1,:])

    #sim_len = 650 recomended

    ##MCMC: mcmc_steps = 200 suggestion

    OP, E, AP, bool_flag_1, bool_flag_2, Iterate_ended = MCMC(4, .18, .5, nuc, roi, mcmc_steps, mcmc_temp, sigmaD, sigmaF, 0, 1, 0, bounds_d) #.18 is timestep
    #bool_flag_1 and bool_flag_2 and interate ends are for debugging, if either is true then the simulation has gone wrong


    lo_mejor = Error.argmin() #index of parameter optimal
    los_mejores = AP[:, lo_mejor] #best parameters
    epx_time = int(exp_time/.18) #translate time to steps

    stuck_in_roi, roi_pre = simulate(los_mejores[0], los_mejores[1], 0.5, nuc, roi, exp_time)

    stuck_norm = stuck_in_roi / roi_pre
    stuck_time = np.arange(sim_len+1) * 0.18 #converts array indices into seconds

    fig, ax = plt.subplots(1)
    ax[0].plot(stuck_time, stuck_norm, ".", label = "Simulation")
    ax[0].plot(data[0,:], data_norm, ".", label = "Data")
    ax[0].legend()
    ax[0].set_xlabel("Time (s)")
    ax[0].set_ylabel("Fraction of Proteins Bound/Baseline")
    ax[1].hist(AP[0,:],20)
    ax[2].hist(AP[1,:],20)
    ax[3].plot(Error)
    #figure to return
    plt.savefig(fit_plot_name)
    #CSV to return to the user, as in simulated data results
    newfile = open("MCMC_results.txt", 'w')
    for i in range(len(E)):
        text = "{0} {1} {2}\n".format(E[i], AP[0,i], AP[1,i])
        newfile.write(fit_data_name)


    N = 50
    f_bleached = .46 #user input about the color used in experiment
    fmin = 0
    fmax = 1
    dmin = 0
    D, F, ret_error = rand_sam(f_bleached, nuc, roi, N, fmin, fmax, dmin, bound_d)


	return fit_plot_name, fit_data_name, D_final, F_final, ret_error



def main(args):

    wrapper(data_file, roi_file, mask_file, bounds_f, bounds_d, exp_time, sigmaD, sigmaF temp, offset, mcmc_steps)
    return 0

if __name__ == '__main__':
    import sys
    sys.exit(main(sys.argv))
