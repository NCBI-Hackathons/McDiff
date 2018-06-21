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




def main(args):
    return 0

if __name__ == '__main__':
    import sys
    sys.exit(main(sys.argv))
