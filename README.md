![Packagist](https://img.shields.io/packagist/l/doctrine/orm.svg)

# McDiff
A Monte Carlo Approach for Estimating Diffusion Coefficients (add example use case)

# Abstract
Cell biologists can study the recruitment of DNA repair proteins to sites of DNA damage in live cells by using laser micro-irradiation to induce damage, which is known as Fluorescence Accumulation after DNA Damage (FADD). By monitoring the time-dependent accumulation of proteins at the sites of damage, or region of interest (ROI), the biologists then aim to calculate the coefficient of free diffusion (D) of each molecule of interest within the nucleus. This code simulates particles freely diffusing within a cell nucleus using a Monte Carlo model, becoming trapped at the ROI. Our code then seeks to optimize the best Diffusion constant (D) and fraction of protein that accumulates (F), and creates a heatmap for the best fit by using r-squared values comparing a fit to the data.


## What's the problem?
Before the implementation of this program, biologists would have to run either a mathematica or matlab scripts where each single simulation is run with manually entered coefficients. Quality of the fit to the experimental data was evaluated by visual inspection by the biologist. If the simulation was not aligned with the experimental data, then the biologist would change their coefficient values, then re-run the script and see if those values fit better, a tedious and labor-intensive process. Our new implementation runs simulations with multiple coefficients and calculates an r-squared error value that represents a best fit sample. Then the biologist can use that simulation, and be given the proper coefficients with the proper error.

# Keywords
 Free Diffusion Constant, Nucleus, Fluorescence Accumulation after DNA Damage (FADD)  

## Website (if applicable)


# What is McDiff?

![mcdiff monte carlo diffusion 1](https://user-images.githubusercontent.com/30201026/41800946-0bad625e-7635-11e8-9a17-afce8f03d06c.png)
![mcdiff monte carlo diffusion 2](https://user-images.githubusercontent.com/30201026/41800952-0fce1f22-7635-11e8-8ab3-c5c58d4ffb62.png)
![mcdiff monte carlo diffusion 3](https://user-images.githubusercontent.com/30201026/41800955-11ce0652-7635-11e8-89ba-8d115953c38b.png)
![mcdiff monte carlo diffusion 4](https://user-images.githubusercontent.com/30201026/41800960-13e82648-7635-11e8-8c7f-f38291615545.png)

# How to use McDiff
`$python3 sims.py`

# Software Workflow Diagram

![initial_flowchart_mcdiff](https://user-images.githubusercontent.com/23224399/41737510-beb45eb4-754c-11e8-816c-8720f1ae12e1.png)

## Methods

Programing Methods:  
  We started by getting comfortable with the experiments that were running in the nucleus. We then started going over the matlab and mathematica script that originally ran this program, with the lack of simulation.   
  We started in on programing by working on our parsing functions.  
  Then we had the outline of the simulated nucleus and populated the nucleus with the simulated particles. 
  We ran into a roadblock when the simulation of movement within the nucleus was taking a long time ~ 10 minutes. We found that the we could speed it up by almost x100 by using a vector library within the shapely library. Shapely is used to outline the nucleus and simulate the walls.   
  Once we had a simulation of the nucleus we started work on plotting the nucleus simulation results.   
  Then we worked on importing data from the experiments, and fitting simulated curves to the experimental curves. From there we could get an r-squared value and see which simulation gave us the best fitted curve with adjusted D and F values.   This is done by running each simulation with different D’s and F’s and finding the best simulation with the most fit curve.

# Installation options:

McDiff should be installed directly from Github.

### Installing McDiff from Github

### Step 1:
```
# Clone the repo
`$git clone https://github.com/NCBI-Hackathons/McDiff.git`
```

That's it! You should be ready to run.

### Configuration and dependencies
 
System requirements: python 3.6.5  
  Libraries: shapely, matplotlib, descartes, numpy, random


# Testing

We tested McDiff with the sample files found in the main directory. These files are:  
1. [Test Data](https://github.com/NCBI-Hackathons/McDiff/blob/master/test_files/1.31.18_GFPP1_Hela_1min_002.csv "Test Data")    

2. [Outline of Nucleus](https://github.com/NCBI-Hackathons/McDiff/blob/master/test_files/1.31.18_GFPP1_Hela_1min_002NuclMask.txt "Outline of Nucleus")    

3. [Region of Interest](https://github.com/NCBI-Hackathons/McDiff/blob/master/test_files/1.31.18_GFPP1_Hela_1min_002ROI.txt "Region of Intrest")     
 
# License
This project is licensed under the MIT License - see the LICENSE file for details

### Website

There is also a Docker image for hosting the main website. This should only be used for debug purposes.

  1. `$git clone https://github.com/NCBI-Hackathons/<this software>.git`
  2. `$cd Website`
  3. `$docker build --rm -t <this software>/website .`
  4. `$docker run -t -i <this software>/website`
  
