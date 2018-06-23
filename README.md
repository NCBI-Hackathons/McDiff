![Packagist](https://img.shields.io/badge/python-3.6.5-orange.svg)
![Packagist](https://img.shields.io/packagist/l/doctrine/orm.svg)
![Packagist](https://img.shields.io/badge/Hackathon-in--progress-orange.svg)

# McDiff
A Monte Carlo Approach for Estimating Diffusion Coefficients   
This would be used by a cell biologist to get the coeficeint of free diffusion, and fraction of protein that accumulates at the region of interest within a cell.  
`$python3 run_sims.py mask_file ROI_file data_file`  
  
  Where `mask_file` is the file of the outline of the nucleus `ROI_file` is the file for region of interest and `data_file` is the file of the outputs from your FADD experiment.

# Abstract
Cell biologists can study the recruitment of DNA repair proteins to sites of DNA damage in live cells by using laser micro-irradiation to induce damage, which is known as Fluorescence Accumulation after DNA Damage (FADD). By monitoring the time-dependent accumulation of proteins at the sites of damage, or region of interest (ROI), the biologists then aim to calculate the coefficient of free diffusion (D) of each molecule of interest within the nucleus. This code simulates particles freely diffusing within a cell nucleus using a Monte Carlo model, becoming trapped at the ROI. Our code then seeks to optimize the best Diffusion constant (D) and fraction of protein that accumulates (F), and creates a heatmap for the best fit by using r-squared values comparing a fit to the data.


## What's the problem?
Before the implementation of this program, biologists would have to run either a mathematica or matlab scripts where each single simulation is run with manually entered coefficients. Quality of the fit to the experimental data was evaluated by visual inspection by the biologist. If the simulation was not aligned with the experimental data, then the biologist would change their coefficient values, then re-run the script and see if those values fit better, a tedious and labor-intensive process. Our new implementation runs simulations with multiple coefficients and calculates an r-squared error value that represents a best fit sample. Then the biologist can use that simulation, and be given the proper coefficients with the proper error.

# Keywords
 Free Diffusion Constant, Nucleus, Fluorescence Accumulation after DNA Damage (FADD)  

## Website (if applicable)



# How to use McDiff
`$python3 run_sims.py mask_file ROI_file data_file` 

  Where `mask_file` is the file of the oultine of the nucleus `ROI_file` is the file for region of intrest and `data_file` is the file of the outputs from your FADD experiment.  
  
# Output
You will get a .png of the best fit simulation over your experimental data, the D and F coeficients and the r-squared error value for the simulation that was best.
****Get Example output

# Assumption

# Software Workflow Diagram

![initial_flowchart_mcdiff](https://user-images.githubusercontent.com/23224399/41737510-beb45eb4-754c-11e8-816c-8720f1ae12e1.png)

## Methods

**Warm-up:**  

We started by getting comfortable with the experiments that were being performed in the nucleus. We then reviewed the existing Matlab and Mathematica scripts that originally ran this program, although with the lack of simulation to find optimal parameters.  

**Parsing/ Setup functions:** 

We uploaded the outline of the simulated nucleus. We then populated the nucleus with a gaussian distributed collection of simulated particles.  This is done by creating a square region around the nucleus. Then 36,000 particles are placed in the square; they are checked for localization within or without the nucleus and the first 12,000 points localized within the nucleus are used for the simulation. This implementation is ~10x faster than placing individual points and checking for their inclusion in the nucleus.   

**Simulation:**  

The simulation starts with all of the point in the nucleus, then the region of interest (ROI) is placed within the nucleus. We then calculate the percent of particles that are considered to be immobile, as for some cells, there are some proteins that are not free to move. Then within the ROI we calculate the percent that will not be read by the camera in the experiment because they will be bleached.  Once this setup is complete we can start to have the mobile particles move.  

The way that we calculate movement was given to us by **P-chem Guy** as he created the original program, and knows how the simulations should move. Each particle is moved either in the positive/negative x,y directions, with lengths of each having a random gaussian movement. 
The program then moves the particles for each timestep, and the number of particles within the ROI is recorded, with each timestep.  

We ran into a roadblock when the simulation of movement within the nucleus was taking a long time (~10 minutes). We found that the we could speed it up by almost 100x by using a vector library within the shapely library. Shapely is used to outline the nucleus and simulate the walls.  

Once we had a simulation of particles moving in the nucleus, we started work on plotting the nucleus simulation results.   

Then we worked on importing data from the experiments, and fitting simulated curves to the experimental curves. From there we could get an r-squared value and see which simulation gave us the best fitted curve with adjusted D and F values. This is done by running each simulation with different D’s and F’s and finding the best simulation with the most fit curve. 

**Output:**  

The user is then given a graph in .png form of the best fit graph, along with the D,F, and r-squared values for reporting.

# Installation options:

McDiff should be installed directly from Github.

### Installing McDiff from Github

### Step 1:
```
# Clone the repo
$git clone https://github.com/NCBI-Hackathons/McDiff.git
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
  
