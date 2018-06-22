# McDiff
A Monte Carlo Approach for Estimating Diffusion Coefficients

# Abstract
Cell biologists can study the recruitment of DNA repair proteins to sites of DNA damage in live cells by using laser micro-irradiation to induce damage, which is known as Fluorescence Accumulation after DNA Damage (FADD). By monitoring the time-dependent accumulation of proteins at the sites of damage, or region of interest (ROI), the biologists then aim to calculate the coefficient of free diffusion (D) of each molecule of interest within the nucleus. This code simulates particles freely diffusing within a cell nucleus using a Monte Carlo model, becoming trapped at the ROI. Our code then seeks to optimize the best Diffusion constant (D) and fraction of protein that accumulates (F), and creates a heatmap for the best fit by using r-squared values comparing a fit to the data.
Keywords
 Free Diffusion Constant, Nucleus, Fluorescence Accumulation after DNA Damage (FADD)

# Introduction
Before the implementation of this program, biologists would have to run either a mathematica or matlab scripts where each single simulation is run with manually entered coefficients. Quality of the fit to the experimental data was evaluated by visual inspection by the biologist. If the simulation was not aligned with the experimental data, then the biologist would change their coefficient values, then re-run the script and see if those values fit better, a tedious and labor-intensive process. Our new implementation runs simulations with multiple coefficients and calculates an r-squared error value that represents a best fit sample. Then the biologist can use that simulation, and be given the proper coefficients with the proper error.

# Keywords
 Free Diffusion Constant, Nucleus, Fluorescence Accumulation after DNA Damage (FADD)

# GoodDoc -- A Template for Simple and Clear Documentation of Bioinformatics Code!

## Please cite our work -- here is the ICMJE Standard Citation:

### ...and a link to the DOI:

## Awesome Logo

### You can make a free DOI with zenodo <link>

## Website (if applicable)

## Intro statement

## What's the problem?

## Why should we solve it?

# What is <this software>?

Overview Diagram

# How to use <this software>

# Software Workflow Diagram

![initial_flowchart_mcdiff](https://user-images.githubusercontent.com/23224399/41737510-beb45eb4-754c-11e8-816c-8720f1ae12e1.png)

# File structure diagram 
#### _Define paths, variable names, etc_

# Installation options:

We provide two options for installing <this software>: Docker or directly from Github.

### Docker

The Docker image contains <this software> as well as a webserver and FTP server in case you want to deploy the FTP server. It does also contain a web server for testing the <this software> main website (but should only be used for debug purposes).

1. `docker pull ncbihackathons/<this software>` command to pull the image from the DockerHub
2. `docker run ncbihackathons/<this software>` Run the docker image from the master shell script
3. Edit the configuration files as below

### Installing <this software> from Github

1. `git clone https://github.com/NCBI-Hackathons/<this software>.git`
2. Edit the configuration files as below
3. `sh server/<this software>.sh` to test
4. Add cron job as required (to execute <this software>.sh script)

### Configuration

```Examples here```

# Testing

We tested four different tools with <this software>. They can be found in [server/tools/](server/tools/) . 

# Additional Functionality

### DockerFile

<this software> comes with a Dockerfile which can be used to build the Docker image.

  1. `git clone https://github.com/NCBI-Hackathons/<this software>.git`
  2. `cd server`
  3. `docker build --rm -t <this software>/<this software> .`
  4. `docker run -t -i <this software>/<this software>`
  
### Website

There is also a Docker image for hosting the main website. This should only be used for debug purposes.

  1. `git clone https://github.com/NCBI-Hackathons/<this software>.git`
  2. `cd Website`
  3. `docker build --rm -t <this software>/website .`
  4. `docker run -t -i <this software>/website`
  
