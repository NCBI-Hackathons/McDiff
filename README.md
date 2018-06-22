# McDiff
A Monte Carlo Approach for Estimating Diffusion Coefficients

![initial_flowchart_mcdiff](https://user-images.githubusercontent.com/23224399/41737510-beb45eb4-754c-11e8-816c-8720f1ae12e1.png)

# Abstract
Cell biologists can study the recruitment of DNA repair proteins to sites of DNA damage in live cells by using laser micro-irradiation to induce damage, which is known as Fluorescence Accumulation after DNA Damage (FADD). By monitoring the time-dependent accumulation of proteins at the sites of damage, or region of interest (ROI), the biologists then aim to calculate the coefficient of free diffusion (D) of each molecule of interest within the nucleus. This code simulates particles freely diffusing within a cell nucleus using a Monte Carlo model, becoming trapped at the ROI. Our code then seeks to optimize the best Diffusion constant (D) and fraction of protein that accumulates (F), and creates a heatmap for the best fit by using r-squared values comparing a fit to the data.
Keywords
 Free Diffusion Constant, Nucleus, Fluorescence Accumulation after DNA Damage (FADD)

# Introduction
Before the implementation of this program, biologists would have to run either a mathematica or matlab scripts where each single simulation is run with manually entered coefficients. Quality of the fit to the experimental data was evaluated by visual inspection by the biologist. If the simulation was not aligned with the experimental data, then the biologist would change their coefficient values, then re-run the script and see if those values fit better, a tedious and labor-intensive process. Our new implementation runs simulations with multiple coefficients and calculates an r-squared error value that represents a best fit sample. Then the biologist can use that simulation, and be given the proper coefficients with the proper error.
