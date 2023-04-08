# SwitchingBsubtilis

This is part of MSc Systems Biology course, University of Cambridge, 2022/2023, carried out over 6 months.

The file details are listed below

Method1local: Data processing 
1. #Part 1: basic data processing - import data into two data frames, clean, and visualisation 
2. #Part 2: focused data processing - plotting the movement of the front with several front thresholds and maximum signals of YFP or RFP

Method2local: single cell IRL-ODE model
1. #Part 1: create model, using Catalyst.jl
2. #Part 2: run both deterministic and stochastic ODE and visualisation
3. #Part 3: IRL-ODE parameterisation - scanning binding constants using analytical solution for finding RSS values and eigenvalues (with additional R codes for plotting bifurcation diagram) 

Method2HPC: Scanning the effect of noises on % of switch with either highR0 or lowR0 initial condition set and visualisation using histrogram. Used with HPC.  

Method3local: biofilm level NuBac-PDE model
1. #Part 1: create model, using MethodOfLines.jl
2. #Part 2: run deterministic PDE and visualisation
3. #Part 3: PDE parameterisation - scanning the diffusion coefficient of bacteria (Db) using loss function of the front movement difference

Method3HPC: Scanning Db (diffusion coefficient of bacteria in NuBac-PDE), nr and ny (fluorescence coefficients of RFP and YFP, respectively) and visualisation using heatmap. Used with HPC.

Method4HPC: Creating the lookup table for %Matrix phase, created by a given combination of the complexing constants. Note that to prevent time-limit exceeding, spliting cases into several files are recommendend

Method4local: Using the NuBac-PDE model prediction with the lookup table (created by stoIRL ODE runs), to form a loss function g(cf) that calculate the sum of squared deviation between the two. It is then optimised using Optim.jl and visualisation.



