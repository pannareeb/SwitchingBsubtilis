# 🦠SwitchingBsubtilis

This is part of MSc Systems Biology course, University of Cambridge, 2022/2023, carried out over 6 months.
____________________________________________________ 

🦠The folders are:

1. RelevantDocument - Two presentations and the write-up of this project
2. ExperimentalData - biofilm data of YFP reporter (motile cells) and RFP reporter (matrix cells), described in the write-up
3. SupplementaryVideo - video of data or the simulation outputs, described in the write-up
4. DataSimOut-ForOptimise - normalised YFP and RFP signals and the simulation output of nutrient and bacterial density given by the fitted NuBacPDE, and the look-up tables from stoIRL run with specific complexing constants (i.e. the outputs of Method 1-4, which are the inputs of Method 5)
5. FinalSimOut  (Pending)
____________________________________________________

🦠The details of all code files are listed below.

(Note that the names of files follow the section "Materials and Methods" in the write-up of this work (In Doc folder). Hence, the Method 4 in the write-up, which is the parameterisation, was part of the file Method2 and Method3 in this Github). 
____________________________________________________
🧫Method1_local: Data processing 
1. #Part 1: basic data processing - import data into two data frames, clean, and visualisation 
2. #Part 2: focused data processing - plotting the movement of the front with several front thresholds and maximum signals of YFP or RFP
____________________________________________________
🧫Method2_local: single cell IRL-ODE model
1. #Part 1: create ODE model of SinI-SinR-SlrR (IRL), using Catalyst.jl
2. #Part 2: run both deterministic and stochastic ODE and visualisation
3. #Part 3: IRL-ODE parameterisation - scanning binding constants using analytical solution for finding RSS values and eigenvalues (with additional R codes for plotting bifurcation diagram) 

🧫Method2_HPC: TestNoise4(di,dl) is a function that test the effect of noises on %switch with either highR0 or lowR0 initial condition set and visualisation using histrogram. It can be used for scanning the effect of noises as a function of di and dl. Used with HPC.  
____________________________________________________
🧫Method3_local: biofilm level NuBac-PDE model
1. #Part 1: create PDE model, using MethodOfLines.jl
2. #Part 2: run deterministic PDE and visualisation
3. #Part 3: PDE parameterisation - scanning the diffusion coefficient of bacteria (Db) using loss function of the front movement difference, followed by the scanning of nr and ny (fluorescence coefficients of RFP and YFP, respectively) using all points in biofilms, and visualisation using heatmap.
____________________________________________________
(Method4 in the write-upwere part 3 of Method2local and Method3local files and also include the whole code in Method2HPC and Method3HPC)
____________________________________________________
🧫Method5_HPC: Creating the lookup tables for %Motile phase (and %Matrix phase) with a high R0 ic and low R0 ic, created by a given combination of the complexing constants. Note that to prevent time-limit exceeding, spliting cases into several files are recommendend. The range of each of the complexing constants to create the entries for was no need to be more than 10.

🧫Method5_1_local: 
1. Preparation of data and model output
2. Set up functions for link function and the loss function to optimise (a loss function g(cf) calculate the sum of squared deviation between the the NuBac-PDE model prediction with the lookup tables)
3. It is then optimised using Optim.jl and visualisation.
4. Intial checks of link function for relevant statistics.

🧫Method5_2_local: Using the parameterised link function from Method5.1 to predict the system behaviours of the growing biofilm 
1. Import and create basics
2. Set up main functions for testing of the candidate model (or a list of candidate models)
3. Runs of those functions inoutting model of interest
____________________________________________________

🧫slurm_example: slurm plain text file to be used for running julia in HPC.

