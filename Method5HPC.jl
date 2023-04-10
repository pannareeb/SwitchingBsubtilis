#tosendto HPC
import Pkg
Pkg.add("Catalyst")
Pkg.add("Latexify")
Pkg.add("DifferentialEquations")
Pkg.add("Plots")
Pkg.add("DataFrames")
Pkg.add("StatsPlots")
Pkg.add("ColorSchemes")
Pkg.add("Random")
Pkg.add("Distributions")
Pkg.add("JumpProcesses")
Pkg.add("CSV")
Pkg.add("DelimitedFiles")
Pkg.add("ModelingToolkit")
using ColorSchemes
using Random, Distributions
using Catalyst,Latexify, DifferentialEquations, JumpProcesses 
using Plots, DataFrames, StatsPlots
using CSV,DelimitedFiles

Pkg.add("MethodOfLines")
Pkg.add("DomainSets")
Pkg.add("Symbolics")
using ModelingToolkit
using MethodOfLines, DomainSets
using Symbolics

Pkg.add("Tables")
Pkg.add("Optim")
using Tables
using Optim

simplestIRL = @reaction_network begin
    α_0, 0 --> I #constitutive SinI production 
    β_0, 0 --> R #constitutive SinR production
    γ/(1+(R^2)), 0 --> L #SlrR production is a function of SinR level, hill coeff = 2
    1, I --> 0 #exponential dilution/degradation of SinI
    1, R --> 0 #exponential dilution/degradation of SinR
    1, L --> 0 #exponential dilution/degradation of SlrR
    δ_I, I + R --> 0 #complexing of SinR and SinI (C)
    δ_L, L + R --> 0 # complexing of SinR and SlrR (K)
end α_0 β_0 γ δ_I δ_L

#import the data and the model prediction of nutrients
halfn = Matrix(CSV.read("halfn_forSearch.csv", DataFrame))
#halfn = soln[11:13:1843,after_mid:(end-1)]'
spreadYFPm_g = Matrix(CSV.read("spreadYFPm_g.csv", DataFrame))
spreadRFPm_g = Matrix(CSV.read("spreadRFPm_g.csv", DataFrame))

#create the look-up table for nMotile and nMatrix
dI_span = 0:0.5:3
dL_span = 0:0.5:10
nMotileTable = Matrix{Float64}(undef, length(dI_span),length(dL_span) )
nMatrixTable = similar(nMotileTable)
for i in 1:length(dI_span)
    for l in 1:length(dL_span)
        nsim = 100
        u0 = [:I => 100, :R => 300, :L => 100] # we use either this High R0 ic
        #or u0 = [:I => 100, :R => 100, :L => 100] #this is low R0 ic
        tspan = (0.0,20.0)
        p = Dict(:α_0=> 85.,:β_0=>100., :γ=>125., :δ_I=>dI_span[i],:δ_L=>dL_span[l])
        Solu = Vector{Vector{Vector{Float64}}}(undef, nsim)
        Solt = Vector{Vector{Float64}}(undef,nsim)
        arrayrun = Array{Matrix{Float64}}(undef, nsim, 1)
        lastpoint = Vector{Vector{Float64}}(undef,nsim)
        for k in 1:nsim
            discreteIRLprob = DiscreteProblem(simplestIRL,u0,tspan,p)
            jumpIRLprob = JumpProcesses.JumpProblem(simplestIRL,discreteIRLprob, JumpProcesses.Direct())
            jsol = solve(jumpIRLprob, JumpProcesses.SSAStepper())
            Solu[k] = jsol.u #sol.u is Vector{Vector{Float64}}
            Solt[k] = jsol.t #sol.t is Vector{Float64} 
            arrayrun[k] = mapreduce(permutedims, vcat, Solu[k]) #take Solu_o[i] -> transpose each -> vcat them
            lastpoint[k] = arrayrun[k][end,:]
        end
        dflastpoint = mapreduce(permutedims, vcat, lastpoint)
        count_highR = sum(x -> x > 10, dflastpoint[:,2])
        nMotileTable[i,l] = count_highR
        nMatrixTable[i,l] = nsim - count_highR   
    end 
end

nMotileTable_df = DataFrame(nMotileTable, :auto)
nMatrixTable_df = DataFrame(nMatrixTable, :auto)
#when we use high R0 ic
CSV.write("percentMotileTableHighR0.csv", nMotileTable_df)
CSV.write("percentMatrixTableHighR0.csv", nMatrixTable_df)

