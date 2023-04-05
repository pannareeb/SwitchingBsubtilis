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
Pkg.add("BifurcationKit")
Pkg.add("Setfield")
Pkg.add("SteadyStateDiffEq")
Pkg.add("PolynomialRoots")
Pkg.add("Polynomials")
using PolynomialRoots, Polynomials
using ColorSchemes
using Random, Distributions
using Catalyst,Latexify, DifferentialEquations, JumpProcesses 
using Plots, DataFrames, StatsPlots
using CSV,DelimitedFiles
using Setfield
using BifurcationKit, SteadyStateDiffEq

simplestIRL = @reaction_network begin
    α_0, 0 --> I #constitutive SinI production #take this out as we now include SpoOA-P (AP)
    β_0, 0 --> R #constitutive SinR production
    γ/(1+(R^2)), 0 --> L #SlrR production is a function of SinR level

    1, I --> 0 #exponential dilution/degradation of SinI
    1, R --> 0 #exponential dilution/degradation of SinR
    1, L --> 0 #exponential dilution/degradation of SlrR

    δ_I, I + R --> 0 #complexing of SinR and SinI (C)
    δ_L, L + R --> 0 # complexing of SinR and SlrR (K)
end α_0 β_0 γ δ_I δ_L
pspace = plot(xaxis = false, yaxis = false)

function TestNoise4(di,dl)
    #stochastic run
    nsim = 500
    u0HiR = Dict(:I => 100, :R => 300, :L => 100)
    u0LoR = Dict(:I => 100, :R => 100, :L => 100)
    tspan = (0.0,100.0)
    u0 = u0HiR
        p = Dict(:α_0=> 85,:β_0=>100, :γ=>125, :δ_I=>di,:δ_L=>dl)
        Solu = Vector{Vector{Vector{Float64}}}(undef, nsim)
        Solt = Vector{Vector{Float64}}(undef,nsim)
        for i in 1:nsim
            discreteIRLprob = DiscreteProblem(simplestIRL,u0,tspan,p)
            jumpIRLprob = JumpProcesses.JumpProblem(simplestIRL,discreteIRLprob, JumpProcesses.Direct())
            jsol = solve(jumpIRLprob, JumpProcesses.SSAStepper())
            Solu[i] = jsol.u #sol.u is Vector{Vector{Float64}}
            Solt[i] = jsol.t #sol.t is Vector{Float64}   
        end      
        arrayrun = Array{Matrix{Float64}}(undef, nsim, 1)
        lastpoint = Vector{Vector{Float64}}(undef,nsim)
        for i in 1:nsim
            arrayrun[i] = mapreduce(permutedims, vcat, Solu[i]) #take Solu_o[i] -> transpose each -> vcat them
            lastpoint[i] = arrayrun[i][end,:]
        end
        dflastpoint = mapreduce(permutedims, vcat, lastpoint)
        count_highR = sum(x -> x > 10, dflastpoint[:,2])
        h1 = histogram(dflastpoint[:,2], title = "I0, R0, L0 = $(u0[:I]),$(u0[:R]),$(u0[:L]) | δI = $(di), δL = $(dl)\n no. High Rss (Rss > 10) = $(count_highR)  ",label = false,titlefontsize =8,guidefontsize = 8,xtickfontsize=8, ytickfontsize = 8,normalize =:probability, ylim = (0.,1.0))
        vline!([10], lw = 2, label = false)
        #println("Fin-HiR")
    u0 = u0LoR
        p = Dict(:α_0=> 85,:β_0=>100, :γ=>125, :δ_I=>di,:δ_L=>dl)
        Solu = Vector{Vector{Vector{Float64}}}(undef, nsim)
        Solt = Vector{Vector{Float64}}(undef,nsim)
        for i in 1:nsim
            discreteIRLprob = DiscreteProblem(simplestIRL,u0,tspan,p)
            jumpIRLprob = JumpProcesses.JumpProblem(simplestIRL,discreteIRLprob, JumpProcesses.Direct())
            jsol = solve(jumpIRLprob, JumpProcesses.SSAStepper())
            Solu[i] = jsol.u #sol.u is Vector{Vector{Float64}}
            Solt[i] = jsol.t #sol.t is Vector{Float64}   
        end      
        arrayrun = Array{Matrix{Float64}}(undef, nsim, 1)
        lastpoint = Vector{Vector{Float64}}(undef,nsim)
        for i in 1:nsim
            arrayrun[i] = mapreduce(permutedims, vcat, Solu[i]) #take Solu_o[i] -> transpose each -> vcat them
            lastpoint[i] = arrayrun[i][end,:]
        end
        dflastpoint = mapreduce(permutedims, vcat, lastpoint)
        count_highR = sum(x -> x > 10, dflastpoint[:,2]) 
        h2 = histogram(dflastpoint[:,2], title = "I0, R0, L0 = $(u0[:I]),$(u0[:R]),$(u0[:L]) | δI = $(di), δL = $(dl)\n no. High Rss (Rss > 10) = $(count_highR) ",label = false,titlefontsize =8,guidefontsize = 8,xtickfontsize=8, ytickfontsize = 8,normalize =:probability, ylim = (0.,1.0))
        vline!([10], lw = 2, label = false)
    hh2 = plot(h1,h2, layout = (2,1))
    plot(pspace,hh2, layout = grid(1,2, widths = [0.1,0.9]))
end

TestNoise4(0.4,0.4)
savefig("Stari04l04_real.svg")
TestNoise4(0.4,2.)
savefig("Stari04l20_real.svg")

# TestNoise4(0.8,0.8)
# savefig("i08l08.png")
# TestNoise4(0.8,4.)
# savefig("i08l40.png")

TestNoise4(1.2,1.2)
savefig("Stari12l12_real.svg")
TestNoise4(1.2,6.)
savefig("Stari12l60_real.svg")

TestNoise4(2.,2.)
savefig("Stari20l20_real.svg")
TestNoise4(2.,10.)
savefig("Stari20l100_real.svg")

TestNoise4(4.,4.)
savefig("Stari40l40_real.svg")
TestNoise4(4.,20.)
savefig("Stari40l200_real.svg")

# TestNoise4(8.,8.)
# savefig("Stari80l80.svg")
# TestNoise4(8.,40.)
# savefig("i80l400.png")

TestNoise4(10.,10.)
savefig("Stari100l100_real.svg")
# TestNoise4(10.,20.)
# savefig("Stari100l200.svg")
# TestNoise4(10.,30.)
# savefig("i100l300.png")
# TestNoise4(10.,40.)
# savefig("i100l400.png")
# TestNoise4(10.,50.)
# savefig("i100l500.png")
