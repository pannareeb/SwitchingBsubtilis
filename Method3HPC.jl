#This is not actually in use anymore as the scanning od Db, nr, and ny fpor PDE model can be made faster by using a more corse dx scale.
#Tosendto HPC 
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
Pkg.add("MethodOfLines")
Pkg.add("DomainSets")
Pkg.add("Symbolics")
Pkg.add("ModelingToolkit")
using ModelingToolkit
using MethodOfLines, DomainSets
using Symbolics
using PolynomialRoots, Polynomials
using ColorSchemes
using Random, Distributions
using Catalyst,Latexify, DifferentialEquations, JumpProcesses 
using Plots, DataFrames, StatsPlots
using CSV,DelimitedFiles
using Setfield
using BifurcationKit, SteadyStateDiffEq

spreadRFPnorm = CSV.read("spreadRFPnormtoHPC.csv", DataFrame)
spreadRFPnorm = Matrix(spreadRFPnorm)
spreadYFPnorm = CSV.read("spreadYFPnormtoHPC.csv", DataFrame)
spreadYFPnorm = Matrix(spreadYFPnorm)
Ndistance = size(spreadYFPnorm,1)
Ntimepoint = size(spreadYFPnorm,2)

Db_scan = [0.04  0.06 0.07 0.075 0.08 0.10]
nr_scan = [0.01 0.05 0.10 0.25 0.5 1.0 1.25 2.0]
ny_scan = [0.01 0.05 0.10 0.25 0.5 1.0 1.25 2.0]
lossfuncollect = Array{Float64}(undef, length(Db_scan), length(nr_scan), length(ny_scan))
lossfun(halfb_t, IcProxy) = sum((halfb_t .- IcProxy).^2)
#halfb is the b(t,x) from nubacPDE from the centre; 
#from b(t,x), solb is 1853*421 (dt,dx)
function founder(x, x0, xmax, bstart)
    if x > (xmax-x0)/2 && x < (xmax+x0)/2
        return bstart
    else
        return 0.
    end 
end

#Run nubacPDE to get solb first
for i in 4#length(Db_scan)
    @parameters t x #independent variables
    @variables b(..) n(..) #dependent variables
    #define derivatives
    Dt = Differential(t)
    Dx = Differential(x)
    Dxx = Differential(x)^2
    @register founder(x, x0,xmax,bstart)
    #Space and time domains
    xmax = 11.160121*2 #diameter of biofilm, based on data we collect
    x0 = 3. # diameter of founding populations, based on data we collect
    tmax = 92.6 #based on data we collect
    bstart = 0.001 #(taken from paper)
    nstart = 1. #normalised to 1 (taken from paper)
    domains = [t ∈ Interval(0.0,tmax),x ∈ Interval(0.0,xmax)]
    D_n =  3.36 #Diffusion coeff of nutrient (glycerol conc 0.1-0.5%); unit = mm2 h-1 
    #D_b = 0.001787256032 #Diffusion coeff of bacteria ; unit = mm2 h-1 
    D_b = Db_scan[i]
    kg = 1.90 #maximum growth rate of bact; unit = h-1 
    kd = 0.62 #maximum death rate of bact; unit = h-1 
    K = 1.21 # Half–saturation nutrient for bacterial growth or death; no unit
    θ = 0.1 #Fraction of nutrient released from dead bacteria; no unit
    ω = 0.3 #Yield coefficient of bacteria (how much consumed nu become growth); no unit
    #2D PDE and #we discretize space, leave the time undiscretized #we use MethodOfLines
    eqb = Dt(b(t,x)) ~  D_b * Dxx(b(t,x)) + kg * b(t,x)*(n(t,x)/ (n(t,x)+K)) - kd * b(t,x)*(n(t,x)/ (n(t,x)+K)) 
    eqn = Dt(n(t,x)) ~  D_n * Dxx(n(t,x)) - ω*kg * b(t,x)*(n(t,x)/ (n(t,x)+K)) + θ*ω*kd * b(t,x)*(n(t,x)/ (n(t,x)+K)) 
    eqs = [eqb,eqn]
    ics = [b(0,x) ~ founder(x,x0,xmax,bstart), n(0,x)~ nstart]#Initial conditions
    bcs = [Dx(n(t,xmax))~ 0.,Dx(b(t,xmax))~ 0.,Dx(n(t,0))~ 0.,Dx(b(t,0))~ 0.] #and boundary conditions
    ic_bcs = vcat(ics, bcs)
    @named pde_system = PDESystem(eqs,ic_bcs,domains,[t,x],[b(t,x),n(t,x)])#build PDE system
    dx = 0.05314343333 #define stepsize of x - what if it is 11.160121*2 /420 = 0.05314343333
    discretization = MOLFiniteDifference([x => dx], t, approx_order = 2) #t is continuous
    prob = discretize(pde_system,discretization) #form discretised PDE problem
    #here we can use ODE solver to solve
    sol = solve(prob, Tsit5(), saveat = 0.05)
    discrete_x = sol[x]
    discrete_t = sol[t]
    soln = sol[n(t,x)]
    solb = sol[b(t,x)]
    #CSV.write("solb$(i)Db$(Int(D_b*1000)).csv", DataFrame(solb, :auto))
    #CSV.write("soln$(i)Db$(Int(D_b*1000)).csv", DataFrame(soln, :auto))
    #define ToFindDbhalfb
    after_mid = ceil(Int,size(solb,2)/2)+1
    halfb = solb[11:13:1843,after_mid:(end-1)]
    ToFindDbhalfb = halfb'
    #then compare the signals
    for j in 1:1#length(nr_scan)
        for k in 1:1#length(ny_scan)
            #sum up signals
            IcProxy = similar(spreadRFPnorm)
            nr = nr_scan[j]
            ny = ny_scan[k]
            for tp in 1:Ntimepoint
                for dt in 1:Ndistance
                    IcProxy[dt,tp] = nr*spreadRFPnorm[dt,tp] + ny*spreadYFPnorm[dt,tp]
                end
            end 
            skipdistance = Int(size(IcProxy,1)/size(ToFindDbhalfb,1))
            ToFindnrnyIcProxy = IcProxy[1:skipdistance:end,:]
            #visualise heatmap
            #halfn = soln[11:13:1843,after_mid:(end-1)]
            # TimeSkip = discrete_t[11:13:1843]
            # DistanceSkip = discrete_x[after_mid:(end-1)] .- xmax/2
            # ToFindDb = heatmap(TimeSkip,DistanceSkip,ToFindDbhalfb, xticks = TimeSkip[1:10:end], yticks = round.(DistanceSkip[1:20:end], digits = 1),tick_direction = :out, title = "b(t,x) from nubacPDE", ylab = "Distance unit from centre", xlab = "Time unit") 
            # ToFindnrny =heatmap(TimeSkip,DistanceSkip,ToFindnrnyIcProxy, xticks = TimeSkip[1:10:end], yticks = round.(DistanceSkip[1:20:end], digits = 1),tick_direction = :out,title = "the summation equation of Ir and Iy", ylab = "Distance unit from centre", xlab = "Time unit")
            # ComparebPDEvsSumSig = plot(ToFindDb,ToFindnrny, layout = (2,1),titlefontsize = 10, guidefontsize = 8, size = (600,600))
            # stringname = "hmDb$(i)$(Int(D_b*1000))_nr$(j)$(Int(nr*100))_ny$(k)$(Int(ny*100)).png"
            # png(ComparebPDEvsSumSig,stringname)
            lossfuncollect[i,j,k] = lossfun(ToFindDbhalfb,ToFindnrnyIcProxy)
        end
    end
end

for i in 1:length(Db_scan)
    CSV.write("lossfun$(i)Db$(Int(Db_scan[i]*1000)).csv", DataFrame(lossfuncollect[i,:,:], :auto))
end



