#Method 5.1 - Optimisation pipeline
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
Pkg.add("Polynomials")
using ColorSchemes
using Random, Distributions
using Catalyst,Latexify, DifferentialEquations, JumpProcesses 
using Plots, DataFrames, StatsPlots
using CSV,DelimitedFiles
using Polynomials

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


#import the model prediction from NuBac and the experimental data
halfn = Matrix(CSV.read("/Users/panareeboonyuen/SwitchingBsubtilis/09fminsearch/halfn_forSearch.csv", DataFrame))
halfb = Matrix(CSV.read("/Users/panareeboonyuen/SwitchingBsubtilis/09fminsearch/halfb_forSearch.csv", DataFrame))
spreadYFPnorm = CSV.read("spreadYFPnorm_double.csv", DataFrame)
spreadRFPnorm = CSV.read("spreadRFPnorm_double.csv", DataFrame)
spreadYFPm_g = Matrix(spreadYFPnorm[1:2:end,:])
spreadRFPm_g = Matrix(spreadRFPnorm[1:2:end,:])

#import the tables from stoIRL run which used 100 sims, each with timespan = 20 hr
#20 has proved to be good enough to distinguish between high Rss and low Rss case correctly
#but there are two types of tables for percentMotileTable
#one using a high R0 condition 
#the other using a low R0 condition

#the best range of tables we have got, wide and high resolution
percentMotileTable_HighR0 = Matrix(CSV.read("/Users/panareeboonyuen/SwitchingBsubtilis/10NewOptimisation/PmotfromHighR0Table.csv", DataFrame,header =false))./100
percentMotileTable_LowR0 = Matrix(CSV.read("/Users/panareeboonyuen/SwitchingBsubtilis/10NewOptimisation/PmotfromLowR0Table.csv", DataFrame,header =false))./100
# dI_span and dL_span
dI_span = 0:0.5:3
dL_span = 0:0.5:10 

# #small range, 0.5 resolution
# percentMotileTable_HighR0 = Matrix(CSV.read("/Users/panareeboonyuen/SwitchingBsubtilis/10NewOptimisation/MotileTableHighR0.csv", DataFrame))./100
# percentMotileTable_LowR0 = Matrix(CSV.read("/Users/panareeboonyuen/SwitchingBsubtilis/10NewOptimisation/MotileTableLowR0.csv", DataFrame))./100
# #shortened dI_span and dL_span
# dI_span = 0:0.5:2
# dL_span =0:0.5:4

# #wide range, 1.0 resolution
# percentMotileTable_HighR0 = Matrix(CSV.read("/Users/panareeboonyuen/SwitchingBsubtilis/10NewOptimisation/TestMotileTableHighR0.csv", DataFrame))./100
# percentMotileTable_LowR0 = Matrix(CSV.read("/Users/panareeboonyuen/SwitchingBsubtilis/10NewOptimisation/TestMotileTableLowR0.csv", DataFrame))./100
# #shortened dI_span and dL_span
# dI_span = 1:5
# dL_span = 0:10 


#visualise the percentage motile tables
motile_HighR0 = contourf(dL_span,dI_span,percentMotileTable_HighR0,color=:turbo,xlab = "SinR-SlrR complexing constant (δL)", ylab = "SinR-SinI complexing constant (δI)", title = "Percentage of a cell being in a motile phase\n(a high Rss) from 100 StoIRL runs with a high R0 ic", titlefontsize = 8, guidefontsize = 8, tick_direction = :out)
motile_LowR0 = contourf(dL_span,dI_span,percentMotileTable_LowR0,color=:turbo,xlab = "SinR-SlrR complexing constant (δL)", ylab = "SinR-SinI complexing constant (δI)", title = "Percentage of a cell being in a motile phase\n(a high Rss) from 100 StoIRL runs with a low R0 ic", titlefontsize = 8, guidefontsize = 8, tick_direction = :out)
plot(motile_HighR0,motile_LowR0, layout = (2,1), guidefontsize = 6)

#import for making distanceYRskip and TimeInterval
spreadRFP = CSV.read("/Users/panareeboonyuen/SwitchingBsubtilis/04ExpDataSpread/rfp_tapa_data.csv", DataFrame)
distanceYRskip = Vector(spreadRFP[:,1])[1:2:end]
TimeInterval = collect(0:0.66:92.4)


#define function "dILf(n)" that produced δI and δL from a given nutrient and coefficient
function dIL_f(cf,n)
    #cf is the coefficient vector and n is the nutrient
    dL = cf[1]+ cf[2]*n +cf[3]*n^2 +cf[4]*n^3 +cf[5]*n^4
    dI = cf[6]*dL
    return (dI,dL)
end 
#This link function will need to be in g for the purpose of using Optim.optimisation, but should be expressed globally as well 


#define a function for an input cf and give out g
#g(cf) is the loss function to be minimised by Optim.optimise using Nelder-Mead algo (default)
function g(cf_optim)
    #define function "dILf(n)" that produced δI and δL from a given nutrient and coefficient
    function dIL_f(cf,n)
        #cf is the coefficient vector and n is the nutrient
        dL = cf[1]+ cf[2]*n +cf[3]*n^2 +cf[4]*n^3 +cf[5]*n^4
        dI = cf[6]*dL
        return (dI,dL)
    end 
    #define a loss function for one data-prediction point 
    function gloss(cf) 
        glosstosum = similar(halfn)
        finalPMot_fromHR0 = similar(halfn)
        finalPMot_fromLR0 = similar(halfn)
        finalPMot = similar(halfn) 
        ny = 0.12
        nr = 0.08
        bMotile = 0
        bMatrix = 0
        δI_f = 0
        δL_f = 0
        mdi_ix = 0
        mdl_ix = 0
        for i in 1:size(glosstosum,1)
            for j in 1:size(glosstosum,2)
                #define experimental data
                YFP = spreadYFPm_g[i,j]
                RFP = spreadRFPm_g[i,j]
                #define model prediction - to be used later
                n = halfn[i,j]
                b = halfb[i,j]
                #define δI and δL from the link function, 
                #with the rounding to 0.5, given a predicted n
                δI_f = (ceil(dIL_f(cf,n)[1]) + floor(dIL_f(cf,n)[1]))/2
                δL_f = (ceil(dIL_f(cf,n)[2]) + floor(dIL_f(cf,n)[2]))/2
                # #with the rounding to whole number 
                # δI_f = floor(dIL_f(cf,n)[1])
                # δL_f = floor(dIL_f(cf,n)[2])
                #make  δI_f and δL_f stay within the range for search 
                if (δI_f < minimum(dI_span))
                    δI_f = minimum(dI_span)
                end
                if (δL_f < minimum(dI_span))
                    δL_f = minimum(dI_span)
                end
                if (δI_f > maximum(dI_span))
                    δI_f = maximum(dI_span)
                end
                if (δL_f > maximum(dL_span))
                    δL_f = maximum(dL_span)
                end
                #Searching in the spans to get the indices 
                mdi_ix = searchsorted(dI_span,δI_f)[1] 
                mdl_ix = searchsorted(dL_span,δL_f)[1]
                #use the indices to get the percentage of Motile SS
                #at the first time point (t=1 i.e. index j = 1), need no recalling of the previous state propbability 
                if j == 1 
                    finalPMot_fromHR0[i,j] = percentMotileTable_HighR0[mdi_ix,mdl_ix]
                    finalPMot_fromLR0[i,j] = 0
                    finalPMot[i,j] = finalPMot_fromHR0[i,j]+finalPMot_fromLR0[i,j]
                #at other time points, we calculated the compound probability (the table entry * how much in the previous state the cell ended up in SS corresponding to the type of the tables)
                else 
                    #contribution from HighR0 runs
                    finalPMot_fromHR0[i,j] = percentMotileTable_HighR0[mdi_ix,mdl_ix]*finalPMot[i,j-1]
                    #contribution from LowR0 runs
                    finalPMot_fromLR0[i,j] = percentMotileTable_LowR0[mdi_ix,mdl_ix]*(1-finalPMot[i,j-1])
                    finalPMot[i,j] = finalPMot_fromHR0[i,j]+finalPMot_fromLR0[i,j]
                end
                bMotile = b*finalPMot[i,j]
                bMatrix = b - bMotile
                glosstosum[i,j] = (bMotile-ny*YFP)^2 + (bMatrix-nr*RFP)^2
            end
            #println("fin - $(distanceYRskip[i])")
        end
        return  glosstosum
    end 
    sum_g = sum(gloss(cf_optim))
    return (sum_g)
end

#intial guess of cf vector - random
cf_0 = [1.,1.,1.,1.,1.,1.] .* rand(6)
test = g(cf_0) #loss function from this guess

#Optimisation part

#try testing different algorithms that need no gradiant
#but SimulatedAnnealing used even more than 5000 steps and still not succeed
#res = Optim.optimize(g, cf_0,SimulatedAnnealing(), Optim.Options(g_tol = 1e-8,iterations = 5000)) #actual optimisation process
#so we used the default algo - Nelder-Mead
resNM = Optim.optimize(g, cf_0)
answer = Optim.minimizer(resNM)
g(answer)
eqn = Polynomials.Polynomial(round.(answer[1:5],digits = 3))
plot(eqn, legendposition = :topright, xlim = (-1,1))

#we will do repeated optimisation as the answer seems inconsistant depending on the initial guesses
nOptround = 10
colbynOptround = range(HSL(colorant"red"), stop=HSL(colorant"purple"), length=nOptround)
cf_0_list = Vector{Vector{Float64}}(undef,nOptround)
cf_final_list = similar(cf_0_list)
for i in 1:nOptround
    cf_0_list[i] = rand(6)
    res = Optim.optimize(g,cf_0_list[i]) #use default: Convergence measures≤1.0e-08, Nelder-Mead algo
    cf_final_list[i] = Optim.minimizer(res)
    println("fin-$(i)")
end
cf_0_list #random initial cf for all 10 runs
cf_final_list #final cf for all 10 runs

#visualisation - with equations
plot(title = "fitted-δL",legendposition = :outerbottom, palette = colbynOptround,xlab = "Nutrient (n)", ylab = "dILf(n) output")
for i in 1:nOptround
    answer = cf_final_list[i]
    eqn = Polynomials.Polynomial(round.(answer[1:5],digits = 2))
    plot!(eqn, xlim = (0,1), lw = 2, marker = :hexagon, ms = 2, msw=0.5)
end
fitted_δLplot = current()
plot(title = "fitted-δI",legendposition = :outerbottom, palette = colbynOptround,xlab = "Nutrient (n)", ylab = "dILf(n) output")
for i in 1:nOptround
    answer = cf_final_list[i]
    eqn = Polynomials.Polynomial(round.(answer[1:5]*answer[6],digits = 2))
    plot!(eqn, xlim = (0,1), lw = 2, marker = :circle,ms = 2, msw=0.5)
end
fitted_δIplot = current()
plot(fitted_δLplot,fitted_δIplot, layout = (2,1), titlefontsize = 10, size = (500,900), guidefontsize = 10)
hcat(cf_0_list,cf_final_list)

#visualisation - clearer
n = 0:0.01:1 #nutrient in our system
dLfitted = similar(n)
dIfitted = similar(n)
dIfitted_array = Matrix{Float64}(undef,length(n) ,nOptround)
dLfitted_array = similar(dIfitted_array)

listofmodelID = collect(1:nOptround)
for i in 1:nOptround
    answer = cf_final_list[i]
    dLfitted = answer[1] .+ answer[2].*n .+ answer[3].*(n.^2) .+ answer[4].*(n.^3) .+ answer[5].*(n.^4) ;
    dIfitted = answer[6] .* dLfitted ;
    dLfitted_array[:,i] = dLfitted
    dIfitted_array[:,i] = dIfitted
    if sum(dIfitted.<0) >0 
        println("model $(i) produced some dI<0 ")
        filter!(e->e != i,listofmodelID)
    end
    if sum(dLfitted.<0) >0
        println("model $(i) produced some dL<0 ")
        filter!(e->e != i,listofmodelID)
    end
end
listofmodelID #the final model IDs we will use for finding Detstate

plot(n,dLfitted_array[:,listofmodelID],title = "fitted-δL", palette = colbynOptround,m = :hexagon, ms = 2, xlab = "Nutrient (n)", ylab = "dILf(n) output - δL",legendfontsize = 8, fgcolorlegend = false, bgcolorlegend = false, label = true, msw=0.1, legendposition = :bottomleft)
plot(n,dIfitted_array[:,listofmodelID],title = "fitted-δI", palette = colbynOptround,m = :hexagon, ms = 2, xlab = "Nutrient (n)", ylab = "dILf(n) output - δI",legendfontsize = 8, fgcolorlegend = false, bgcolorlegend = false, label = true, msw=0.1)






