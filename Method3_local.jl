#Method3: biofilm level NuBacPDE model (part 3 -4 is moved to Appendix III in the write-up)
import Pkg
Pkg.add("MethodOfLines")
Pkg.add("DomainSets")
Pkg.add("Symbolics")
using ModelingToolkit
using MethodOfLines, DomainSets
using Symbolics
using CSV


#Part 1: the NuBacPDE model
#1.1  Construct state variables and parameters
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
tmax = 92.6 #growing tome of biofilm, based on data we collect
domains = [t ∈ Interval(0.0,tmax),x ∈ Interval(0.0,xmax)]
#pcf is for the modification of the parameters
pcf_1 = Dict(:bstart => 1, :nstart => 1, :D_b=> 1, :D_n => 1, :kg => 1, :kd =>1, :K => 1,:D_b_fin => 0.2)
pcf_2 = Dict(:bstart => 1, :nstart => 1, :D_b=> 1, :D_n => 1, :kg => 4, :kd =>10, :K => 0.8,:D_b_fin => 0.2) 
pcf_3 = Dict(:bstart => 400/4, :nstart => 4/4, :D_b=> 1/4, :D_n => 1/4, :kg => 1/4, :kd =>1/4, :K => 0.02/4,:D_b_fin => 0.03)
#final pcf
pcf = pcf_3
bstart = 0.001 * pcf[:bstart]  #(taken from paper), unitless
nstart = 1. *pcf[:nstart]   #normalised to 1 (taken from paper), unitless
D_b = pcf[:D_b_fin] #0.075 *pcf[:D_b]  #we will vary this to fit with the data
#exactly from the data; D_b = 0.001787256032 #Diffusion coeff of bacteria ; unit = mm2 h-1 
D_n = 3.36*pcf[:D_n]  #Diffusion coeff of nutrient (glycerol conc 0.1-0.5%); unit = mm2 h-1 
kg = 1.90 *pcf[:kg]  #maximum growth rate of bact; unit = h-1 
kd = 0.62 *pcf[:kd]  #maximum death rate of bact; unit = h-1 
K = 1.21 *pcf[:K]  # Half–saturation nutrient for bacterial growth or death; no unit
θ = 0.3 #Fraction of nutrient released from dead bacteria; no unit
ω = 0.1 #Yield coefficient of bacteria (how much consumed nu become growth); no unit
#1D PDE and #we discretize space, leave the time undiscretized using MethodOfLines
eqb = Dt(b(t,x)) ~  D_b * Dxx(b(t,x)) + kg * b(t,x)*(n(t,x)/ (n(t,x)+K)) - kd * b(t,x)*(n(t,x)/ (n(t,x)+K)) 
eqn = Dt(n(t,x)) ~  D_n * Dxx(n(t,x)) - ω*kg * b(t,x)*(n(t,x)/ (n(t,x)+K)) + θ*ω*kd * b(t,x)*(n(t,x)/ (n(t,x)+K)) 
eqs = [eqb,eqn]
ics = [b(0,x) ~ founder(x,x0,xmax,bstart), n(0,x)~ nstart]#Initial conditions
bcs = [Dx(n(t,xmax))~ 0.,Dx(b(t,xmax))~ 0.,Dx(n(t,0))~ 0.,Dx(b(t,0))~ 0.] #and boundary conditions
ic_bcs = vcat(ics, bcs)
@named pde_system = PDESystem(eqs,ic_bcs,domains,[t,x],[b(t,x),n(t,x)])#build PDE system
latexify(pde_system)|> render
dx = 0.05314343333 #define stepsize of x - what if it is 11.160121*2 /420 = 0.05314343333
discretization = MOLFiniteDifference([x => dx], t, approx_order = 2) #t is continuous
prob = discretize(pde_system,discretization) #form discretised PDE problem
#1.2 call solver
#here we can use a normal ODE solver to solve our discretized PDE
sol = solve(prob, Tsit5(), saveat = 0.05)
#writedlm("directsol.txt",sol)
#sol_coarse = solve(prob, Tsit5(), saveat = 0.5)
#visualise
discrete_x = sol[x]
discrete_t = sol[t]
soln = sol[n(t,x)]
solb = sol[b(t,x)]
after_mid = ceil(Int,size(solb,2)/2)+1
halfb = solb[11:13:1843,after_mid:(end-1)]'
halfn = soln[11:13:1843,after_mid:(end-1)]'
CSV.write("/Users/panareeboonyuen/SwitchingBsubtilis/11NewFitFront/FinalsolbDb.csv", DataFrame(solb, :auto))
CSV.write("/Users/panareeboonyuen/SwitchingBsubtilis/11NewFitFront/FinalsolnDb.csv", DataFrame(soln, :auto))
CSV.write("/Users/panareeboonyuen/SwitchingBsubtilis/11NewFitFront/FinalhalfbDb.csv", DataFrame(halfb, :auto))
CSV.write("/Users/panareeboonyuen/SwitchingBsubtilis/11NewFitFront/FinalhalfnDb.csv", DataFrame(halfn, :auto))
#Actually put them into a function runScannedPDE() for a convenient use in Part 3 in this file (PDE parameterisation)
function runScannedPDE(pcf,fin)
    if fin == 1
        D_b = pcf[:D_b_fin]
    else
        D_b = 0.075*pcf[:D_b]
    end
    xmax = 11.160121*2 #diameter of biofilm, based on data we collect
    x0 = 3. # diameter of founding populations, based on data we collect
    tmax = 92.6 #growing tome of biofilm, based on data we collect
    domains = [t ∈ Interval(0.0,tmax),x ∈ Interval(0.0,xmax)]
    #pcf = Dict(:bstart => 1, :nstart => 1, :D_b=> 1, :D_n => 1, :kg => 4, :kd =>10, :K => 0.8)
    bstart = 0.001 * pcf[:bstart]  #(taken from paper), unitless
    nstart = 1. *pcf[:nstart]   #normalised to 1 (taken from paper), unitless
     #0.075 *pcf[:D_b]  #we will vary this to fit with the data
    #exactly from the data; D_b = 0.001787256032 #Diffusion coeff of bacteria ; unit = mm2 h-1 
    D_n = 3.36*pcf[:D_n]  #Diffusion coeff of nutrient (glycerol conc 0.1-0.5%); unit = mm2 h-1 
    kg = 1.90 *pcf[:kg]  #maximum growth rate of bact; unit = h-1 
    kd = 0.62 *pcf[:kd]  #maximum death rate of bact; unit = h-1 
    K = 1.21 *pcf[:K]  # Half–saturation nutrient for bacterial growth or death; no unit
    θ = 0.3 #Fraction of nutrient released from dead bacteria; no unit
    ω = 0.1 #Yield coefficient of bacteria (how much consumed nu become growth); no unit
    #1D PDE and #we discretize space, leave the time undiscretized using MethodOfLines
    eqb = Dt(b(t,x)) ~  D_b * Dxx(b(t,x)) + kg * b(t,x)*(n(t,x)/ (n(t,x)+K)) - kd * b(t,x)*(n(t,x)/ (n(t,x)+K)) 
    eqn = Dt(n(t,x)) ~  D_n * Dxx(n(t,x)) - ω*kg * b(t,x)*(n(t,x)/ (n(t,x)+K)) + θ*ω*kd * b(t,x)*(n(t,x)/ (n(t,x)+K)) 
    eqs = [eqb,eqn]
    ics = [b(0,x) ~ founder(x,x0,xmax,bstart), n(0,x)~ nstart]#Initial conditions
    bcs = [Dx(n(t,xmax))~ 0.,Dx(b(t,xmax))~ 0.,Dx(n(t,0))~ 0.,Dx(b(t,0))~ 0.] #and boundary conditions
    ic_bcs = vcat(ics, bcs)
    @named pde_system = PDESystem(eqs,ic_bcs,domains,[t,x],[b(t,x),n(t,x)])#build PDE system
    #latexify(pde_system)|> render
    dx = 0.05314343333*10 #define stepsize of x - what if it is 11.160121*2 /420 = 0.05314343333
    discretization = MOLFiniteDifference([x => dx], t, approx_order = 2) #t is continuous
    prob = discretize(pde_system,discretization) #form discretised PDE problem

    #1.2
    #here we can use a normal ODE solver to solve our discretized PDE
    sol = solve(prob, Tsit5(), saveat = 0.05)
    #writedlm("directsol.txt",sol)
    #sol_coarse = solve(prob, Tsit5(), saveat = 0.5)
    #visualise
    discrete_x = sol[x]
    discrete_t = sol[t]
    soln = sol[n(t,x)]
    solb = sol[b(t,x)]
    return(discrete_x,discrete_t,soln,solb)
end
#But the Part 2 code below (2.1-2.3) are methods for a prediction from a single PDE model of interest

#Part 2: visualisation of a single model
#2.1 checking b with the whole biofilm signal
#from Method 1
spreadYFPnorm_sum = allRingISum(spreadYFPm_g,distanceYR[1:2:end])
spreadRFPnorm_sum = allRingISum(spreadRFPm_g,distanceYR[1:2:end])
#sum our model the same way as allRingISum
sumEachTime = sum(solb,dims=2)
sumEachTimeNew = sum(halfb,dims =1)
sumEachTime_n = sum(soln,dims=2)
comparebn = plot(discrete_t,sumEachTime,label = "b(t)",ylab = "b(t)",xlab = "Time (hrs)", lw = 2,legendposition = :topright)
plot!(twinx(),discrete_t, sumEachTime_n,label = "n(t)", color = :red,ylab = "n(t)", lw = 2,legendposition = :bottomright, title = "pcf = $(values(pcf))")
compare_b_wholesig = plot(TimeInterval,[spreadYFPnorm_sum',spreadRFPnorm_sum'],lw = 2, color = [:goldenrod :red4], label = ["YFP-Motile" "RFP-Matrix"], legendposition = :topleft,title = "Whole-biofilm signals vs Model prediction of bacteria density (b(t))\npcf = $(values(pcf))", xlab = "Time (hrs)", ylab = "Intensity", xticks = TimeInterval[1:20:end], titlefontsize = 8)
plot!(twinx(),discrete_t,sumEachTime, label = "b(t)", legendposition = :bottomright,lw = 2, ylab = "b(t)")
compare_b_sumsig = plot(TimeInterval,[0.35.*spreadYFPnorm_sum' .+ 0.20.*spreadRFPnorm_sum'],lw = 2, color = :green, label = "Scaled sum", legendposition = :topleft,title = "Whole-biofilm signals vs Model prediction of bacteria density (b(t))\npcf = $(values(pcf)), nr=0.4, ny=0.5", xlab = "Time (hrs)", ylab = "Intensity", xticks = TimeInterval[1:20:end], titlefontsize = 10)
plot!(TimeInterval,sumEachTimeNew', label = "b(t)", legendposition = :bottomright,lw = 2, ylab = "b(t)")

#2.2 make a gif of PDE model predictions
yannob = maximum(solb)/2.
yannon = maximum(soln)/2.
mycolt = range(HSL(colorant"red"), stop=HSL(colorant"blue"), length=length(discrete_t));
anim_bn = @animate for i in 1:8:length(discrete_t)
    p1 = plot(discrete_x .-xmax/2,solb[i, 1:end], ylim = (minimum(solb),maximum(solb)+1),title="Bacteria at Time = $(round(discrete_t[i], digits = 1)) hr", ylab = "Bacteria density (b)",xlab = "Distance from centre (mm)",markershape = :circle,lw=2, legend = false, color = mycolt[i], markersize = 1, markerstrokewidth = 0, xlim = (0,(discrete_x .-xmax/2)[end]))
    annotate!([(10,yannob, ("b0 = $(bstart)\nx0 = $(x0)\nDb = $(round(D_b,digits = 3))\nK = $(K)\nkd = $(kd)\nkg = $(kg)", 8, :blue))], yticks = 0:2.5:12.5,xticks = 0.:2.:12.0)
    p2 = plot(discrete_x .-xmax/2,soln[i, 1:end], ylim = (minimum(soln),maximum(soln)),title= "Nutrient at Time = $(round(discrete_t[i], digits = 1)) hr", ylab = "Nutrient level (n)", xlab = "Distance from centre (mm)",markershape = :diamond,lw=2, legend = false, color = mycolt[i], markersize = 1, markerstrokewidth = 0, xlim = (0,(discrete_x .-xmax/2)[end]))
    annotate!([(10,yannon, ("n0 = $(nstart)\nDn = $(D_n)\nθ = $(θ)\nω = $(ω)", 8, :blue))], xticks = 0.:2.:12.0)
    plot(p1,p2, layout = (3,1), size = (600,600))
end
gif(anim_bn,fps = 16)
#gif(anim_bn,fps = 16, "NuBacSim_bn.gif")
#2.3 visualisation of rate of growth
after_mid = ceil(Int,size(solb,2)/2)+1
halfb = solb[11:13:1843,after_mid:(end-1)]'
halfbdis = distanceYR[1:2:end]
colbytime_halfb = range(HSL(colorant"red"), stop=HSL(colorant"purple"), length=Int((size(halfb,1))))
labelN = Array{String63}(undef, lastindex(halfbdis[1:20:end]))
for i in 1:lastindex(halfbdis[1:20:end])
    labelN[i] = string(Int(round(halfbdis[1:20:end][i],sigdigits=2)))
end
plot(TimeInterval,halfb[1:20:end,:]', palette = colbytime_halfb[1:20:end], m = :circle, legendfontsize = 5, ms = 3, markerstrokewidth = 0, label = permutedims(labelN), ylabel = "Bacteria density (b(t,x))", xlabel = "Time (t, hrs)", leg_title = "Distance (x, μm)", legendtitlefontsize = 6, legend_title_font_halign = :right,legend_position = :best, fgcolorlegend = false, bgcolorlegend = false, guidefontsize = 10, xticks = TimeInterval[1:20:end]) 
#2.4 visualisation of nutrient from a PDE  run (halfn is the result)
function drawHalfn(halfn, t)
    TimeInterval = 0:0.66:92.4
    index = searchsortedfirst(TimeInterval,t)
    return halfn[:,index]
end
nu_plot = plot(distanceYR[1:2:end],drawHalfn(halfn,0), palette = :darktest, lw = 2, label = "t = 0")
plot!(distanceYR[1:2:end],drawHalfn(halfn,10), palette = :rainbow, lw = 2, label = "t = 10")
plot!(distanceYR[1:2:end],drawHalfn(halfn,15), lw = 2, label = "t = 15")
plot!(distanceYR[1:2:end],drawHalfn(halfn,25), lw = 2, label = "t = 25")
plot!(distanceYR[1:2:end],drawHalfn(halfn,40), lw = 2, label = "t = 40")
plot!(distanceYR[1:2:end],drawHalfn(halfn,55), lw = 2, label = "t = 55", legendposition = :outertopright, xlab = "Distance from centre (μm)", ylab = "Predicted nutrient level (n(t))")


#Part 3: PDE parameterisation (here we scan Db first)
#we construct 2 functions to find the front of biofilm from the expeimental data or from the PDE model
function FindFrontData(spreadRFPnorm,spreadYFPnorm,thre)
    maxSigR = maximum(spreadRFPnorm)
    maxSigY = maximum(spreadYFPnorm)
    Ntimepoint = size(spreadYFPnorm,2)
    Ndistance = size(spreadYFPnorm,1)
    bfront = Vector{Int64}(undef,Ntimepoint)
    Ybfrontdis = Vector{Float64}(undef,Ntimepoint)
    Rbfrontdis = Vector{Float64}(undef,Ntimepoint)
    #2.1.1 front threshold is 1,1.5,2,5,or 10% of the max, below is 0.1 (10%)
    for j in 1:Ntimepoint
        for i in reverse(1:Ndistance)
            if spreadYFPnorm[i,j] .> thre*maxSigY 
                bfront[j]= i
                Ybfrontdis[j] = distanceYR[i]
            break 
            else 
                bfront[j]= 0
                Ybfrontdis[j] = 0
            end
        end
    end
    for j in 1:Ntimepoint
        for i in reverse(1:Ndistance)
            if spreadRFPnorm[i,j] .>thre*maxSigR
                bfront[j]= i
                Rbfrontdis[j] = distanceYR[bfront[j]]
            break
            else 
                bfront[j]= 0
                Rbfrontdis[j] = 0
            end
        end
    end
    Ybfrontdis_tocf = Ybfrontdis./1000.
    Rbfrontdis_tocf = Rbfrontdis./1000.
    Realfrontdis = similar(Rbfrontdis_tocf)
    for i in 1:length(Realfrontdis)
        if (Rbfrontdis_tocf[i]>Ybfrontdis_tocf[i])
            Realfrontdis[i] = Rbfrontdis_tocf[i] 
        else
            Realfrontdis[i] = Ybfrontdis_tocf[i]
        end 
    end
    return (Ybfrontdis_tocf,Rbfrontdis_tocf,Realfrontdis)
end
function FindFrontModel(xmax,solb, discrete_t,discrete_x, thre)
    after_mid = ceil(Int,size(solb,2)/2)+1
    halfb = solb[11:13:1843,after_mid:(end-1)] #dt is row and dx is col
    ToFindDbhalfb = halfb' #dt is col and dx is row
    TimeSkip = discrete_t[11:13:1843]
    FrontbPDE = thre*maximum(solb)
    DistanceSkip = discrete_x[after_mid:(end-1)] .- xmax/2
    bbfrontindex = Vector{Int64}(undef,length(TimeSkip));
    bbfrontdis = Vector{Float64}(undef,length(TimeSkip));
    for j in 1:length(TimeSkip)
        for k in reverse(1:length(DistanceSkip))
            if ToFindDbhalfb[k,j] .> FrontbPDE
                bbfrontindex[j]= k
                bbfrontdis[j] = DistanceSkip[k]
            break
            else 
                bbfrontindex[j]= 0
                bbfrontdis[j] = 0
            end
        end
    end
    return (ToFindDbhalfb,TimeSkip,DistanceSkip,bbfrontindex,bbfrontdis)
end
Db_scan = 0.01:0.01:0.10
xmax = 11.160121*2
thre = 0.1
Realfront = FindFrontData(spreadRFPnorm,spreadYFPnorm,thre)
TimeSkip = Vector{Float64}(undef,141);
DistanceSkip = []#Vector{Float64}(undef,210);
bbfrontindex = []#Vector{Int64}(undef,length(TimeSkip));
bbfrontdis = []#Vector{Float64}(undef,length(TimeSkip));
#3.1 calculate loss function using front pf biofilm (difference between model and real front)
bbfrontdisCollect = Matrix{Float64}(undef,length(TimeSkip),length(Db_scan));
bbfrontindexCollect = similar(bbfrontdisCollect);
lossfunfront = Vector{Float64}(undef,length(Db_scan))
for i in 1:length(Db_scan)
    pcf = pcf_3
    pcf[:D_b_fin] = Db_scan[i]
    #Part 1: construct the NuBacPDE model 
    @parameters t x #independent variables
    @variables b(..) n(..) #dependent variables
    #define derivatives
    Dt = Differential(t)
    Dx = Differential(x)
    Dxx = Differential(x)^2
    @register founder(x, x0,xmax,bstart)
    #visualise
    resultPDE = runScannedPDE(pcf,1)
    discrete_x = resultPDE[1]
    discrete_t = resultPDE[2]
    soln = resultPDE[3]
    solb = resultPDE[4]
    #CSV.write("11NewFitFront/solb_newfitfront$(D_b).csv", DataFrame(solb, :auto))
    #CSV.write("11NewFitFront/soln_newfitfront$(D_b)).csv", DataFrame(soln, :auto))   
    Front = FindFrontModel(xmax,solb, discrete_t,discrete_x, thre)
    ToFindDbhalfb = Front[1]
    TimeSkip = Front[2]
    DistanceSkip = Front[3]
    bbfrontindex = Front[4]
    bbfrontdis = Front[5]
    bbfrontdisCollect[:,i] = bbfrontdis
    bbfrontindexCollect[:,i] = bbfrontindex
    println("fin with Db = $(Db_scan[i])")
    #calculate loss fun for front
    Realfront = FindFrontData(spreadRFPnorm,spreadYFPnorm,thre)[3]
    lossfunfront[i] = sum((Realfront.- bbfrontdis).^2)
end
#3.2 visualise front from PDE alone
mycolDb =  range(HSL(colorant"red"), stop=HSL(colorant"purple"), length=length(Db_scan))
labelname = Vector{String255}(undef,length(Db_scan));
for i in 1:length(Db_scan)  
    labelname[i] = "Db = $(round.(Db_scan[i], digits = 3))" 
end
FrontPDEbplot = plot(TimeSkip,bbfrontdisCollect, ylab = "Distance from the centre (mm)", xlab = "Time (hrs)", lw = 2,palette = mycolDb, label = permutedims(labelname), xticks = TimeSkip[1:14:end], yticks = round.(DistanceSkip,digits = 1),title = "Position of the front of PDE b(t,x), thre = $(thre)") 
#3.3 plot loss function
SSplot = plot(Db_scan,(lossfunfront),title = "Sum of square deviation from Realfront at different Db\n $(round.((lossfunfront)))", marker = :circle, titlefontsize = 8, xlab = "Db", ylab = "\nSS (lossfunfront)", xticks = Db_scan[1:1:end], tickfontsize = 8, label = false, ylim = (0,700)) 
#3.4 plot the real YR front again 
threh = thre
Ybfrontdis_tocf = FindFrontData(spreadRFPnorm,spreadYFPnorm,threh)[1]
Rbfrontdis_tocf = FindFrontData(spreadRFPnorm,spreadYFPnorm,threh)[2]
Realfront = FindFrontData(spreadRFPnorm,spreadYFPnorm,threh)[3]
YRmain = plot(TimeInterval,[Ybfrontdis_tocf Rbfrontdis_tocf], ylab = "Distance from centre (μm)", xlab = "Time (hrs)", label = ["YFP-motile" "RFP-matrix"], color = [:goldenrod :red4], ls= :dash, marker = :circle, markersize = 2,markerstrokewidth = 0)
plot!(TimeInterval,Realfront, label = "bfront", color = :green, ls = :dash,markerstrokewidth = 0, lw = 4, linealpha = 0.5,title = "Position of the front of motile and matrix bacteria, thre = $(threh)", xticks = TimeInterval[1:14:end],yticks = round.(DistanceSkip,digits = 1), titlefontsize = 10, ylim = (0,xmax/2))
#3.5 plot both FrontPDEbplot and YRmain
Mixplot =
plot(TimeSkip,bbfrontdisCollect, ylab = "Distance from the centre (mm)", xlab = "Time (hrs)", lw = 2,palette = mycolDb, label = permutedims(labelname),yticks = round.(DistanceSkip,digits = 1), xticks = TimeSkip[1:14:end],title = "Position of the front of PDE b(t,x)") 
plot!(TimeInterval,Realfront, label = "bfront", color = :green, ls = :dash,markerstrokewidth = 0, lw = 4, linealpha = 0.5)
title!("Position of the front of PDE b(t,x) and of the signals, thre = $thre %",titlefontsize = 10, legendfontsize = 5)
#substitute the best Db from here back to the final parameterised NuBac 
#then we will use this to optimise the fluorescence coefficient nr ny


#Part 4: After getting the best Db, we find the best fluorescence coefficient nr (for RFP) and ny (for YFP)
#4.1run PDE model with the best Db
OptDb = 0.03
pcf = pcf_3
pcf[:D_b_fin] = OptDb
@parameters t x #independent variables
@variables b(..) n(..) #dependent variables
#define derivatives
Dt = Differential(t)
Dx = Differential(x)
Dxx = Differential(x)^2
@register founder(x, x0,xmax,bstart)   
resultPDE = runScannedPDE(pcf,1)
discrete_x = resultPDE[1]
discrete_t = resultPDE[2]
soln = resultPDE[3]
solb = resultPDE[4] 
#4.2 Get bacterial density from the PDE model
#We ran FindFrontModel. We are not actually using only front, but this function can give us the formatted solb for point-to-point comparision with the data spreadRFP and spreadYFP.
Front = FindFrontModel(xmax,solb, discrete_t,discrete_x, thre)
ToFindnrny_halfb = Front[1] #this is to be compared with the sum of YFP and RFP signals
TimeSkip = Front[2]
DistanceSkip = Front[3]
bbfrontindex = Front[4]
bbfrontdis = Front[5]
println("fin with Db = $(pcf[:D_b_fin]), next find ptp loss function")
#4.3Creating the span of nr and ny we want to scan
nr_scan = 0.2:0.05:0.5
ny_scan = 0.2:0.05:0.6
lossfuncollect = Array{Float64}(undef, length(nr_scan), length(ny_scan))
#4.4 Scanning over all the combination
for j in 1:length(nr_scan)
    for k in 1:length(ny_scan)
        #sum up signals 
        ScaledSum = similar(spreadYFPm_g)
        nr = nr_scan[j]
        ny = ny_scan[k]
        ScaledSum = nr*spreadRFPm_g .+ ny*spreadYFPm_g
        #the dim of ToFindnrny_halfb and ToFindnrny_ScaledSum are not exactly the same, so we pick subset of ScaledSum (a larger matrix)
        skipdistance = Int(size(ScaledSum,1)/size(ToFindnrny_halfb,1))
        ToFindnrny_ScaledSum = ScaledSum[1:skipdistance:end,:]

        #visualise heatmap       
        plotHalfbwFront = heatmap(TimeSkip,DistanceSkip,ToFindnrny_halfb, xticks = TimeSkip[1:10:end], yticks = round.(DistanceSkip[1:2:end], digits = 1),tick_direction = :out, title = "Formatted b(t,x) from nubacPDE: Db=$(pcf[:D_b_fin]))", ylab = "Distance unit from centre", xlab = "Time (hrs)", size = (600,300), xlim = (0,92.6))
        #plot!(TimeSkip,bbfrontdis, lw = 2, color = :blue, label = "b(t,x) Front")      
        plotScaledsum =heatmap(TimeSkip,DistanceSkip,ToFindnrny_ScaledSum, xticks = TimeSkip[1:10:end], yticks = round.(DistanceSkip[1:2:end], digits = 1),tick_direction = :out,title = "Scaled summation of YFP and RFP:nr=$(nr), ny=$(ny)", ylab = "Distance unit from centre", xlab = "Time (hrs)", size = (600,300), xlim = (0,92.6))
        #plot!(TimeSkip,Realfront, lw = 2, color = :blue, label = "Expermental Front")
        
        ComparebPDEvsSumSig = plot(plotHalfbwFront,plotScaledsum, layout = (2,1),titlefontsize = 10, guidefontsize = 8, size = (600,600), plot_title ="$(pcf[:D_b_fin]))_nr$(nr)_ny$(ny)" )
        stringname = "4newhmDb$(pcf[:D_b_fin])_nr$(nr)_ny$(ny).svg"
        #savefig(ComparebPDEvsSumSig,stringname)
        #calculate using sqrt of sum of square deviation 
        lossfuncollect[j,k] = sqrt(sum((ToFindnrny_halfb.-ToFindnrny_ScaledSum).^2))
    end
end  
lossfuncollect
#4.5 find the optimised ny and nr pair
minLoss = "$(round.(minimum(lossfuncollect) ))"
optNr = 0
optNy = 0
for j in 1:length(nr_scan)
    for k in 1:length(ny_scan)
        if lossfuncollect[j,k] == minimum(lossfuncollect)
            println("$j $k  nr=$(nr_scan[j]),ny=$(ny_scan[k]) gives the minimum lossfun, using every point pair from the model withDb=$(pcf[:D_b_fin]) and the data $(lossfuncollect[j,k])")
            optNr = j
            optNy = k
        end
        if lossfuncollect[j,k] < minimum(lossfuncollect)+minimum(lossfuncollect)*0.01
            println("$j $k  nr=$(nr_scan[j]),ny=$(ny_scan[k]) gives the loss <180, $(lossfuncollect[j,k]) ratio = $(nr_scan[j]/ny_scan[k]) ")
        end
    end
end
#4.6 visualise the heatmap of loss function annotated with the best combination of nr and ny
heatmap(ny_scan,nr_scan,lossfuncollect, xlab = "ny", ylab = "nr",title = "nr = $nr_scan, ny = $ny_scan\nmost minimised loss = $(minLoss) at nr=$(nr_scan[optNr]), ny=$(ny_scan[optNy])", titlefontsize = 10, color = :turbo, xticks = ny_scan, yticks = nr_scan)
annotate!([(ny_scan[optNy],nr_scan[optNr],("Min",6, 0.0,:red))], clim = (175,325), tick_direction = :out)
annotate!([(ny_scan[7],nr_scan[15],("Opt",6, 0.0,:red))], clim = (175,325), tick_direction = :out, title = "smallest loss = $(minLoss), chosen loss = $(lossfuncollect[15,7])")
#Summary of parameters
#Db = 0.03, nr = 0.05, ny = 0.45
#However, this make RFP a really low compare to YFP
#so we go with other pair with a more equal nr and ny => (0.15,0.40) 

#4.7 write a function to visualise the Sum of scaled YFP and RFP, bacterial densities, and Ratios for a combination of nr and ny of interest 
function CalSumScaled(nr,ny)
    SumScaled2 = (ny*spreadYFPm_g) .+ (nr*spreadRFPm_g)
    ratioMo_data2 = (ny*spreadYFPm_g)./SumScaled2
    ratioMa_data2 = (nr*spreadRFPm_g)./SumScaled2
    j0 = heatmap(TimeSkip,distanceYR[1:2:end]./1000,SumScaled2, color = :turbo, title = "Sum of scaled YFP and RFP\nny=$ny nr=$nr", tick_direction = :out,titlefontsize = 10, xlab = "Time (hrs)", ylab = "Distance from centre (mm)", clim = (0,13))
    plot!(TimeSkip,Realfront, lw = 2, color = :limegreen, label = "Exp Front")
    j3 = heatmap(TimeSkip,distanceYR[1:2:end]./1000,ny*spreadYFPm_g, color = :turbo, title = "Motile density\nny=$ny nr=$nr", tick_direction = :out,titlefontsize = 10, xlab = "Time (hrs)", ylab = "Distance from centre (mm)", clim = (0,13))
    plot!(TimeSkip,Realfront, lw = 2, color = :limegreen, label = "Exp Front")
    j4 = heatmap(TimeSkip,distanceYR[1:2:end]./1000,nr*spreadRFPm_g, color = :turbo, title = "Matrix density\nny=$ny nr=$nr", tick_direction = :out,titlefontsize = 10, xlab = "Time (hrs)", ylab = "Distance from centre (mm)", clim = (0,13))
    plot!(TimeSkip,Realfront, lw = 2, color = :limegreen, label = "Exp Front")
    j1 = heatmap(TimeSkip,distanceYR[1:2:end]./1000,ratioMa_data2, color = :turbo, title = "Ratio of Matrix to total\nny=$ny nr=$nr", tick_direction = :out,titlefontsize = 10, xlab = "Time (hrs)", ylab = "Distance from centre (mm)", clim = (0,1))
    plot!(TimeSkip,Realfront, lw = 2, color = :limegreen, label = "Exp Front")
    j2 = heatmap(TimeSkip,distanceYR[1:2:end]./1000,ratioMo_data2, color = :turbo, title = "Ratio Motile to total\nny=$ny nr=$nr", tick_direction = :out, titlefontsize = 10, xlab = "Time (hrs)", ylab = "Distance from centre (mm)", clim = (0,1))
    plot!(TimeSkip,Realfront, lw = 2, color = :limegreen, label = "Exp Front")
    j5 = heatmap(TimeSkip,distanceYR[1:2:end]./1000,ratioMa_data2./ratioMo_data2, color = :turbo, title = "Ratio of Motile to Matrix\nny=$ny nr=$nr", tick_direction = :out,titlefontsize = 10, xlab = "Time (hrs)", ylab = "Distance from centre (mm)")
    plot!(TimeSkip,Realfront, lw = 2, color = :limegreen, label = "Exp Front")
    
    plot(j3,j4,j2,j1,j0,j5,size = (800,1200), layout = (3,2) )
    
end
CalSumScaled(0.05,0.45)
CalSumScaled(0.15,0.40) #pick this one


#Summary of parameters
#Db = 0.03, nr = 0.15, ny = 0.40
