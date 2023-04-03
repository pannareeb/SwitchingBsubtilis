#Method3: biofilm level NuBacPDE model
import Pkg
Pkg.add("MethodOfLines")
Pkg.add("DomainSets")
Pkg.add("Symbolics")
using ModelingToolkit
using MethodOfLines, DomainSets
using Symbolics
using CSV

#Part 1: construct the NuBacPDE model 
@parameters t x #independent variables
@variables b(..) n(..) #dependent variables
#define derivatives
Dt = Differential(t)
Dx = Differential(x)
Dxx = Differential(x)^2
function founder(x, x0, xmax, bstart)
    if x > (xmax-x0)/2 && x < (xmax+x0)/2
        return bstart
    else
        return 0.
    end 
end
@register founder(x, x0,xmax,bstart)
#Space and time domains
xmax = 11.160121*2 #diameter of biofilm, based on data we collect
x0 = 3. # diameter of founding populations, based on data we collect
tmax = 92.6 #growing tome of biofilm, based on data we collect
bstart = 0.001 #(taken from paper), unitless
nstart = 1. #normalised to 1 (taken from paper), unitless
domains = [t ∈ Interval(0.0,tmax),x ∈ Interval(0.0,xmax)]
D_b = 0.075 #we will vary this to fit with the data
#exactly from the data; D_b = 0.001787256032 #Diffusion coeff of bacteria ; unit = mm2 h-1 
D_n =  3.36 #Diffusion coeff of nutrient (glycerol conc 0.1-0.5%); unit = mm2 h-1 
kg = 1.90 #maximum growth rate of bact; unit = h-1 
kd = 0.62 #maximum death rate of bact; unit = h-1 
K = 1.21 # Half–saturation nutrient for bacterial growth or death; no unit
θ = 0.1 #Fraction of nutrient released from dead bacteria; no unit
ω = 0.3 #Yield coefficient of bacteria (how much consumed nu become growth); no unit
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

#Part 2: PDE run and visualisation
#here we can use a normal ODE solver to solve our discretized PDE
sol = solve(prob, Tsit5(), saveat = 0.05)
sol_coarse = solve(prob, Tsit5(), saveat = 0.5)
#visualise
discrete_x = sol[x]
discrete_t = sol[t]
soln = sol[n(t,x)]
solb = sol[b(t,x)]
yannob = maximum(solb)/2.
yannon = maximum(soln)/2.
mycolt = range(HSL(colorant"red"), stop=HSL(colorant"blue"), length=length(discrete_t));
anim_bn = @animate for i in 1:4:length(discrete_t)
    p1 = plot(discrete_x,solb[i, 1:end], ylim = (minimum(solb),maximum(solb)),title="Bacteria at Time = $(round(discrete_t[i], digits = 1)) hr", ylab = "Bacteria number (b)",xlab = "Distance from one edge (x)",markershape = :circle,linewidth=2, legend = false, color = mycolt[i])
    annotate!([(19,yannob, ("b0 = $(bstart), x0 = $(x0)\nDb = $(round(D_b,digits = 2)), K = $(K)\nkd = $(kd), kg = $(kg)", 8, :blue))])
    p2 = plot(discrete_x,soln[i, 1:end], ylim = (minimum(soln),maximum(soln)),title= "Nutrient at Time = $(round(discrete_t[i], digits = 1)) hr", ylab = "Nutrient concentration (n)", xlab = "Distance from one edge (x)",markershape = :diamond,linewidth=2, legend = false, color = mycolt[i])
    annotate!([(19,yannon, ("n0 = $(nstart)\nDn = $(D_n), θ = $(θ)\nω = $(ω)", 8, :blue))])
    plot(p1,p2, layout = (2,1), size = (500,500))
end
gif(anim_bn,fps = 16)
#gif(anim_bn,fps = 16,"05ResultMoreDepth/anim_bnDbFromGuess05.gif")

#Part 3: PDE parameterisation
Db_scan = [0.055 0.06 0.065 0.07 0.075 0.08 0.085]
TimeSkip = Vector{Int64}(undef,141);
DistanceSkip = Vector{Int64}(undef,210);
bbfrontindex = Vector{Int64}(undef,length(TimeSkip));
bbfrontdis = Vector{Float64}(undef,length(TimeSkip));
#3.1 find front from b(t,x) from PDE for all Db's
bbfrontdisCollect = Matrix{Float64}(undef,length(TimeSkip),length(Db_scan));
bbfrontindexCollect = similar(bbfrontdisCollect);
for i in 1:length(Db_scan)
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
    D_b = Db_scan[i] #Diffusion coeff of bacteria ; unit = mm2 h-1
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
    dx = 0.5314343333 #define stepsize of x - what if it is 11.160121*2 /420 = 0.05314343333
    discretization = MOLFiniteDifference([x => dx], t, approx_order = 2) #t is continuous
    prob = discretize(pde_system,discretization) #form discretised PDE problem
    #here we can use ODE solver to solve
    sol = solve(prob, Tsit5(), saveat = 0.05)
    discrete_x = sol[x]
    discrete_t = sol[t]
    soln = sol[n(t,x)]
    solb = sol[b(t,x)]
    #CSV.write("solbLargeStepDb$(Int(D_b*1000))", DataFrame(solb, :auto))
    #CSV.write("solnLargeStepDb$(Int(D_b*1000))", DataFrame(soln, :auto))
    #define ToFindDbhalfb
    after_mid = ceil(Int,size(solb,2)/2)+1
    halfb = solb[11:13:1843,after_mid:(end-1)] #dt is row and dx is col
    ToFindDbhalfb = halfb' #dt is col and dx is row
    TimeSkip = discrete_t[11:13:1843]
    FrontbPDE = 0.1*maximum(solb)
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
    bbfrontdisCollect[:,i] = bbfrontdis
    bbfrontindexCollect[:,i] = bbfrontindex
    println("fin with $(i)")
end
#3.2 visualise
mycolDb =  range(HSL(colorant"red"), stop=HSL(colorant"purple"), length=length(Db_scan))
labelname = Vector{String255}(undef,length(Db_scan))
for i in 1:length(Db_scan)  
    labelname[i] = "Db = $(round.(Db_scan[i], digits = 3))" 
end
FrontPDEbplot = plot(TimeSkip,bbfrontdisCollect, ylab = "Distance from the centre (mm)", xlab = "Time (hrs)", lw = 2,palette = mycolDb, label = [labelname[1] labelname[2] labelname[3] labelname[4] labelname[5] labelname[6] labelname[7]], xticks = TimeSkip[1:14:end], yticks = round.(DistanceSkip,digits = 1),title = "Position of the front of PDE b(t,x), thre = 10%") 
#3.3 plot the real YR front again (refer to Method 1)
Ybfrontdis_tocf = Ybfrontdis./1000.
Rbfrontdis_tocf = Rbfrontdis./1000.
YRmain = plot(TimeInterval,Ybfrontdis_tocf, ylab = "Distance from centre (μm)", xlab = "Time (hrs)", label = "YFP-motile", color = :goldenrod, ls= :dash, marker = :circle, markersize = 1.5,markerstrokewidth = 0)
plot!(TimeInterval,Rbfrontdis_tocf, ylab = "Distance from centre (μm)", xlab = "Time (hrs)", label = "RFP-matrix", color = :red4, ls = :dash, marker = :circle, markersize = 1.5,markerstrokewidth = 0)
plot!(title = "Position of the front of motile and matrix bacteria, thre = 10%", xticks = TimeInterval[1:14:end],yticks = round.(DistanceSkip,digits = 1), titlefontsize = 10)
#3.4 plot both FrontPDEbplot and YRmain
Mixplot = 
plot(TimeSkip,bbfrontdisCollect, ylab = "Distance from the centre (mm)", xlab = "Time (hrs)", lw = 2,palette = mycolDb, label = [labelname[1] labelname[2] labelname[3] labelname[4] labelname[5] labelname[6] labelname[7]],yticks = round.(DistanceSkip,digits = 1), xticks = TimeSkip[1:14:end],title = "Position of the front of PDE b(t,x)")
plot!(TimeInterval,Ybfrontdis_tocf, ylab = "Distance from centre (μm)", xlab = "Time (hrs)", label = "YFP-motile", color = :goldenrod, ls= :dash, marker = :circle, markersize = 1.5,markerstrokewidth = 0)
plot!(TimeInterval,Rbfrontdis_tocf, ylab = "Distance from centre (μm)", xlab = "Time (hrs)", label = "RFP-matrix", color = :red4, ls = :dash, marker = :circle, markersize = 1.5,markerstrokewidth = 0)
title!("Position of the front of PDE b(t,x) and of the signals, thre = 10%",titlefontsize = 10)
#3.5 calculate loss function
Realfront = similar(Rbfrontdis_tocf)
for i in 1:length(TimeSkip)
    if (Rbfrontdis_tocf[i]>Ybfrontdis_tocf[i])
        Realfront[i] = Rbfrontdis_tocf[i] 
    else
        Realfront[i] = Ybfrontdis_tocf[i]
    end 
end
lossfunfront = similar(Db_scan)
for i in 1:length(Db_scan)
    lossfunfront[i] = sum(Realfront.-bbfrontdisCollect[:,i]).^2
end
SSplot = plot(Db_scan[1:6],lossfunfront[1:6], yscale=:log10, title = "Sum of square deviation from Realfront at different Db", yticks = [0,10,100,1000,10000], marker = :circle, label = false, titlefontsize = 12, xlab = "Db", ylab = "SS (lossfunfront)")


