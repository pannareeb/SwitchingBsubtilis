#Method2: single cell ODE model

import Pkg
Pkg.add("Catalyst")
Pkg.add("Latexify")
Pkg.add("DifferentialEquations")
Pkg.add("JumpProcesses")
Pkg.add("Plots")
Pkg.add("DataFrames")
Pkg.add("StatsPlots")
Pkg.add("ColorSchemes")
Pkg.add("Random")
Pkg.add("Distributions")
Pkg.add("CSV")
Pkg.add("DelimitedFiles")
Pkg.add("LinearAlgebra")
Pkg.add("PolynomialRoots")
Pkg.add("Polynomials")
using Catalyst,Latexify, DifferentialEquations, JumpProcesses 
using Plots, DataFrames, StatsPlots
using ColorSchemes
using Random, Distributions
using CSV,DelimitedFiles
using LinearAlgebra
using PolynomialRoots, Polynomials
using Setfield, BifurcationKit

#Part 1:Create the IRL model
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
latexify(convert(ODESystem,simplestIRL))|> render
Graph(simplestIRL)

#Part 2: ODE run 
#2.1 Deterministic Run 
u0 = [:I => 100, :R => 200, :L => 100] #initial conditions of proteins
tspan = (0.0,5) #timespan set to ensure the steadt state is reached 
p = Dict(:α_0=> 85,:β_0=>100, :γ=>125, :δ_I=>5,:δ_L=>2)
myIRLprob = ODEProblem(simplestIRL,u0,tspan,p) #create ODE problem
sol = solve(myIRLprob) #solve sol.u is the I,R,L and sol.t is the time
sslevelODE = round.(sol.u[end], digits=3) #collect the ss level of each protein

#2.2 visulisation - Deterministic Run
#one time-evolution plot 
TimeEvoDet2 = plot(sol, xaxis = "Time", yaxis = "Concentration", w = 2, title = "Plot of simplest IRL system \n α0, β0, γ, δI, δL = $(p[:α_0]), $(p[:β_0]), $(p[:γ]), $(p[:δ_I]), $(p[:δ_L]) \n initial conc. = $(u0)\nSS level of I, R, L = $(sslevelODE)",legend=:topleft,titlefontsize =10, size = (600,500),guidefontsize = 15,xtickfontsize=12, ytickfontsize = 12, lw = 4, legendfontsize = 10)
TimeEvoDetzoom = plot(sol, xaxis = "Time", yaxis = "Concentration", w = 2, title = "Plot of simplest IRL system \n α0, β0, γ, δI, δL = $(p[:α_0]), $(p[:β_0]), $(p[:γ]), $(p[:δ_I]), $(p[:δ_L]) \n initial conc. = $(u0)\nSS level of I, R, L = $(sslevelODE)",legend=:topleft,titlefontsize =10, size = (600,500),guidefontsize = 15,xtickfontsize=12, ytickfontsize = 12, lw = 4, legendfontsize = 10, ylim = (0,50)) 
# one phase plane plot
dfall = mapreduce(permutedims, vcat, sol.u) #change the vector of the vector sol.u to a matrix
PhasePlanDetIRL = plot(title = "Phase plot - I or L and R \n δI = $(p[:δ_I]), δL = $(p[:δ_L])", xaxis = "R(t)", yaxis = "I(t) or L(t)", dfall[:,2], [dfall[:,1],dfall[:,3]], label = ["I(t)" "L(t)"], arrow = :head,lw = 2,  marker = :circle, markersize = 2, markerstrokewidth = 0.2, titlefontsize =10, color = [1 3], yticks = 0:10:100, xticks = 0:20:300)
PhasePlanDetIRLzoom = plot(title = "Phase plot - I or L and R \n δI = $(p[:δ_I]), δL = $(p[:δ_L])", xaxis = "R(t)", yaxis = "I(t) or L(t)", dfall[200:end,2], [dfall[200:end,1],dfall[200:end,3]], label = ["I(t)" "L(t)"], arrow = :head,lw = 2,  marker = :circle, markersize = 2, markerstrokewidth = 0.2, titlefontsize =10, color = [1 3])
PhasePlanDetIR = plot(dfall[:,1],dfall[:,2], title = "Phase plot - I and L \n δ_I = $(p[:δ_I]), δ_L = $(p[:δ_L])", xaxis = "I(t)", yaxis = "R(t)", label = false, arrow = :head,lw = 4, marker = :circle, markersize = 1)
PhasePlanDetRL = plot(title = "Phase plot - I and L \n δ_I = $(p[:δ_I]), δ_L = $(p[:δ_L])", xaxis = "R(t)", yaxis = "L(t)", dfall[:,2],dfall[:,3], label = false, arrow = :head,lw = 4,  marker = :circle, markersize = 1)
PhasePlanDetIL = plot(title = "Phase plot - I and L \n δ_I = $(p[:δ_I]), δ_L = $(p[:δ_L])", xaxis = "I(t)", yaxis = "L(t)", dfall[:,1],dfall[:,3], label = false, arrow = :head,lw = 4,  marker = :circle, markersize = 1)
plot(PhasePlanDetIR,PhasePlanDetRL,PhasePlanDetIL, layout = (1,3), title = "")
#Several traces in a phase plane plot
ILrange = [10]
Rrange = [0 0.1 0.5 1 2.5 5 7.5 10 15 20]
di_sp = 0:10
dl_sp = 0:10
ntrace = (length(ILrange),length(Rrange),length(ILrange))
solu_det = Array{Vector{Vector{Float64}}}(undef, ntrace)
solt_det = Array{Vector{Float64}}(undef,ntrace)
icslist = Array{Vector{Pair{Symbol, Float64}}}(undef,ntrace)
dfallphase = Array{Matrix{Float64}}(undef,ntrace)
collectHighRss =[]
for i in 1:lastindex(ILrange)
    for j in 1:lastindex(Rrange)
        for k in 1:lastindex(ILrange)
            icslist[i,j,k] = [:I => ILrange[i], :R => Rrange[j], :L => ILrange[k]]
            u0 =  icslist[i,j,k]#initial conditions of proteins
            tspan = (0.0,10.0) #timespan set to ensure the steadt state is reached 
            p = Dict(:α_0=> 85,:β_0=>100, :γ=>125, :δ_I=>5.,:δ_L=>2.)
            myIRLprob = ODEProblem(simplestIRL,u0,tspan,p) #create ODE problem
            sol = solve(myIRLprob)
            solu_det[i,j,k] = sol.u
            solt_det[i,j,k] = sol.t
            dfallphase[i,j,k]  = mapreduce(permutedims, vcat, sol.u)
            if sol.u[end][2] > 10
                push!(collectHighRss,[i,j,k])
            end
        end
    end
end
colphase = range(HSL(colorant"red"), stop=HSL(colorant"green"), length=prod(ntrace))
#Several traces in a time evolution
plot(palette =colphase )
for i in 1:lastindex(ILrange)
    for j in 1:lastindex(Rrange)
        for k in 1:lastindex(ILrange)
            if [i,j,k] in collectHighRss
                plot!(solt_det[i,j,k], dfallphase[i,j,k][:,2],label = true, lw = 2)
                
            else
                plot!(solt_det[i,j,k],dfallphase[i,j,k][:,2], label = true , lw = 2) #"$(icslist[i,j,k])"

            end
        end
    end
end
plot!(xlab = "Time (hrs)", ylab = "SinR (R) level", title = "SinR level, High Rss for $(lastindex(collectHighRss)) out of $(prod(ntrace)) runs", titlefontsize = 12, guidefontsize = 10, tickfontsize = 10, ylim = (0,20),xlim = (0,10), tickdirection = :out, legendposition = :outertopright, yticks = 0:2:16)
#Several traces in a 3D phase plane
plot()
for i in 1:lastindex(ILrange)
    for j in 1:lastindex(Rrange)
        for k in 1:lastindex(ILrange)
            if [i,j,k] in collectHighRss
                plot!(dfallphase[i,j,k][:,1],  dfallphase[i,j,k][:,3], dfallphase[i,j,k][:,2],label = false, color = :blue,arrow = :head,lw = 3)
                scatter!([ILrange[i]],[ILrange[k]],[Rrange[j]], mc = :blue, label = false, markershape = :square, markeralpha = 0.5)
            else
                plot!(dfallphase[i,j,k][:,1],  dfallphase[i,j,k][:,3], dfallphase[i,j,k][:,2], label = false , color = :red, arrow = :head) #"$(icslist[i,j,k])"
                scatter!([ILrange[i]],[ILrange[k]],[Rrange[j]],  mc = :red, label = false, markeralpha = 0.5)
            end
        end
    end
end
plot!(xlab = "SinI (I) level", zlab = "SinR (R) level", ylab = "SlrR (L) level",title = "SinI-SinR-SlrR, High Rss for $(lastindex(collectHighRss)) out of $(prod(ntrace)) runs", titlefontsize = 12, size = (600,600), guidefontsize = 8)
#Several traces in 2D L-R phase plane
plot()
for i in 1:lastindex(ILrange)
    for j in 1:lastindex(Rrange)
        for k in 1:lastindex(ILrange)
            if [i,j,k] in collectHighRss
                plot!(dfallphase[i,j,k][:,3], dfallphase[i,j,k][:,2], label = false, color = :blue,arrow = :closed, lw = 3)
                scatter!([ILrange[k]],[Rrange[j]], mc = :blue, label = false, markershape = :square, markeralpha = 0.5)
            else
                plot!(dfallphase[i,j,k][:,3], dfallphase[i,j,k][:,2], label = false , color = :red, arrow = :closed) #"$(icslist[i,j,k])"
                scatter!([ILrange[k]],[Rrange[j]], mc = :red, label = false, markeralpha = 0.5)
            end
        end
    end
end
plot!(xlab = "SlrR (L) level", ylab = "\nSinR (R) level",title = "SlrR-SinR, High Rss for $(lastindex(collectHighRss)) out of $(prod(ntrace)) runs", titlefontsize = 10, size = (400,600), yticks = Rrange)
#Several traces in 2D I-R phase plane
plot()
for i in 1:lastindex(ILrange)
    for j in 1:lastindex(Rrange)
        for k in 1:lastindex(ILrange)
            if [i,j,k] in collectHighRss
                plot!(dfallphase[i,j,k][:,1], dfallphase[i,j,k][:,2], label = false, color = :blue,arrow = :closed, lw = 3)
                scatter!([ILrange[i]],[Rrange[j]], mc = :blue, label = false, markershape = :square, markeralpha = 0.5)
            else
                plot!(dfallphase[i,j,k][:,1], dfallphase[i,j,k][:,2], label = false , color = :red, arrow = :closed) #"$(icslist[i,j,k])"
                scatter!([ILrange[i]],[Rrange[j]], mc = :red, label = false, markeralpha = 0.5)
            end
        end
    end
end
plot!(xlab = "SinI (I) level", ylab = "\nSinR (R) level",title = "SinI-SinR, High Rss for $(lastindex(collectHighRss)) out of $(prod(ntrace)) runs", titlefontsize = 10, size = (400,600), yticks = Rrange)
#Plot null clines
p = Dict(:α_0=> 85,:β_0=>100, :γ=>125, :δ_I=>2,:δ_L=>10)
function nullclinesIRL(pars,sR)
    a = pars[:α_0]
    b = pars[:β_0]
    ga = pars[:γ]
    di = pars[:δ_I]
    dl = pars[:δ_L]
    null_sI = a ./(1 .+di .*sR)
    null_sR = b ./ (1 .+ di.*(a ./(1 .+di .*sR)) .+ dl.*(ga ./((1 .+dl .*sR).*(1 .+sR.^2))))
    null_sL = ga ./((1 .+dl .*sR).*(1 .+sR.^2))
    return null_sI,null_sR,null_sL
end
sR = 0:1:70
nullclinesIRL(p,sR)
plot(sR,nullclinesIRL(p,sR)[1], ylab = "I", xlim = (0,80), ylim = (0,100), lw = 2)
plot!(nullclinesIRL(p,sR)[3],sR ,xlab = "L",xlim = (0,80), ylim = (0,100),lw = 2)

plot(sR,nullclinesIRL(p,sR)[2], ylab = "R", lw = 2)
plot!(nullclinesIRL(p,sR)[1],sR ,xlab = "I", lw = 2)
plot(sR,sR, ylab = "R")
plot!(nullclinesIRL(p,sR)[3],sR ,xlab = "L")
plot(nullclinesIRL(p,sR)[1])

#2.3 Stochastic Run 
nsim = 100
u0 = [:I => 100, :R => 100, :L => 100] 
tspan = (0.0,10.0)
p = Dict(:α_0=> 85,:β_0=>100, :γ=>125, :δ_I=>5,:δ_L=>2)
Solu = Vector{Vector{Vector{Float64}}}(undef, nsim)
Solt = Vector{Vector{Float64}}(undef,nsim)
for i in 1:nsim
    discreteIRLprob = DiscreteProblem(simplestIRL,u0,tspan,p)
    jumpIRLprob = JumpProcesses.JumpProblem(simplestIRL,discreteIRLprob, JumpProcesses.Direct())
    jsol = solve(jumpIRLprob, JumpProcesses.SSAStepper())
    Solu[i] = jsol.u #sol.u is Vector{Vector{Float64}}
    Solt[i] = jsol.t #sol.t is Vector{Float64} 
    println("Fin- $(i)")
end
#2.4 find statistics and visualisation - Stochastic Run
arrayrun = Array{Matrix{Float64}}(undef, nsim, 1)
lastpoint = Vector{Vector{Float64}}(undef,nsim)
avgtenlast = Matrix{Float64}(undef,nsim,3)
avg100last = similar(avgtenlast)
for i in 1:nsim
    arrayrun[i] = mapreduce(permutedims, vcat, Solu[i]) #take Solu_o[i] -> transpose each -> vcat them
    lastpoint[i] = arrayrun[i][end,:]
    avgtenlast[i,:] = mean(arrayrun[i][(end-9):end,:], dims = 1)
    avg100last[i,:] = mean(arrayrun[i][(end-99):end,:], dims = 1)
end
dflastpoint = mapreduce(permutedims, vcat, lastpoint)
ssIsd = round(std((dflastpoint[:,1])), digits = 3)
ssRsd = round(std((dflastpoint[:,2])), digits = 3)
ssLsd = round(std((dflastpoint[:,3])), digits = 3)
#plot R at steady state (RSS) of each simulation and count the number of runs that reach high RSS
co = 8
count_ZeroR = sum(x -> x == 0, dflastpoint[:,2])
count_highR = sum(x -> x >= co, dflastpoint[:,2])
count_lowR = sum(x -> x < co, dflastpoint[:,2])
histogram(dflastpoint[:,2], title = "% simulations vs Rss range |\n Initial conc: $(u0) |\n % High Rss (Rss > $(co)) = $(count_highR) ", xlabel = "Rss", ylabel = "% simulations", label = false, titlefontsize = 10)
vline!([co], label = "Cut-off")
#plot highlighted runs
simHighR = []
simLowR = []
for i in 1:length(dflastpoint[:,2])
    if dflastpoint[:,2][i] > co
        push!(simHighR, i)
    else
        push!(simLowR, i)
    end
end
maxSS = maximum(dflastpoint)
if (length(simHighR)>=2) & (length(simLowR)>=2)
    h1 = Int(simHighR[end]);
    h2 = Int(simHighR[end-1]);
    l1 = Int(simLowR[end]);
    l2 = Int(simLowR[end-1]);
    p6_1 = plot(palette =:Spectral_4, arrayrun[[h1;h2;l1;l2],1],layout= (3,1), title = ["SinI, ss = $(mean(dflastpoint[:,1]))+/- $(ssIsd)" "SinR, ss = $(mean(dflastpoint[:,2]))+/- $(ssRsd)" "SlrR, ss = $(mean(dflastpoint[:,3]))+/- $(ssLsd)" ], titlefontsize = 10)
    plot!(ylim = (0,maxSS), legend = false)
else
    p6_1 = plot(palette =:Spectral_4, arrayrun[1:4,1],layout= (3,1), title = ["SinI, $(u0[1]), ss = $(mean(dflastpoint[:,1]))+/- $(ssIsd)" "SinR, $(u0[2]), ss = $(mean(dflastpoint[:,2]))+/- $(ssRsd)" "SlrR, $(u0[3]), ss = $(mean(dflastpoint[:,3]))+/- $(ssLsd)" ])
end
#plot all HighR or LowR runs separately
plot(arrayrun[simHighR],layout= (3,1),palette =colbynsim,legend = false, ylim = (0,50), title = ["SinI, $(u0[1]), ss = $(mean(dflastpoint[:,1]))+/- $(ssIsd)" "SinR, $(u0[2]), ss = $(mean(dflastpoint[:,2]))+/- $(ssRsd)" "SlrR, $(u0[3]), ss = $(mean(dflastpoint[:,3]))+/- $(ssLsd)" ], titlefontsize = 8, plot_titlefontsize = 10,plot_title = "$(length(simHighR)) runs reached motile SS in $(tspan[2]) hrs",xlim = (5000,8000), tickdirection = :out)
plot(arrayrun[simLowR],layout= (3,1),palette =colbynsim,legend = false, ylim = (0,8), title = ["SinI, $(u0[1]), ss = $(mean(dflastpoint[:,1]))+/- $(ssIsd)" "SinR, $(u0[2]), ss = $(mean(dflastpoint[:,2]))+/- $(ssRsd)" "SlrR, $(u0[3]), ss = $(mean(dflastpoint[:,3]))+/- $(ssLsd)" ], titlefontsize = 8, plot_titlefontsize = 10,plot_title = "$(length(simLowR)) runs reached matrix SS in $(tspan[2]) hrs",xlim = (9000,10000), tickdirection = :out)
#plot many runs
runstoplot = 12:10:100
#timetoplot = range(0, stop = tspan[2], length = size(arrayrun[4,1],1))
colbynsim = range(HSL(colorant"red"), stop=HSL(colorant"purple"), length= length(runstoplot))
plot(arrayrun[runstoplot],layout= (3,1),palette =colbynsim,legend = false, ylim = (0,20), title = ["SinI, $(u0[1]), ss = $(mean(dflastpoint[:,1]))+/- $(ssIsd)" "SinR, $(u0[2]), ss = $(mean(dflastpoint[:,2]))+/- $(ssRsd)" "SlrR, $(u0[3]), ss = $(mean(dflastpoint[:,3]))+/- $(ssLsd)" ], titlefontsize = 8, plot_titlefontsize = 10,plot_title = "$(length(simLowR)) runs reached matrix SS in $(tspan[2]) hrs",xlim = (0,10000), tickdirection = :out)
scatter(avgtenlast[:,2], label = false, title = "$tspan", titlefontsize = 10)
hline!([8], label = "cut-off")
scatter(avg100last[:,2], label = false, title = "avg100last $tspan",titlefontsize = 10)
hline!([8], label = "cut-off", ylim = (0,10))
sum( x->x>=8,avg100last[:,2])
sum( x->x>=8,avgtenlast[:,2])
sum( x->x>=8,dflastpoint[:,2])


#Part 3: IRL-ODE parameterisation
#only try to change δI and δL using two-way sacnning
di_sp = 0:0.01:10
dl_sp = 0:0.01:10
#pars =Dict(:α_0=>85,:β_0=>100, :γ=>125, :δ_I=>2,:δ_L=>10) #default 
pars =Dict(:α_0=>85,:β_0=>100, :γ=>125, :δ_I=>2,:δ_L=>10) #default 
#3.1 Finding RSS values for a certain parameter range
est_Rss = Array{Vector{ComplexF64}}(undef,length(di_sp),length(dl_sp))
for i in 1:length(di_sp)
    for j in 1:length(dl_sp)
        pars =Dict(:α_0=>85,:β_0=>100, :γ=>125,:δ_I=>di_sp[i],:δ_L=>dl_sp[j]) 
        a = pars[:α_0]
        b = pars[:β_0]
        ga = pars[:γ]
        di = pars[:δ_I]
        dl = pars[:δ_L] 
        c0 = b
        c1 = di*b + dl*b -1 -di*a - dl*ga
        c2 = di*dl*b + b - di - dl - di*dl*a - di*dl*ga
        c3 = di*b + dl*b - di*dl -1 -di*a 
        c4 = di*dl*b - di -dl - di*dl*a
        c5 = -di*dl
        coeff = round.([ c0 c1 c2 c3 c4 c5],digits = 3)
        p = Polynomials.Polynomial([c0,c1,c2,c3,c4,c5])
        est_Rss[i,j]  = PolynomialRoots.roots(coeffs(p))
    end
end
#3.2 Plotting the number of real positive RSS
didl_PosRss_out = copy(est_Rss)
didl_RealRss_out = copy(est_Rss)
didl_nPosSS_out = Array{Int64}(undef,length(di_sp),length(dl_sp))
for i in 1:length(di_sp)
    for j in 1:length(dl_sp) 
        realR = filter(x -> abs(imag(x)) < 1e-5, est_Rss[i,j])
        real_PosR  = filter(x -> real(x) > 0, realR)
        didl_RealRss_out[i,j] = realR
        didl_PosRss_out[i,j] = real_PosR
        didl_nPosSS_out[i,j] = length(real_PosR)
    end
end
#3.3 Summary
est_Rss #Rss values
didl_RealRss_out #Rss with no or negligible imag part
didl_PosRss_out #Rss with no or negligible imag part and positive real part 
didl_nPosSS_out #number of real positive Rss
#3.4 visualise
#3.4.1 plot number of real positive Rss
mycmap = [:black, :grey]
HeatMapnPosSS = heatmap(dl_sp,di_sp,didl_nPosSS_out,size = (600,600),xlabel = "δL", ylabel = "δI", title = "Number of real positive Rss \n α0, β0, γ = $(pars[:α_0]), $(pars[:β_0]), $(pars[:γ])\n scanning δI over 0.0 - $(di_sp[end]) \nand δL over 0.0 - $(dl_sp[end])", c= cgrad(mycmap, 2, categorical=true, rev = false),guidefontsize = 15,xtickfontsize=10, ytickfontsize = 10, clims = (0,4), yticks = di_sp[1:100:end], xticks = dl_sp[1:100:end], tick_direction = :out)
vline!([2], label = false, color = :white)
hline!([2], label = false,color = :white)
vline!([1], label = false, color = :white)
hline!([1], label = false,color = :white)
vline!([4], label = false, color = :white)
hline!([4], label = false,color = :white)
vline!([5], label = false, color = :white)
hline!([5], label = false,color = :white)
#3.4.2 plot the High Rss
didl_maxRss = Array{Float64}(undef,length(di_sp),length(dl_sp))
didl_minRss = similar(didl_maxRss)
didl_monoRss = similar(didl_maxRss)
for i in 1:length(di_sp)
    for j in 1:length(dl_sp) 
        if (didl_nPosSS_out[i,j] == 3)
            didl_maxRss[i,j] = maximum(real(didl_PosRss_out[i,j])) 
            didl_minRss[i,j] = minimum(real(didl_PosRss_out[i,j])) 
            didl_monoRss[i,j] = -1
        else
            didl_maxRss[i,j] = -1
            didl_minRss[i,j] = -1
            didl_monoRss[i,j] = maximum(real(didl_PosRss_out[i,j]))
        end
    end
end

cupb = ceil(maximum(didl_maxRss[didl_maxRss .>0]))
clob = floor(minimum(didl_maxRss[didl_maxRss .>0])) 
colforhm = range(HSL(colorant"red"), stop=HSL(colorant"green"), length=500)
mycmap2 = vcat(:black,colforhm)
HeatMapHighRss = heatmap(dl_sp,di_sp,didl_maxRss, clims = (clob-1,cupb), color = mycmap2, size = (600,600),xlabel = "δL", ylabel = "δI", title = "The value of High Rss \n α0, β0, γ = $(pars[:α_0]), $(pars[:β_0]), $(pars[:γ])\n scanning δI over 0.0 - $(di_sp[end]) \nand δL over 0.0 - $(dl_sp[end])",guidefontsize = 15,xtickfontsize=10, ytickfontsize = 10, yticks = di_sp[1:100:end], xticks = dl_sp[1:100:end], tick_direction = :out)
contour!(dl_sp,di_sp,didl_maxRss,color = :white, lw = 0.5)

cupblowR = (maximum(didl_minRss[didl_minRss .>0]))
HeatMapLowRss = heatmap(dl_sp,di_sp,didl_minRss, color = mycmap2, clims = (0, cupblowR), size = (600,600),xlabel = "δL", ylabel = "δI", title = "The value of Low Rss \n α0, β0, γ = $(pars[:α_0]), $(pars[:β_0]), $(pars[:γ])\n scanning δI over 0.0 - $(di_sp[end]) \nand δL over 0.0 - $(dl_sp[end])",guidefontsize = 15,xtickfontsize=10, ytickfontsize = 10, yticks = di_sp[1:100:end], xticks = dl_sp[1:100:end], tick_direction = :out)
contour!(dl_sp,di_sp,didl_minRss,color = :white, lw = 0.5)

cupbmono = maximum(didl_monoRss)
clobmono = minimum(didl_monoRss[didl_monoRss .>0])
HeatMapmonoRss = heatmap(dl_sp,di_sp,didl_monoRss, color = mycmap2,size = (600,600),xlabel = "δL", ylabel = "δI", title = "The value of mono Rss \n α0, β0, γ = $(pars[:α_0]), $(pars[:β_0]), $(pars[:γ])\n scanning δI over 0.0 - $(di_sp[end]) \nand δL over 0.0 - $(dl_sp[end])",guidefontsize = 15,xtickfontsize=10, ytickfontsize = 10, yticks = di_sp[1:100:end], xticks = dl_sp[1:100:end], tick_direction = :out)
contour!(dl_sp,di_sp,didl_monoRss,color = :white, lw = 0.5)
HeatMapmonoRsszoom = heatmap(dl_sp,di_sp,didl_monoRss, color = mycmap2,size = (600,600),xlabel = "δL", ylabel = "δI", title = "The value of mono Rss \n α0, β0, γ = $(pars[:α_0]), $(pars[:β_0]), $(pars[:γ])\n scanning δI over 0.0 - $(di_sp[end]) \nand δL over 0.0 - $(dl_sp[end])",guidefontsize = 15,xtickfontsize=10, ytickfontsize = 10, yticks = di_sp[1:10:end], xticks = dl_sp[1:100:end], tick_direction = :out, ylims = (0,0.5), )
contour!(dl_sp,di_sp,didl_monoRss,color = :white, lw = 0.5)

#3.5 Analytical solution for eigenvalues 
#3.5.1 Find Eigenvalues and their sign for every fixed point
findEigen = function (Rss, pars)
    a = pars[:α_0]
    b = pars[:β_0]
    ga = pars[:γ]
    di = pars[:δ_I]
    dl = pars[:δ_L]
    Iss = a/(1 +di*Rss) 
    Lss = ga/((1 +dl*Rss)*(1+Rss^2))
    a11 = 1 - di*Rss
    a12 = -di*Iss
    a13 = 0
    a21 = -di*Rss
    a22 = -1-(di*Iss)-(dl*Lss)
    a23 = -dl*Rss
    a31 = 0
    a32 = ((-2*ga*Rss)/((1 +Rss^2)^2))-(dl*Lss) 
    a33 = -1-(dl*Rss)
    JacobM = [a11 a12 a13; a21 a22 a23; a31 a32 a33] 
    return eigvals(JacobM)
end
#due to a high computational demand, we reduced the resolution of di_sp and dl_sp
di_sp = 0:0.1:10
dl_sp = 0:0.1:10 #already rerun the didl_nPosSS_outto have the same dims
didl_eigenval_out = Array{Vector{ComplexF64}}(undef,length(di_sp),length(dl_sp)) #for monostable system
didl_eigenval_highR = similar(didl_eigenval_out) #for high Rss in a bistable system
didl_eigenval_unstable = similar(didl_eigenval_out) #for unstable Rss in a bistable system
didl_eigenval_lowR = similar(didl_eigenval_out) #for low Rss in a bistable system
for i in 1:length(di_sp)
    for j in 1:length(dl_sp)
        pars =Dict(:α_0=>85,:β_0=>100, :γ=>125, :δ_I=>di_sp[i],:δ_L=>dl_sp[j])
        if length((didl_PosRss_out)[i,j]) == 1
            Rss = real(didl_PosRss_out)[i,j][1] #get the real values of stationationary points  
            didl_eigenval_out[i,j]  = findEigen(Rss,pars)
            didl_eigenval_highR[i,j]  = [0,0,0]
            didl_eigenval_unstable[i,j]  = [0,0,0]
            didl_eigenval_lowR[i,j]  = [0,0,0]
            #print("for $(i), $(j), it is monostable")
        else
            didl_eigenval_out[i,j] = [0,0,0]
            highRss, unstableR, lowRss = real(didl_PosRss_out)[i,j]
            didl_eigenval_highR[i,j]   = findEigen(highRss,pars)
            didl_eigenval_unstable[i,j]  = findEigen(unstableR,pars)
            didl_eigenval_lowR[i,j]  = findEigen(lowRss,pars)
            #print("for $(i), $(j), it is bistable")
        end
    end
    println("fin for all dl at di = $(di_sp[i])")
end
#3.5.2 find the number of positive eigenvals of each fixed point
didl_nPosEigen_out = Array{Int64}(undef,length(di_sp),length(dl_sp))
didl_nPosEigen_highR = similar(didl_nPosEigen_out)
didl_nPosEigen_unstable = similar(didl_nPosEigen_out)
didl_nPosEigen_lowR = similar(didl_nPosEigen_out)
for i in 1:length(di_sp)
    for j in 1:length(dl_sp) 
        if length((didl_PosRss_out)[i,j]) == 1
            PosEigen = filter(x -> real(x)>0, didl_eigenval_out[i,j])
            didl_nPosEigen_out[i,j] = length(PosEigen) 
        else
            didl_nPosEigen_out[i,j] = -1
        end
        if length((didl_PosRss_out)[i,j]) == 3
            didl_nPosEigen_highR[i,j] = length(filter(x -> real(x)>0, didl_eigenval_highR[i,j]))
            didl_nPosEigen_unstable[i,j] = length(filter(x -> real(x)>0, didl_eigenval_unstable[i,j]))
            didl_nPosEigen_lowR[i,j] = length(filter(x -> real(x)>0, didl_eigenval_lowR[i,j]))
        else
            didl_nPosEigen_highR[i,j] = -1
            didl_nPosEigen_unstable[i,j] = -1
            didl_nPosEigen_lowR[i,j] = -1
        end
    end
end
nPosEigen1 = heatmap(dl_sp,di_sp,didl_nPosEigen_out, clims = (-1,1),xlabel = "δL", ylabel = "δI", title = "nPosEigen of mono Rss \n α0, β0, γ = $(pars[:α_0]), $(pars[:β_0]), $(pars[:γ])\n scanning δI over 0.0 - $(di_sp[end]) \nand δL over 0.0 - $(dl_sp[end])",guidefontsize = 15,xtickfontsize=10, ytickfontsize = 10,tick_direction = :out, size = (600,600))
nPosEigen2 = heatmap(dl_sp,di_sp,didl_nPosEigen_highR, clims = (-1,1),xlabel = "δL", ylabel = "δI", title = "nPosEigen at highRss \n α0, β0, γ = $(pars[:α_0]), $(pars[:β_0]), $(pars[:γ])\n scanning δI over 0.0 - $(di_sp[end]) \nand δL over 0.0 - $(dl_sp[end])",guidefontsize = 15,xtickfontsize=10, ytickfontsize = 10,tick_direction = :out, size = (600,600))
nPosEigen3 = heatmap(dl_sp,di_sp,didl_nPosEigen_unstable, clims = (-1,1),xlabel = "δL", ylabel = "δI", title = "nPosEigen at unstable Rss \n α0, β0, γ = $(pars[:α_0]), $(pars[:β_0]), $(pars[:γ])\n scanning δI over 0.0 - $(di_sp[end]) \nand δL over 0.0 - $(dl_sp[end])",guidefontsize = 15,xtickfontsize=10, ytickfontsize = 10,tick_direction = :out, size = (600,600))
nPosEigen4 = heatmap(dl_sp,di_sp,didl_nPosEigen_lowR, clims = (-1,1),xlabel = "δL", ylabel = "δI", title = "nPosEigen at lowRss \n α0, β0, γ = $(pars[:α_0]), $(pars[:β_0]), $(pars[:γ])\n scanning δI over 0.0 - $(di_sp[end]) \nand δL over 0.0 - $(dl_sp[end])",guidefontsize = 15,xtickfontsize=10, ytickfontsize = 10,tick_direction = :out, size = (600,600))
plotnPosEigen = plot(nPosEigen1,nPosEigen2,nPosEigen3,nPosEigen4, titlefontsize = 10)

#3.6 Bifurcation diagram in R (easier), and also in Julia
#3.6.1 it is more easily implementated in R using package deBif due to a better description 
#Here is the code for Bf plotting in R
# library(deBif)
# #Create
# state <- c(I = 100., R = 300. , L = 100.)
# parms <- c(a=85., b=100. , g=125. , di =2. , dl=10. )
# IRLsys <-function(t,state,parms){
#   with(as.list(c(state,parms)), {
#     dI = a - I - di*I*R 
#     dR = b - R - di*I*R - dl*L*R
#     dL = g/(1+R^2) - L - dl*L*R
    
#     return(list(c(dI, dR, dL)))
#   })
# }
# #implement
# bifurcation(IRLsys, state, parms)
# #I=2 R=15 L=1
# and I=60 R=1 L=40
#to find 2 branches when scanning dl

#3.6.2 in Julia - even though, we use the same initial values, the branches are not completed - so we have to ignore this method
# #step1- set bifur parameter and its span, and variable we want to plot 
# u0 = [:I => 2, :R => 15, :L => 1] # give bfhighR
# u0 = [:I => 60, :R => 1, :L => 40] #give bfunstableR when scanning δ_L from 2 to 10
# pbf =Dict(:α_0=> 85,:β_0=>100, :γ=>125, :δ_I=>2,:δ_L=>10)
# bif_par = :δ_L
# p_span = (2.,10.)
# plot_var = :R 
# #step2- start parameter phase-space with plotting at bif_par = p_span[1]
# p_bstart = copy(pbf)
# p_bstart[bif_par] = p_span[1]
# #step3- extract the ODE derivative function and its jacobian in a form that BifurcationKit can use
# oprob = ODEProblem(simplestIRL, u0, (0.0,0.0), p_bstart; jac = true)
# F = (u,p) -> oprob.f(u, p, 0)
# J = (u,p) -> oprob.f.jac(u, p, 0)
# #step4- tell the index in oprob.p of our bifurcation parameter, and the index in oprob.u0 of the variable we wish to plot
# # get δ_I and I as a symbolic variables
# @unpack α_0,β_0,γ,δ_I,δ_L,R,I,L = simplestIRL
# # find their indices in oprob.p and oprob.u0 respectively
# bif_idx  = findfirst(isequal(δ_L), parameters(simplestIRL))
# plot_idx = findfirst(isequal(R), species(simplestIRL))
# #step5- bundle the information we have compiled so far into a BifurcationProblem
# bprob = BifurcationProblem(F, oprob.u0, oprob.p, (@lens _[bif_idx]);
#                            recordFromSolution = (x, p) -> x[plot_idx], J = J)
# #step6- specify the input options for the pseudo-arclength continuation method (PACM) which produces the diagram
# bopts = ContinuationPar(dsmax = 0.05,          # Max arclength in PACM.
#                         dsmin = 1e-7,          # Min arclength in PACM.
#                         ds=0.001,              # Initial (positive) arclength in PACM.
#                         maxSteps = 1000000,     # Max number of steps.
#                         pMin = p_span[1],      # Min p-val (if hit, the method stops).
#                         pMax = p_span[2],      # Max p-val (if hit, the method stops).
#                         detectBifurcation = 3) # Value in {0,1,2,3} as least to most accurate method to calculate
# #step7- now ready to compute the bifurcation diagram:
# #gr(size = (700,600))
# bf = bifurcationdiagram(bprob, PALC(), 2, (args...) -> bopts)

# bfunstableR = plot(bf, xlabel = string(bif_par), ylabel = string(plot_var), title = "Bf - $(u0)\n α0, β0, γ, δI, δL = $(pbf[:α_0]), $(pbf[:β_0]), $(pbf[:γ]), $(pbf[:δ_I]), $(pbf[:δ_L])", size = (600,600),guidefontsize = 15,xtickfontsize=12, ytickfontsize = 12, legendfontsize = 12, ylim = (0,20), xlim = (0,10)) 
# bfhighR


