#Method 1: Data processing steps

#Part 1: Data processing - general
#1.1 import data into two data frames (same dimensions)
spreadRFP = CSV.read("/Users/panareeboonyuen/SwitchingBsubtilis/04ExpDataSpread/rfp_tapa_data.csv", DataFrame)
spreadYFP = CSV.read("/Users/panareeboonyuen/SwitchingBsubtilis/04ExpDataSpread/yfp_hag_data.csv", DataFrame)
#selelct only the first column for the distance
distanceYR = Vector(spreadYFP[:,1])
#the rest of the data frames are the intensity, each col is at different timepoints, and each row is at different distance from the centre
spreadYFPm = Matrix(spreadYFP[:,2:142])
spreadRFPm = Matrix(spreadRFP[:,2:142])
#for further processing
Ndistance = size(spreadYFPm,1)
Ntimepoint = size(spreadYFPm,2)
colbydis = range(HSL(colorant"red"), stop=HSL(colorant"green"), length=Ndistance);
colbytime = range(HSL(colorant"red"), stop=HSL(colorant"green"), length=Ntimepoint);
TimeInterval = 0.0:0.66:92.4 #known from data collection description (collect every 40 mins)
Bandw= 26.5715 #the width of rings over which each reported intensity is averaged (unit is microns)
#1.2 visualise 1
spreadYFPplot = 
plot(distanceYR,spreadYFPm,legend = true, legendfontsize = 1, ylabel = "Intensity", xlabel = "Distance from centre (microns)", palette = colbytime, legendposition = :outerright, title = "YFP signal over time", titlefontsize = 10, ylim = (0,40))
scatter!(distanceYR,[spreadYFPm[:,1] spreadYFPm[:,end]],markersize = 1, marker = :cross, markerstrokewidth = 1, color = :black)
spreadRFPplot = 
plot(distanceYR,spreadRFPm,legend = true, legendfontsize = 1, ylabel = "Intensity", xlabel = "Distance from centre (microns)", palette = colbytime, legendposition = :outerright, title = "RFP signal over time", titlefontsize = 10, ylim = (0,40))
scatter!(distanceYR,[spreadRFPm[:,1] spreadRFPm[:,end]],markersize = 1, marker = :cross, markerstrokewidth = 1, color = :black)

#1.3 cleaning process
#1.3.1 : Normalised with the end values (the furthest from the centre) of corresponding timepoint
BasalYFP = Vector(spreadYFP[end,2:end])
BasalRFP = Vector(spreadRFP[end,2:end])
spreadYFPnorm1 = spreadYFPm.-BasalYFP'
spreadRFPnorm1 = spreadRFPm.-BasalRFP'
#1.3.2 visualise 2
spreadYFPnormplot = 
plot(distanceYR,spreadYFPnorm1,legend = true, legendfontsize = 1, ylabel = "Intensity", xlabel = "Distance from centre (microns)", palette = colbytime, legendposition = :outerright, title = "YFP signal (end values removed) over time", titlefontsize = 10, ylim = (0,40))
scatter!(distanceYR,[spreadYFPnorm1[:,1] spreadYFPnorm1[:,end]],markersize = 1, marker = :cross, markerstrokewidth = 1, color = :black)
spreadRFPnormplot = 
plot(distanceYR,spreadRFPnorm1,legend = true, legendfontsize = 1, ylabel = "Intensity", xlabel = "Distance from centre (microns)", palette = colbytime, legendposition = :outerright, title = "RFP signal (end values removed) over time", titlefontsize = 10, ylim = (0,40))
scatter!(distanceYR,[spreadRFPnorm1[:,1] spreadRFPnorm1[:,end]],markersize = 1, marker = :cross, markerstrokewidth = 1, color = :black)

#1.4: Normalised with bg noise common for all time points
#1.4.1 inspect at the 1st time slice when we know that RFP (matrix) should be 0, but it is not, so we normalise them further 
nn = 1
inc = 10
pbefore = plot()
plot!(distanceYR,spreadYFPnorm1[:,nn], label = "YFP at $(nn)", color = :goldenrod, ls = :dash)
plot!(distanceYR,spreadRFPnorm1[:,nn], label = "RFP at $(nn)", color = :red4, ls = :dash)
plot!(distanceYR,spreadYFPnorm1[:,nn+inc], label = "YFP at $(nn+inc)", color = :goldenrod2)
plot!(distanceYR,spreadRFPnorm1[:,nn+inc], label = "RFP at $(nn+inc)",color = :red3)
#1.4.2 decide the cutoffs and plot on the same graph
linetonormaliseRFP = (-1.2/10000).*distanceYR .+ 1.2 #this is of RFP
linetonormaliseYFP = 0.2 #this is of YFP

#1.4.3 therefore, we normalised the data again with these lines
spreadYFPnorm = spreadYFPnorm1 .- linetonormaliseYFP #remove noise below YFP cutoff
spreadRFPnorm = spreadRFPnorm1 .- linetonormaliseRFP #remove noise below RFP cutoff
#make negative signal to zero
spreadYFPnorm[spreadYFPnorm .<0] .= 0.
spreadRFPnorm[spreadRFPnorm .<0] .= 0.

#1.4.4 visualise the doubly normalised signals
spreadYFPnormplot2 = 
plot(distanceYR,spreadYFPnorm,legend = true, legendfontsize = 1, ylabel = "Intensity", xlabel = "Distance from centre (microns)", palette = colbytime, legendposition = :outerright, title = "YFP signal (end + base removed) over time", titlefontsize = 10, ylim = (0,40))
scatter!(distanceYR,[spreadYFPnorm[:,1] spreadYFPnorm[:,end]],markersize = 1, marker = :cross, markerstrokewidth = 1, color = :black)
spreadRFPnormplot2 = 
plot(distanceYR,spreadRFPnorm,legend = true, legendfontsize = 1, ylabel = "Intensity", xlabel = "Distance from centre (microns)", palette = colbytime, legendposition = :outerright, title = "RFP signal (end + base removed) over time", titlefontsize = 10, ylim = (0,40))
scatter!(distanceYR,[spreadRFPnorm[:,1] spreadRFPnorm[:,end]],markersize = 1, marker = :cross, markerstrokewidth = 1, color = :black)
#1.4.5 Summary all plots
ppraw = 
plot(spreadYFPplot,spreadRFPplot, layout = (2,1), size = (800,800),legend =false, plot_title = "Experimental Dual Reporter Data") 
ppbefore = 
plot(spreadYFPnormplot,spreadRFPnormplot, layout = (2,1), size = (800,800),legend =false, plot_title = "Experimental Dual Reporter Data, after end cut") 
ppafter = 
plot(spreadYFPnormplot2,spreadRFPnormplot2, layout = (2,1), size = (800,800),legend =false, plot_title = "Experimental Dual Reporter Data, after end and base cut") 
ppraw
pcompare = plot(ppraw ,ppafter, layout = (1,2), plot_title = "Comparison", title = " ")

#1.4.6 animate the doubly normalised signals over time
anim_doublenormYR = @animate for i in 1:Ntimepoint 
    plot(distanceYR,spreadYFPnorm[:, i], ylim = (0,30), lw = 2, label = "YFP - Motility", color = :goldenrod)
    plot!(distanceYR,spreadRFPnorm[:, i], ylim = (0,30),  lw = 2, label = "RFP - Matrix", color = :red4)
    plot!(title= "Doubly-normalised signal at time = $(round(TimeInterval[i], digits = 1)) hr", ylab = "Intensity", xlabel = "Distance from centre (microns)", legendposition = :outerbottomright)
end
gif(anim_doublenormYR,fps = 6)

#1.5 We can calculate the total signals for each time over all the distance (whole biofilm) for each signal
function allRingISum(Iv,Rv)
    eachRingA = 2*pi*Bandw .*Rv
    eachRingA[eachRingA.<0] .= 0
    withI = Iv .* eachRingA
    withI[withI .<0] .= 0
    allRingSum = sum(withI, dims = 1)
    return allRingSum
end
spreadYFPnorm_sum = allRingISum(spreadYFPnorm,distanceYR)
spreadRFPnorm_sum = allRingISum(spreadRFPnorm,distanceYR)
psigwhole = plot(TimeInterval,[spreadYFPnorm_sum', spreadRFPnorm_sum'],lw = 2, color = [:goldenrod :red4], label = ["YFP-motile" "RFP-matrix"], title = "whole biofilm signals", xlab = "Time (hrs)", ylab = "Intensity", xticks = TimeInterval[1:20:end])
#we can see the log phase of both population and see the stationary phases after around 40th hr for YFP and around 50th hr of RFP population

#Part 2: Data processing - plotting front and max
#2.1 plotting front
maxSigR = maximum(spreadRFPnorm)
maxSigY = maximum(spreadYFPnorm)
#plot Ybfrontdis(t) and Rbfrontdis(t)
bfront = Vector{Int64}(undef,Ntimepoint);
Ybfrontdis = Vector{Float64}(undef,Ntimepoint);
Rbfrontdis = Vector{Float64}(undef,Ntimepoint);
#2.1.1 front threshold is 1,1.5,2,5,or 10% of the max, below is 0.1 (10%)
for j in 1:Ntimepoint
    for i in reverse(1:Ndistance)
        if spreadYFPnorm[i,j] .> 0.1*maxSigY 
            bfront[j]= i
            Ybfrontdis[j] = distanceYR[bfront[j]]
        break 
        else 
            bfront[j]= 0
            Ybfrontdis[j] = 0
        end
    end
end
Ymain = plot(TimeInterval,Ybfrontdis, ylab = "Distance from one edge", xlab = "Timepoint", title = "Position of the front of Y bacteria", xticks = TimeInterval[1:20:end])
for j in 1:Ntimepoint
    for i in reverse(1:Ndistance)
        if spreadRFPnorm[i,j] .> 0.1*maxSigR
            bfront[j]= i
            Rbfrontdis[j] = distanceYR[bfront[j]]
        break
        else 
            bfront[j]= 0
            Rbfrontdis[j] = 0
        end
    end
end
Rmain = plot(TimeInterval,Rbfrontdis, ylab = "Distance from one edge", xlab = "Timepoint", title = "Position of the front of R bacteria", xticks = TimeInterval[1:20:end])
#2.1.2 plot front
pfrontmove = 
plot(TimeInterval,Ybfrontdis, ylab = "Distance from centre (μm)", xlab = "Time (hrs)", label = "YFP-motile", color = :goldenrod, xticks = TimeInterval[1:20:end],lw = 1, marker = :circle, ms = 2, markeralpha = 1)
plot!(TimeInterval,Rbfrontdis, ylab = "Distance from centre (microns)", xlab = "Time (hrs)", label = "RFP-matrix", color = :red4, xticks = TimeInterval[1:20:end],lw = 1, marker = :circle, ms = 2, markeralpha = 1, title = "The movement of the signal front, thre = 1%max")

#2.2 plotting max movement 
bmax = Vector{Int64}(undef,Ntimepoint);
Ybbmaxdis = Vector{Float64}(undef,Ntimepoint);
for j in 1:Ntimepoint
    for i in 1:Ndistance
        if spreadYFPnorm[i,j] ==  maximum(spreadYFPnorm[:,j])
            bmax[j]= i
            Ybbmaxdis[j] = distanceYR[bmax[j]]
        break
        end
    end
end
Rbbmaxdis = Vector{Float64}(undef,Ntimepoint);
for j in 1:Ntimepoint
    for i in 1:Ndistance
        if spreadRFPnorm[i,j] ==  maximum(spreadRFPnorm[:,j])
            bmax[j]= i
            Rbbmaxdis[j] = distanceYR[bmax[j]]
        break
        end
    end
end
pmaxsigmove = plot(TimeInterval,[Ybbmaxdis Rbbmaxdis], color = [:goldenrod :red4], title = "The movement of the maximum signal", xlab = "Time (hrs)", ylab = "Distance from the centre", label = ["max YFP" "max RFP"], lw = 1, marker = :circle, ms = 2, markeralpha = 1,xticks = TimeInterval[1:20:end])
vline!([TimeInterval[14]], label = "t = $(round(TimeInterval[14], digits = 1))", color = :goldenrod, ls = :dash)
vline!([TimeInterval[106]], label = "t = $(round(TimeInterval[106], digits = 1))", color = :red4, ls = :dash)
#2.3 plotting max values 
Rsignalmax = Vector{Float64}(undef,Ntimepoint);
for ti in 1:Ntimepoint
    Rsignalmax[ti] = maximum(spreadRFPnorm[:,ti])
end
Ysignalmax = Vector{Float64}(undef,Ntimepoint);
for ti in 1:Ntimepoint
    Ysignalmax[ti] = maximum(spreadYFPnorm[:,ti])
end
pmaxsigvalue = plot(TimeInterval,[Ysignalmax Rsignalmax], color = [:goldenrod :red4], label = ["YFP-motile" "RFP-matrix"], title = "Maximum signals", xlab = "Time (hrs)", ylab = "Intensity", xticks = TimeInterval[1:20:end], ylims = (0,30),lw = 1, marker = :circle, ms = 2, markeralpha = 1)
plot(pmaxsigmove,pmaxsigvalue, layout = (2,1))
