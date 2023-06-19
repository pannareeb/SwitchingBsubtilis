#Method 5.2 - using the optimised link function to find Detstate and bMotile and bMatrix

#Part 1: import and create basics
#1.1 model IRL
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
#1.2 result from 5.1
cf_final_list #the list of the optimised cf vectors for dILf(n) 
listofmodelID #the final model IDs we will use 
#1.3 ODE using either highR0 ic or LowR0 ic
function detstateHighR0(di,dl)
    u0 = Dict(:I => 100, :R => 300, :L => 100)
    tspan = (0.0,20.0) #timespan set to ensure the steady state is reached 
    p = Dict(:α_0=> 85,:β_0=>100, :γ=>125, :δ_I=>di,:δ_L=>dl)
    myIRLprob = ODEProblem(simplestIRL,u0,tspan,p) #create ODE problem
    sol = solve(myIRLprob) #solve sol.u is the I,R,L and sol.t is the time
    sslevelODE = round.(sol.u[end], digits=3) #collect the ss level of each protein
    return sslevelODE
end
function detstateLowR0(di,dl)
    u0 = Dict(:I => 100, :R => 100, :L => 100)
    tspan = (0.0,20.0) #timespan set to ensure the steadt state is reached 
    p = Dict(:α_0=> 85,:β_0=>100, :γ=>125, :δ_I=>di,:δ_L=>dl)
    myIRLprob = ODEProblem(simplestIRL,u0,tspan,p) #create ODE problem
    sol = solve(myIRLprob) #solve sol.u is the I,R,L and sol.t is the time
    sslevelODE = round.(sol.u[end], digits=3) #collect the ss level of each protein
    return sslevelODE
end
#1.4 define function "dILf(n)" that produced δI and δL from a given nutrient and coefficient
function dIL_f(cf,n)
    #cf is the coefficient vector and n is the nutrient
    dL = cf[1]+ cf[2]*n +cf[3]*n^2 +cf[4]*n^3 +cf[5]*n^4
    dI = cf[6]*dL
    return (dI,dL)
end 

#1.5 function that make dIdLTable for all time and distance points
function makedIdLTable(answer)
    dItable = Matrix{Float64}(undef, size(halfn))
    dLtable = similar(dItable)
    for i in 1:size(halfn,1) 
        for j in 1:size(halfn,2) 
            n = halfn[i,j]
            δI_f = dIL_f(answer,n)[1] 
            δL_f = dIL_f(answer,n)[2] 

            if (δI_f < minimum(dI_span))
                δI_f = minimum(dI_span) 

            end
            if (δL_f < minimum(dL_span))
                δL_f = minimum(dL_span)

            end
            if (δI_f > maximum(dI_span))
                δI_f = maximum(dI_span)

            end
            if (δL_f > maximum(dL_span))
                δL_f = maximum(dL_span)

            end
            δI_f = floor(δI_f)
            δL_f = floor(δL_f)
            dItable[i,j] = δI_f
            dLtable[i,j] = δL_f

        end
    end
    #plotting if wanted
    # heatmap(dLtable, title = "δL from the model $(nn)\n$(round.(cf_final_list[nn],digits = 2))", titlefontsize = 10)
    # savefig("10NewOptimisation/From10Rounds/dLtable_plot_$(nn).svg")
    # heatmap(dItable, title = "δI from the model $(nn)\n$(round.(cf_final_list[nn],digits = 2))", titlefontsize = 10)
    # savefig("10NewOptimisation/From10Rounds/dItable_plot_$(nn).svg")
    return (dItable,dLtable)
end
#1.6 function that run for the time t = 1 across all distance points
function intiRssfromHR0(answer)
    Rss_HR0_t1 = Vector{Float64}(undef,lastindex(distanceYRskip))
    for i in 1:lastindex(Rss_HR0_t1)
        dItable = makedIdLTable(answer)[1]
        dLtable = makedIdLTable(answer)[2]
        di = dItable[i,1]
        dl = dLtable[i,1]
        Rss_HR0_t1[i] = detstateHighR0(di,dl)[2]
    end
    return Rss_HR0_t1
end
#Reason: For each run, at first t - all distance points using highR0 [100,300,100], and then at t = 2, depending on t=1 outcome (Rss > 10 or < 10) 
#it may begin with [100,300,100] or [100,100,100], so for each detIRL run, we change u0, di, dl

#1.7 function that take two functions above and give out the whole model prediction for bMatrix(t,x) and bMotile(t,x)
function detRss_variedR0_allt(answer)
    Rss_variedR0_allt = Matrix{Float64}(undef,lastindex(distanceYRskip),lastindex(TimeInterval))
    Rss_variedR0_allt[:,1] = intiRssfromHR0(answer)
    dItable = makedIdLTable(answer)[1]
    dLtable = makedIdLTable(answer)[2]
    for i in 1:size(Rss_variedR0_allt,1)
        for j in 2:size(Rss_variedR0_allt,2)
            di = dItable[i,j]
            dl = dLtable[i,j]
            if Rss_variedR0_allt[i,j-1] >= 8
                #u0 = u0HiR
                Rss_variedR0_allt[i,j] = detstateHighR0(di,dl)[2]
            else 
                #u0 = u0LoR
                Rss_variedR0_allt[i,j] = detstateLowR0(di,dl)[2]
            end
        end
    end
    #plotting if wanted
    # Rss_variedR0_plot = heatmap(TimeInterval, distanceYRskip,Rss_variedR0_allt, xlab = "Time (hrs)", ylab = "Distance from centre (μm)", title = "Rss, from the model $(nn)")
    # savefig("10NewOptimisation/From10Rounds/Rss_variedR0_plot$(nn).svg")
    # DetState_plot = heatmap(TimeInterval, distanceYRskip,Rss_variedR0_allt.>10, xlab = "Time (hrs)", ylab = "Distance from centre (μm)", title = "DetState, from the model $(nn)\nyellow = Motile, black = Matrix", titlefontsize = 10)
    # savefig("10NewOptimisation/From10Rounds/DetState_plot$(nn).svg")
    return Rss_variedR0_allt
end

#Part 2: main functions for testing of the candidate model (or a list of candidate models)
#2.1 Deterministic results
function testonecf01(repbestcf,listornot)
    if listornot == 1
        rep_bestcf = repbestcf[1][5]
        println("take the model from round 5")
    else
        rep_bestcf = repbestcf
        println("test this given cf")
    end
    # makedIdLTable(rep_bestcf)[1] is dI
    # makedIdLTable(rep_bestcf)[2] is dL
    detRssRepbestmodel = detRss_variedR0_allt(rep_bestcf) #is the det Rss
    p1 = heatmap(TimeInterval, distanceYRskip./1000,makedIdLTable(rep_bestcf)[2], title = "δL from the model\n$(round.(rep_bestcf,digits = 2))", xlab = "Time (hrs)", ylab = "Distance from centre (μm)",titlefontsize = 10, color = :turbo)
    plot!(TimeSkip,Realfront, lw = 2, color = :limegreen, label = "Exp Front")
    p2 = heatmap(TimeInterval, distanceYRskip./1000,makedIdLTable(rep_bestcf)[1], title = "δI from the model\n$(round.(rep_bestcf,digits = 2))", xlab = "Time (hrs)", ylab = "Distance from centre (μm)",titlefontsize = 10, color = :turbo)
    plot!(TimeSkip,Realfront, lw = 2, color = :limegreen, label = "Exp Front")
    p3 = heatmap(TimeInterval, distanceYRskip./1000,detRssRepbestmodel, xlab = "Time (hrs)", ylab = "Distance from centre (μm)", title = "DetRss, from the model\n$(round.(rep_bestcf,digits = 2))", color = :turbo)
    plot!(TimeSkip,Realfront, lw = 2, color = :limegreen, label = "Exp Front")
    #p4 = heatmap(TimeInterval, distanceYRskip,detRssRepbestmodel.>10, xlab = "Time (hrs)", ylab = "Distance from centre (μm)", title = "DetState, from the model $(nn)\nyellow = Motile, black = Matrix", titlefontsize = 10)
    p4 = heatmap(TimeInterval, distanceYRskip./1000,halfn,xlab = "Time (hrs)", ylab = "Distance from centre (μm)",titlefontsize = 10, tick_direction = :out, title = "n(t,x) from PDE", color = :turbo)
    plot!(TimeSkip,Realfront, lw = 2, color = :limegreen, label = "Exp Front")
    pdet = plot(p1,p2,p3,p4, size = (800,800),tick_direction = :out,color = :turbo)
    return pdet
end
#2.2 Stochastic results
function testonecf02(repbestcf, listornot)
    if listornot == 1
        rep_bestcf = repbestcf[1][5]
        println("take the model from round 5")
    else
        rep_bestcf = repbestcf
    end
    bMotile_out = Matrix{Float64}(undef,lastindex(distanceYRskip),lastindex(TimeInterval))
    bMatrix_out = similar(bMotile_out)
    finalPMot_fromHR0 = similar(bMotile_out)
    finalPMot_fromLR0 = similar(bMotile_out)
    finalPMot = similar(bMotile_out) 
    dI_span = 0:1
    dL_span = 0:1 
    bMotile = 0
    bMatrix = 0
    δI_f = 0
    δL_f = 0
    mdi_ix = 0
    mdl_ix = 0
    for i in 1:lastindex(distanceYRskip)
        for j in 1:lastindex(TimeInterval)
            #define model prediction - to be used later
            n = halfn[i,j]
            b = halfb[i,j]
            #define δI and δL from the link function, 
            #with the rounding to 0.5, given a predicted n
            δI_f = ((dIL_f(rep_bestcf,n)[1]))#+ floor(dIL_f(cf,n)[1]))/2
            δL_f = ((dIL_f(rep_bestcf,n)[2]))# + floor(dIL_f(cf,n)[2]))/2
            # #with the rounding to whole number 
            # δI_f = floor(dIL_f(cf,n)[1])
            # δL_f = floor(dIL_f(cf,n)[2])
            #make  δI_f and δL_f stay within the range for search 
            if (δI_f < minimum(dI_span))
                δI_f = minimum(dI_span) 

            end
            if (δL_f < minimum(dL_span))
                δL_f = minimum(dL_span)

            end
            if (δI_f > maximum(dI_span))
                δI_f = maximum(dI_span)

            end
            if (δL_f > maximum(dL_span))
                δL_f = maximum(dL_span)

            end
            δI_f = floor(δI_f)
            δL_f = floor(δL_f)
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
            bMotile_out[i,j] = bMotile
            bMatrix_out[i,j] = bMatrix
        end
        #println("fin - $(distanceYRskip[i])")
    end
    ratio_bmotile = bMotile_out./(bMotile_out.+bMatrix_out)
    # heatmap(bMotile_out)
    # heatmap(bMatrix_out)
    # heatmap(ratio_bmotile)
    # heatmap(halfn)
    # heatmap(halfb)
    p5 = heatmap(TimeInterval, distanceYRskip./1000,bMotile_out,xlab = "Time (hrs)", ylab = "Distance from centre (mm)", title = "Density of motile bacteria, from the model\n$(round.(rep_bestcf,digits = 2))", titlefontsize = 10, color = :turbo,tick_direction = :out, clim = (0,13))
    plot!(TimeSkip,Realfront, lw = 2, color = :limegreen, label = "Exp Front")

    p6 = heatmap(TimeInterval, distanceYRskip./1000,bMatrix_out,xlab = "Time (hrs)", ylab = "Distance from centre (mm)", title = "Density of matrix bacteria, from the model\n$(round.(rep_bestcf,digits = 2))", titlefontsize = 10, color = :turbo,tick_direction = :out, clim = (0,13))
    plot!(TimeSkip,Realfront, lw = 2, color = :limegreen, label = "Exp Front")

    p7 = heatmap(TimeInterval, distanceYRskip./1000,ratio_bmotile,xlab = "Time (hrs)", ylab = "Distance from centre (mm)", title = "The ratio of Motile bacteria to total, from the model\n$(round.(rep_bestcf,digits = 2))", titlefontsize = 10, color = :turbo,tick_direction = :out, clim = (0.8,1))
    plot!(TimeSkip,Realfront, lw = 2, color = :limegreen, label = "Exp Front")

    p8 = heatmap(TimeInterval, distanceYRskip./1000,(1 .- ratio_bmotile),xlab = "Time (hrs)", ylab = "Distance from centre (mm)", title = "The ratio of Matrix bacteria to total, from the model\n$(round.(rep_bestcf,digits = 2))", titlefontsize = 10, color = :turbo,tick_direction = :out, clim = (0,0.2))
    plot!(TimeSkip,Realfront, lw = 2, color = :limegreen, label = "Exp Front")

    p9 = heatmap(TimeInterval, distanceYRskip./1000,halfb,xlab = "Time (hrs)", ylab = "Distance from centre (mm)",titlefontsize = 10, tick_direction = :out, color = :turbo, title = "b(t,x) from the model", clim = (0,13))
    plot!(TimeSkip,Realfront, lw = 2, color = :limegreen, label = "Exp Front")

    p10 = heatmap(TimeInterval, distanceYRskip./1000, ratio_bmotile ./(1 .- ratio_bmotile) ,xlab = "Time (hrs)", ylab = "Distance from centre (mm)",titlefontsize = 10, tick_direction = :out, color = :turbo, title = "Ratio of Motile to Matrix, from the model")    
    plot!(TimeSkip,Realfront, lw = 2, color = :limegreen, label = "Exp Front")

    psto = plot(p5,p6,p7,p8,p9,p10,tick_direction = :out,size = (800,1200), layout = (3,2))
    return psto
end

#2.3 common stats
function testonecf03(repbestcf, listornot)
    if listornot == 1
        rep_bestcf = repbestcf[1][5]
        println("take the model from round 5")
    else
        rep_bestcf = repbestcf
    end
    out = gloss_global(rep_bestcf)
    minall = round(minimum(out[3]),digits = 1)
    maxall = round(maximum(out[3]),digits = 1)
    minYFP = round(minimum(out[4]),digits = 1)
    minRFP = round(minimum(out[5]),digits = 1)
    maxYFP = round(maximum(out[4]),digits = 1)
    maxRFP = round(maximum(out[5]),digits = 1)
    minYtoall = round(minimum(out[4]./out[3]),digits = 1)
    maxYtoall = round(maximum(out[4]./out[3]),digits = 1)
    h1 = heatmap(TimeInterval,distanceYRskip,out[3], title = "Loss matrix from YFP and RFP,\n($minall-$maxall)", titlefontsize = 10) #loss matrix
    h2 = heatmap(TimeInterval,distanceYRskip,out[4], title = "Loss matrix only from YFP,\n($minYFP-$maxYFP)", titlefontsize = 10) #loss matrix from YFP
    h3 = heatmap(TimeInterval,distanceYRskip,out[5], title = "Loss matrix only from RFP,\n($minRFP-$maxRFP)", titlefontsize = 10) #loss matrix from RFP
    h4 = heatmap(TimeInterval,distanceYRskip,out[4]./out[3], title = "Contribution from YFP to total loss,\n($minYtoall-$maxYtoall)", titlefontsize = 10) #
    plot(h1,h2,h3,h4,size = (800,800), plot_title = "$(round.(rep_bestcf)), $(round.(sum(out[3])))")
end


#part 3: Real Test of candidate
#test each candidate
testonecf01(candidate1,0)
testonecf02(candidate1,0)
testonecf03(candidate1,0)
ncan = 1
testonecf01(cf_final_list[ncan],0)
testonecf02(cf_final_list[ncan],0)





