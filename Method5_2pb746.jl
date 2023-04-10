#Method 5.2 - using the optimised link function to find Detstate and bMotile and bMatrix

#5.2.1. import and create basics
#model IRL
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
#result from 5.1
cf_final_list #the list of the optimised cf vectors for dILf(n) 
listofmodelID #the final model IDs we will use 

#ODE using either highR0 ic or LowR0 ic
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
#define function "dILf(n)" that produced δI and δL from a given nutrient and coefficient
function dIL_f(cf,n)
    #cf is the coefficient vector and n is the nutrient
    dL = cf[1]+ cf[2]*n +cf[3]*n^2 +cf[4]*n^3 +cf[5]*n^4
    dI = cf[6]*dL
    return (dI,dL)
end 


#5.2.2. Begin detIRL runs

#function that make dIdLTable for all time and distance points
function makedIdLTable(answer)
    dItable = Matrix{Float64}(undef, size(halfn))
    dLtable = similar(dItable)
    for i in 1:size(halfn,1) 
        for j in 1:size(halfn,2) 
            n = halfn[i,j]
            dIf = dIL_f(answer,n)[1] 
            dLf = dIL_f(answer,n)[2] 
            dItable[i,j] = dIf
            dLtable[i,j] = dLf
        end
    end
    #plotting if wanted
    # heatmap(dLtable, title = "δL from the model $(nn)\n$(round.(cf_final_list[nn],digits = 2))", titlefontsize = 10)
    # savefig("10NewOptimisation/From10Rounds/dLtable_plot_$(nn).svg")
    # heatmap(dItable, title = "δI from the model $(nn)\n$(round.(cf_final_list[nn],digits = 2))", titlefontsize = 10)
    # savefig("10NewOptimisation/From10Rounds/dItable_plot_$(nn).svg")
    return (dItable,dLtable)
end
#function that run for the time t = 1 across all distance points
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

#function that take two functions above and give out the whole model prediction for bMatrix(t,x) and bMotile(t,x)
function detRss_variedR0_allt(answer)
    Rss_variedR0_allt = Matrix{Float64}(undef,lastindex(distanceYRskip),lastindex(TimeInterval))
    Rss_variedR0_allt[:,1] = intiRssfromHR0(answer)
    dItable = makedIdLTable(answer)[1]
    dLtable = makedIdLTable(answer)[2]
    for i in 1:size(Rss_variedR0_allt,1)
        for j in 2:size(Rss_variedR0_allt,2)
            di = dItable[i,j]
            dl = dLtable[i,j]
            if Rss_variedR0_allt[i,j-1] > 10
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

dItable_array = Array{Matrix{Float64}}(undef, length(cf_final_list))
dLtable_array = similar(dItable_array)
detRss_variedR0_allt_array = Array{Matrix{Float64}}(undef, length(cf_final_list))
#Real running for all 10 models
for nn in listofmodelID #only run in the model the we are interested
    answer = cf_final_list[nn]
    #is the cf vector of the fitted dILf(n) learned using Optim.optimize
    function makedIdLTable(answer)
        dItable = Matrix{Float64}(undef, size(halfn))
        dLtable = similar(dItable)
        for i in 1:size(halfn,1) 
            for j in 1:size(halfn,2) 
                n = halfn[i,j]
                dIf = dIL_f(answer,n)[1] 
                dLf = dIL_f(answer,n)[2] 
                dItable[i,j] = dIf
                dLtable[i,j] = dLf
            end
        end
        #plotting if wanted
        # heatmap(dLtable, title = "δL from the model $(nn)\n$(round.(cf_final_list[nn],digits = 2))", titlefontsize = 10)
        # savefig("10NewOptimisation/From10Rounds/dLtable_plot_$(nn).svg")
        # heatmap(dItable, title = "δI from the model $(nn)\n$(round.(cf_final_list[nn],digits = 2))", titlefontsize = 10)
        # savefig("10NewOptimisation/From10Rounds/dItable_plot_$(nn).svg")
        return (dItable,dLtable)
    end

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
    
    function detRss_variedR0_allt(answer)
        Rss_variedR0_allt = Matrix{Float64}(undef,lastindex(distanceYRskip),lastindex(TimeInterval))
        dItable = makedIdLTable(answer)[1]
        dLtable = makedIdLTable(answer)[2]
        Rss_variedR0_allt[:,1] = intiRssfromHR0(answer)
        for i in 1:size(Rss_variedR0_allt,1)
            for j in 2:size(Rss_variedR0_allt,2)
                di = dItable[i,j]
                dl = dLtable[i,j]
                if Rss_variedR0_allt[i,j-1] > 10
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
    #call function makedIdLTable() to collect the tables
    dItable_array[nn] = makedIdLTable(answer)[1]
    dLtable_array[nn] = makedIdLTable(answer)[2]
    #call function detRss_variedR0_allt(), which calls the function intiRssfromHR0(), which calls the function  makedIdLTable to collect the matrix

    detRss_variedR0_allt_array[nn] = detRss_variedR0_allt(answer)
    #finish the DetIRL runs
    println("fin $(nn)")
end
#Checking each model for detIRL output
nn = 9#only consider nn in listofmodelID = [1,2,5,6,9]
p1 = heatmap(TimeInterval, distanceYRskip,dItable_array[nn], title = "δL from the model $(nn)\n$(round.(cf_final_list[nn],digits = 2))", xlab = "Time (hrs)", ylab = "Distance from centre (μm)",titlefontsize = 10)
p2 = heatmap(TimeInterval, distanceYRskip,dLtable_array[nn], title = "δI from the model $(nn)\n$(round.(cf_final_list[nn],digits = 2))", xlab = "Time (hrs)", ylab = "Distance from centre (μm)",titlefontsize = 10)
p3 = heatmap(TimeInterval, distanceYRskip,detRss_variedR0_allt_array[nn], xlab = "Time (hrs)", ylab = "Distance from centre (μm)", title = "DetRss, from the model $(nn)")
p4 = heatmap(TimeInterval, distanceYRskip,Rss_variedR0_allt.>10, xlab = "Time (hrs)", ylab = "Distance from centre (μm)", title = "DetState, from the model $(nn)\nyellow = Motile, black = Matrix", titlefontsize = 10)
plot(p1,p2,p3,p4, size = (800,800),tick_direction = :out)

#Begin stoIRL runs
#list of what needed
dI_span = 0:0.5:3
dL_span = 0:0.5:10 
bMotile_out_array = Array{Matrix{Float64}}(undef,length(cf_final_list))
bMatrix_out_array = similar(bMotile_out_array)
RationbMotile_out_array = similar(bMotile_out_array)
for nn in listofmodelID #model the we are interested
    answer = cf_final_list[nn]
    #start the stoIRL runs
    bMotile_out = Matrix{Float64}(undef,lastindex(distanceYRskip),lastindex(TimeInterval))
    bMatrix_out = similar(bMotile_out)
    glosstosum = similar(bMotile_out)
    finalPMot_fromHR0 = similar(bMotile_out)
    finalPMot_fromLR0 = similar(bMotile_out)
    finalPMot = similar(bMotile_out) 
    ny = 0.12
    nr = 0.08
    dI_span = 0:0.5:3
    dL_span = 0:0.5:10 
    bMotile = 0
    bMatrix = 0
    δI_f = 0
    δL_f = 0
    mdi_ix = 0
    mdl_ix = 0
    for i in 1:lastindex(distanceYRskip)
        for j in 1:lastindex(TimeInterval)
            #define experimental data
            YFP = spreadYFPm_g[i,j]
            RFP = spreadRFPm_g[i,j]
            #define model prediction - to be used later
            n = halfn[i,j]
            b = halfb[i,j]
            #define δI and δL from the link function, 
            #with the rounding to 0.5, given a predicted n
            δI_f = (ceil(dIL_f(answer,n)[1]) + floor(dIL_f(answer,n)[1]))/2
            δL_f = (ceil(dIL_f(answer,n)[2]) + floor(dIL_f(answer,n)[2]))/2
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
            bMotile_out[i,j] = bMotile
            bMatrix_out[i,j] = bMatrix
        end
        #println("fin - $(distanceYRskip[i])")
    end
    ratio_bmotile = bMotile_out./(bMotile_out.+bMatrix_out)
    bMotile_out_array[nn] = bMotile_out
    bMatrix_out_array[nn] = bMatrix_out
    RationbMotile_out_array[nn] = ratio_bmotile
    #plotting and exporting if wanted
    #writedlm("10NewOptimisation/From10Rounds/bMotile_out_$(nn)",bMotile_out)
    #writedlm("10NewOptimisation/From10Rounds/bMatrix_out_$(nn)",bMatrix_out)
    # p1 = heatmap(TimeInterval, distanceYRskip,bMotile_out,xlab = "Time (hrs)", ylab = "Distance from centre (μm)", title = "Density of motile bacteria, from the model $(nn)", titlefontsize = 10, tick_direction = :out, color = :rainbow, clim = (0,5))
    # p2 = heatmap(TimeInterval, distanceYRskip,bMatrix_out,xlab = "Time (hrs)", ylab = "Distance from centre (μm)", title = "Density of matrix bacteria, from the model $(nn)", titlefontsize = 10, tick_direction = :out, color = :rainbow,clim = (0,5))
    # plot(p1,p2, layout = (1,2), size = (800,500))
    # savefig("10NewOptimisation/From10Rounds/0DensitySepOut_$(nn).svg")
    #heatmap(TimeInterval, distanceYRskip,ratio_bmotile,xlab = "Time (hrs)", ylab = "Distance from centre (μm)", title = "The ratio of Motile bacteria to total, from the model $(nn)", titlefontsize = 10, tick_direction = :out, color = :rainbow)
    #savefig("10NewOptimisation/From10Rounds/NewDensityRatio_$(nn).svg")
end

#checking stoIRL result from each model 
nn = 9#only consider nn in listofmodelID = [1,2,5,6,9] 
p5 = heatmap(TimeInterval, distanceYRskip,bMotile_out_array[nn],xlab = "Time (hrs)", ylab = "Distance from centre (μm)", title = "Density of motile bacteria, from the model $(nn)", titlefontsize = 10, tick_direction = :out, color = :rainbow, clim = (0,5))
p6 = heatmap(TimeInterval, distanceYRskip,bMatrix_out_array[nn],xlab = "Time (hrs)", ylab = "Distance from centre (μm)", title = "Density of matrix bacteria, from the model $(nn)", titlefontsize = 10, tick_direction = :out, color = :rainbow,clim = (0,5))
p7 = heatmap(TimeInterval, distanceYRskip,RationbMotile_out_array[nn],xlab = "Time (hrs)", ylab = "Distance from centre (μm)", title = "The ratio of Motile bacteria to total,\nfrom the model $(nn)", titlefontsize = 10, tick_direction = :out, color = :rainbow)
plot(p5,p6,p7,tick_direction = :out,size = (800,800) )


