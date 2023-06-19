#R code for analysis of IRL model (instead of using Julia's Bifurcationkit.jl)
library(deBif)
#Create
state <- c(I = 100., R = 300. , L = 100.)
parms <- c(a=85., b=100. , g=125. , di =2. , dl=10. )
IRLsys <-function(t,state,parms){
  with(as.list(c(state,parms)), {
    dI = a - I - di*I*R 
    dR = b - R - di*I*R - dl*L*R
    dL = g/(1+R^2) - L - dl*L*R
    
    return(list(c(dI, dR, dL)))
  })
}
#implement
bifurcation(IRLsys, state, parms)
phaseplan(IRLsys, state, parms)
