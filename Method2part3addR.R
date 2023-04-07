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
#I=2.429 R=17 L=0.003 
#or I=62.714 R=0.178 L=43.638  to find 2 branches when scanning di
#similar when scanning dl