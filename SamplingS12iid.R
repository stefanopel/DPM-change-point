###############################################
### function for sampling S[1,] given S[2,] ###
###############################################
SamplingS1iid = function(S2=S[2,],theta1,theta2,m1,m2,T,y,r,b){
  
  # forward step
  p = matrix(0,nr=m1+1,nc=T)
  p[1,1] = 1;
  
  ff = sapply(2:T,function(t)
    dnorm(y[t],theta1,theta2[S2[t]],log=TRUE))
  
  for(t in 2:T){
    for(j in 1:m1) p[j+1,t] = log(b/(r+b)*p[j+1,t-1]+r/(r+b)*p[j,t-1])
    p[1,t] = log(b/(r+b)*p[1,t-1])
    p[,t] = p[,t] + ff[,t-1]
    p[,t] = exp(p[,t]-max(p[,t]))
  }
  
  S1 = rep(1,T)
  
  # backward
  S1[T] = sample(m1+1,1,prob=p[,T])
  if(S1[T]==1) return(S1=S1)
  
  for(t in (T-1):2){
    active = c(S1[t+1]-1,S1[t+1]) 
    p[active[2],t] = p[active[2],t]*b/(r+b)
    p[active[1],t] = p[active[1],t]*r/(r+b)
    S1[t]  = sample(active,1,prob=p[active,t])
    if(S1[t]==1) return(S1=S1)}
  
  return(S1=S1)
}


###############################################
### function for sampling S[2,] given S[1,] ###
###############################################
SamplingS2iid = function(S1=S[1,],theta1,theta2,m1,m2,T,y,r,b){
  
  # forward step
  p = matrix(0,nr=m2+1,nc=T)
  p[1,1] = 1;
  
  ff = sapply(2:T,function(t)
    dnorm(y[t],theta1[S1[t]],theta2,log=TRUE))
  
  for(t in 2:T){
    for(j in 1:m2) p[j+1,t] = log(b/(r+b)*p[j+1,t-1]+r/(r+b)*p[j,t-1])
    p[1,t] = log(b/(r+b)*p[1,t-1])
    p[,t] = p[,t] + ff[,t-1]
    p[,t] = exp(p[,t]-max(p[,t]))
  }
  
  S2 = rep(1,T)
  
  # backward
  S2[T] = sample(m2+1,1,prob=p[,T])
  if(S2[T]==1) return(S2=S2)
  
  for(t in (T-1):2){
    active = c(S2[t+1]-1,S2[t+1]) 
    p[active[2],t] = p[active[2],t]*b/(r+b)
    p[active[1],t] = p[active[1],t]*r/(r+b)
    S2[t]  = sample(active,1,prob=p[active,t])
    if(S2[t]==1) return(S2=S2)}
  
  return(S2=S2)
}