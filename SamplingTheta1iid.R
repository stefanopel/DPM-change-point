##############################################
### function for sampling theta1, iid case ###
##############################################
SamplingTheta1iid = function(S,theta1,theta2,m1,m2,tau,mu0,lam0,y,T,alpha1,alpha2){
  
  alpha = alpha1
  
  # fix change points from S
  tau   = which(diff(S[1,])==1)+1
  m     = length(tau)
  
  # compute current distinct values and clustering
  theta = theta1[1:(m+1)]
  theta.star = unique(theta)
  n.star = length(theta.star)
  conf = rep(0,m+1)
  for(i in 1:n.star) conf[theta==theta.star[i]] = i
  
  # if there are no regimes
  if(m==0){
    st = theta2[S[2,]]
    var.post  = 1/(1/lam0^2+sum(1/st^2))
    mean.post = (mu0/lam0^2+sum(y/st^2))*var.post
    theta     = rnorm(1,mean.post,sqrt(var.post))
    return(theta1=theta)}
  
  # sample clustering configuration
  for(i in 1:(m+1)){
    
    # cancel theta if not associated with other observations
    mi = sum(conf[-i] == conf[i])
    if(mi == 0) theta[i] = NA
    
    # identify observations
    if(i == 1){Tcluster = 1:(tau[1]-1)
    }  else if(i == m+1){Tcluster = tau[m]:T
    } else{Tcluster = tau[i-1]:(tau[i]-1)}
    Ti = length(Tcluster)
    
    st = theta2[S[2,Tcluster]]    
    Pi = sapply(1:n.star,function(j){
      mi = sum(conf[-i]==j)
      Pi = log(mi/(m+alpha))+sum(dnorm(y[Tcluster],theta.star[j],st,log=TRUE))})
    var.post  = (1/lam0^2+sum(1/st^2))^(-1)
    mean.post = (mu0/lam0^2+sum(y[Tcluster]/st^2))*var.post
    Pi = c(Pi,log(alpha/(m+alpha)*sqrt(var.post/2/pi)/lam0)-sum(log(st)) - (Ti+1)/2*log(2*pi) +
             0.5*(mean.post^2/var.post-(sum(y[Tcluster]^2/st^2)+mu0^2/lam0^2)))    
    Pi = exp(Pi-max(Pi))
    
    sampled = sample(n.star+1,1,prob=Pi)
    conf[i] = sampled
    if(sampled==n.star+1) {
      theta[i] = rnorm(1,mean.post,sqrt(var.post))
    } else{theta[i] = theta.star[sampled]}
    
    theta.star = unique(theta)
    n.star = length(theta.star)
    for(j in 1:n.star) conf[theta==theta.star[j]] = j
  }
  
  # sample unique values
  c.star = unique(conf)
  theta.star = sapply(1:length(c.star),function(i){
    
    # identify observations
    clusters = which(conf==c.star[i])
    Tclusters = c()
    for(j in 1:length(clusters)){
      if(clusters[j]==1){
        Tclusters = c(Tclusters,1:(tau[1]-1))
      } else if(clusters[j]==m+1){
        Tclusters = c(Tclusters,tau[m]:T)
      } else{Tclusters = c(Tclusters,tau[clusters[j]-1]:(tau[clusters[j]]-1))}}
    
    st        = theta2[S[2,Tclusters]]
    Ti        = length(Tclusters)
    var.post  = (1/lam0^2+sum(1/st^2))^(-1)
    mean.post = (mu0/lam0^2+sum(y[Tclusters]/st^2))*var.post
    th = rnorm(1,mean.post,sqrt(var.post))
    return(th)})
  
  # recover theta
  for(i in 1:length(c.star)) theta[conf==c.star[i]] = theta.star[i]
  
  return(theta1=theta) 
}
