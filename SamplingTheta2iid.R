##############################################
### function for sampling theta2, iid case ###
##############################################
SamplingTheta2iid = function(S,theta1,theta2,m1,m2,tau,a.sig,b.sig,y,T,alpha1,alpha2){
  
  alpha = alpha2
  
  # fix change points from S
  tau   = which(diff(S[2,])==1)+1
  m     = length(tau)
  
  # compute current distinct values and clustering
  theta = theta2[1:(m+1)]
  theta.star = unique(theta)
  n.star = length(theta.star)
  conf = rep(0,m+1)
  for(i in 1:n.star) conf[theta==theta.star[i]] = i
  
  # if there are no regimes
  if(m==0){
    mm = theta1[S[1,]]
    a.post = a.sig + T/2
    b.post = b.sig + sum((y-mm)^2)/2
    theta    = rinvgamma(1,a.post,b.post)
    return(theta2=sqrt(theta))}
  
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
    
    mm = theta1[S[1,Tcluster]]
    Pi = sapply(1:n.star,function(j){
      mi = sum(conf[-i]==j)
      Pi = log(mi/(m+alpha))+sum(dnorm(y[Tcluster],mm,theta.star[j],log=TRUE))})
    a.post = a.sig + Ti/2
    b.post = b.sig + 0.5*sum((y[Tcluster]-mm)^2)
    Pi = c(Pi,log(alpha/(m+alpha)*b.sig^a.sig)+lgamma(a.post)-lgamma(a.sig)-log((2*pi)^(Ti/2))-
             a.post*log(b.post))
    Pi = exp(Pi-max(Pi))
    
    sampled = sample(n.star+1,1,prob=Pi)
    conf[i] = sampled
    if(sampled==n.star+1) {
      theta[i] = sqrt(rinvgamma(1,a.post,b.post))
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
    
    mm     = theta1[S[1,Tclusters]]
    Ti     = length(Tclusters)
    a.post = a.sig + Ti/2
    b.post = b.sig + 0.5*sum((y[Tclusters]-mm)^2)
    th = sqrt(rinvgamma(1,a.post,b.post))
    return(th)})
  
  # recover theta
  for(i in 1:length(c.star)) theta[conf==c.star[i]] = theta.star[i]
  
  return(theta2=theta)
}