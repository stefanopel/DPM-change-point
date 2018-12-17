#############################################################
#### function for sampling alpha (Escobar and West 1995) ####
#############################################################
SamplingAlpha = function(alpha,theta,a.alpha,b.alpha){
  
  # count values
  if(is.null(nrow(theta))){m = length(theta)
  } else {m = ncol(theta)}
  
  # count distinct values
  if(is.null(nrow(theta))){k = length(unique(theta))
  } else {k = ncol(t(unique(t(theta))))}  
  
  # auxiliary var for alpha
  eta.alpha = rbeta(1,alpha+1,m+1)         
  
  # sampling
  weight = (a.alpha + k - 1)/((m+1)*(b.alpha-log(eta.alpha)))
  weight = weight/(1+weight) 
  ind.alpha = rbinom(1,1,weight)
  alpha  = rgamma(1,a.alpha+k,b.alpha-log(eta.alpha))*ind.alpha + 
    rgamma(1,a.alpha+k-1,b.alpha-log(eta.alpha))*(1-ind.alpha)
  
  return(alpha=alpha)
}