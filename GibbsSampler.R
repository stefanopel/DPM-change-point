GibbsSampler = function(y,N,m1,m2,r=1,b=1,mu0=0,lam0=1,a.sig=1,b.sig=1,
                        a.alpha=0.5,b.alpha=0.0001,lambda=NULL,sigma2=NULL,S0=NULL,
                        alpha10=NULL,alpha20=NULL){

  # Inputs
  # y       = T vec of data
  # N       = scalar num of iterations
  # m1      = scalar max num of location change points
  # m2      = scalar max num of variance change points
  # r       = scalar initial balls for change point
  # b       = scalar initial balls for not change point 
  # mu0     = scalar location 1st hyperprior
  # lam0    = scalar location 2nd hyperprior
  # a.sig   = scalar variance 1st hyperprior
  # b.sig   = scalar variance 2nd hyperprior
  # a.alpha = scalar 1st alpha hyperprior
  # b.alpha = scalar 2nd alpha hyperprior
  # lambda  = scalar MCMC initial value for location
  # sigma2  = scalar MCMC initial value for variance
  # S0      = 2xT MCMC initial value for S 
  # alpha10 = scalar MCMC initial concentration for location DP
  # alpha20 = scalar MCMC initial concentration for variance DP

  # Outputs
  # alpha.out  = 2xN matrix of chains for alphas
  # S.out      = 2xNxT array if chains for S
  # theta1.out = N list with (m1+1) vectors with theta_1
  # theta2.out = N list with (m2+1) vectors with theta_2
  
  # preliminaries
  t0 = proc.time()
  library("mvtnorm")
  library("MCMCpack")
  library("MASS")
  T <- length(y) # sample size
  
  # MCMC initial values
  if(is.null(lambda))lambda <- rnorm(1,mu0,lam0)
  if(is.null(sigma2))sigma2 <- sqrt(rinvgamma(1,a.sig,b.sig))
  if(is.null(S0))S0 <- matrix(1,2,T)
  if(is.null(alpha10))alpha10 <- 1e3
  if(is.null(alpha20))alpha20 <- 1e3
  theta1 <- lambda
  theta2 <- sigma2
  S      <- S0
  alpha1 <- alpha10
  alpha2 <- alpha20
  
  # import functions for gibbs sampling
  source('SamplingS12iid.R')
  source('SamplingTheta1iid.R')
  source('SamplingTheta2iid.R')
  source('SamplingAlpha.R')
  
  #################################
  ####      Gibbs sampler      ####
  #################################
  S.out      <- array(1,dim=c(2,N,T))
  theta1.out <- list()
  theta2.out <- list()
  alpha.out  <- matrix(0,nr=2,nc=N) 
  
  for(i in 1:N){
    
    # complete missing values of theta1
    n.th <- length(theta1) 
    if(n.th<m1+1) theta1 <- c(theta1,rnorm(m1+1-n.th,mu0,lam0))
    
    # complete missing values of theta2
    n.th <- length(theta2)
    if(n.th<m2+1) theta2 <- c(theta2,sqrt(rinvgamma(m2+1-n.th,a.sig,b.sig)))    
    
    # sampling S
    S[1,] <- SamplingS1iid(S2=S[2,],theta1,theta2,m1,m2,T,y,r,b)
    S[2,] <- SamplingS2iid(S1=S[1,],theta1,theta2,m1,m2,T,y,r,b)

    # sampling theta1
    theta1[1:S[1,T]] <- SamplingTheta1iid(S,theta1,theta2,m1,m2,tau,mu0,lam0,y,T,alpha1,alpha2)
    
    # merge adjacent regimes with same theta1
    if(S[1,T]>1){
      tmp        <- diff(theta1[1:S[1,T]])
      to.merge   <- which(tmp==0)
      n.to.merge <- length(to.merge)
    } else n.to.merge <- 0
    if(n.to.merge>0){
      theta1 <- theta1[(1:S[1,T])[-to.merge]]
      n.th   <- S[1,T] - n.to.merge
      theta1 <- c(theta1,rnorm(m1+1-n.th,mu0,lam0))
      for(ii in n.to.merge:1)
        S[1,S[1,]==(to.merge[ii]+1)] <- to.merge[ii]
    }
    
    # sampling alpha1
    alpha1 <- SamplingAlpha(alpha1,theta1,a.alpha,b.alpha)     
    
    # sampling theta2 and alpha2
    theta2[1:S[2,T]] <- SamplingTheta2iid(S,theta1,theta2,m1,m2,tau,a.sig,b.sig,y,T,alpha1,alpha2)
    
    # merge adjacent regimes with same theta2
    if(S[2,T]>1){
      tmp        <- diff(theta2[1:S[2,T]])
      to.merge   <- which(tmp==0)
      n.to.merge <- length(to.merge)  
    } else n.to.merge <- 0
    if(n.to.merge>0){
      theta2 <- theta2[(1:S[2,T])[-to.merge]]
      n.th   <- S[2,T] - n.to.merge
      theta2 <- c(theta2,sqrt(rinvgamma(m2+1-n.th,a.sig,b.sig)))
      for(ii in n.to.merge:1)
        S[2,S[2,]==(to.merge[ii]+1)] <- to.merge[ii]
    }
    
    # sampling alpha2
    alpha2 <- SamplingAlpha(alpha2,theta2,a.alpha,b.alpha)   
    
    # store values
    alpha.out[1,i]  <- alpha1
    alpha.out[2,i]  <- alpha2
    S.out[,i,]      <- S
    theta1.out[[i]] <- theta1
    theta2.out[[i]] <- theta2
  }  
  t1 = proc.time()-t0
  
  return(list(alpha.out = alpha.out,
              S.out     = S.out,
              theta1.out = theta1.out,
              theta2.out = theta2.out,
              m1=m1,m2=m2,r=r,b=b,mu0=mu0,
              lam0=lam0,a.sig=a.sig,b.sig=b.sig,
              a.alpha=a.alpha,b.alpha=b.alpha,
              lambda=lambda,sigma2=sigma2,t1=t1,
              S0=S0,alpha10=alpha10,alpha20=alpha20))
}


