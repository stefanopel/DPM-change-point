"# DPM-change-point" 

    # generate some data
    y <- rnorm(100,0.5,0.3)   
    y <- c(y,rnorm(22,1,0.3))
    y <- c(y,rnorm(78,1,0.6))
    y <- c(y,rnorm(23,0.25,0.6))
    y <- c(y,rnorm(77,0.25,0.15))
    y <- c(y,rnorm(25,0.75,0.15))
    y <- c(y,rnorm(75,0.75,0.45))
    
    # run MCMC for 1000 iterations, 
    # 4 max change-points for mean                                
    # 4 max change-points for variance    
    out <- GibbsSampler(y,1000,4,4)
    
    # some output
    S.hat <- apply(out$S.out,c(1,3),mean)
    plot(S.hat[1,]) # indicator process estimate for mean regime
    for(i in 1:3) abline(v=i*100)
    
    plot(S.hat[2,]) # indicator process estimate for variance regime
    abline(v=122); abline(v=223); abline(v=325)

    plot(sapply(out$theta1.out,function(x)x[[1]]),type="l") # trace plot for mean in 1st regime
    plot(sapply(out$theta2.out,function(x)x[[2]]),type="l") # trace plot for sd in 2nd regime

    out$t1 # computation time
    
