# Load the library deSolve to use differential equations
# install.packages("deSolve")
library(deSolve)


# Make the function for the confidence limits
"bootstrap"<- function(x,nboot,theta,...,func=NULL){
  
  call<-match.call()
  
  n<-length(x)
  
  bootsam<-matrix(sample(x,size=n*nboot,replace=T),nrow=nboot)
  
  thetastar<-apply(bootsam,1,theta,...)
  
  func.thetastar<-NULL; jack.boot.val<-NULL; jack.boot.se<-NULL;
  if(!is.null(func)){
    match1<-function(bootx,x){duplicated(c(bootx,x))[(length(x)+1):(2*length(x))]}
    matchs<-t(apply(bootsam,1,match1,x))
    func.thetastar<-func(thetastar)
    jack.boot<-function(inout,thetastar,func){func(thetastar[!inout])}
    jack.boot.val<-apply(matchs,2,jack.boot,thetastar,func)
    
    if (sum(is.na(jack.boot.val)>0)) {cat("At least one jackknife influence value for func(theta) is undefined",fill=T)
      cat("Increase nboot and try again",fill=T)
      return()}
    
    if (sum(is.na(jack.boot.val))==0)
    {jack.boot.se<-sqrt( ((n-1)/n)*sum((jack.boot.val-mean(jack.boot.val))^2))
    
    }}
  
  return(list(thetastar=thetastar,func.thetastar=func.thetastar,jack.boot.val=jack.boot.val,jack.boot.se=jack.boot.se,call=call))
  
}

boot <- data.frame



# Make the function for the Maximum Likelihood approach
ML.function <- function(n.sim, time.sim, model, p, tau, tf, n.boot){
  
  # Create an empty storage for the estimates of p 
  estimates <- matrix(0, 1, n.sim)
  
  # Create an empty storage for the 1 and 0 values
  coverage <- matrix(0, 1, n.sim)
  
  # Put the estimation process in a for loop to run it nsim times 
  for(x in 1:n.sim){
    
    # simulate data
    out <- ode(y = c(L = 0), times = time.sim, func = model, parms = c(p,tau))
    
    # extract just the second column from this simulated data set
    data.simulated <- out[,2]  
    
    # add noise
    noise <- rnorm(10,0,p/10) 
    data.sim.noise <- data.simulated + noise
    
    # make the cost function with the Maximum Likelihood approach
    cost.ML <- function(psi){
      
      out <- ode(y = c (L = 0), times = time.sim, func = model, parms = c(p = exp(psi[1]), tau))
      data_pred <- out[ ,2]
      
      if(tf == "none") out1<-ifelse(data_pred<0,0,data_pred)
      if(tf == "arcsinesq") out1<-ifelse(data_pred<0,0,asin(sqrt(data_pred)))
      if(tf == "sqroot") out1<-ifelse(data_pred<0,0,sqrt(data_pred))
      if(tf == "log") out1<-ifelse(data_pred<=0,-10**(-10),log(data_pred))
      if(tf == "qdroot") out1<-ifelse(data_pred<0,0,(data_pred)^(0.25))
      if(tf == "arcsineqd") out1<-ifelse(data_pred<0,0,asin((data_pred)^(0.25)))
      
      out_cost <- (data.fit - out1) 
      
      nobs<-length(data.fit)    
      
      SSR <- sum(out_cost^2)
      Lik.l <- ((2*pi*(exp(psi[2]))^2)^(-nobs/2))*exp(-SSR/(2*(exp(psi[2]))^2))
      minlogl <- -log(Lik.l)
      return(minlogl)
      
    }
    
    # transform the data 
    if(tf == "none") data.fit <- ifelse(data.sim.noise<0,0,data.sim.noise)
    if(tf == "arcsinesq") data.fit <- ifelse(data.sim.noise<0,0,asin(sqrt(data.sim.noise)))
    if(tf == "sqroot") data.fit <- ifelse(data.sim.noise<0,0,sqrt(data.sim.noise))
    if(tf == "log") data.fit <- ifelse(data.sim.noise<=0, -10**(-10),log(data.sim.noise))
    if(tf == "qdroot") data.fit <- ifelse(data.sim.noise<0,0,(data.sim.noise)^(0.25))
    if(tf == "arcsineqd") data.fit <- ifelse(data.sim.noise<0,0,asin((data.sim.noise)^(0.25)))
    
    # optimization
    estim<-nlminb(c(-6,-2),cost.ML,upper=c(0,Inf),lower=c(-15,-Inf))
    
    # get the parameter estimate by taking the exponent 
    prm.ML<-exp(estim$par) 
    
    # store each parameter estimate in the empty matrix we created before 
    estimates[ ,x] <- prm.ML[1] # so you keep only the estimates of p 
    
    ### function  theta.res will compute all the boostrap estimates
    theta.res<-function(residual){
      res1<-residual  
      tim<- time.sim
      tau<-63
      
      # define the new function to optimize with the adding residuals to the data
      cost.res<-function(psi){
        
        out<-ode(y=c(L=0), times=tim, func = model1, parms =c(p=exp(psi[1]),tau=tau))
        pred<-out[,2]
        
        if(tf == "none") out1<-ifelse(pred<0,0,pred)
        if(tf == "arcsinesq") out1<-ifelse(pred<0,0,asin(sqrt(pred)))
        if(tf == "sqroot") out1<-ifelse(pred<0,0,sqrt(pred))
        if(tf == "log") out1<-ifelse(pred<=0,-10^(-10),log(pred))
        if(tf == "qdroot") out1<-ifelse(pred<0,0,(pred)^(0.25))
        if(tf == "arcsineqd") out1<-ifelse(pred<0,0,asin((pred)^(0.25)))
        
        if(tf == "none") data.fit.res<-ifelse(data.sim.noise<0,0,data.sim.noise+res1)
        if(tf == "arcsinesq") data.fit.res<-ifelse(data.sim.noise<0,0,asin(sqrt(data.sim.noise+res1)))
        if(tf == "sqroot") data.fit.res<-ifelse(data.sim.noise<0,0,sqrt(data.sim.noise+res1))
        if(tf == "log") data.fit.res <- ifelse(data.sim.noise<=0, -10**(-10),log(data.sim.noise+res1))
        if(tf == "qdroot") data.fit.res <- ifelse(data.sim.noise<0,0,(data.sim.noise+res1)^(0.25))
        if(tf == "arcsineqd") data.fit.res <- ifelse(data.sim.noise<0,0,asin((data.sim.noise+res1)^(0.25)))
        
        out_cost <- (data.fit.res - out1) 
        
        nobs<-length(data.fit.res)    
        
        SSR <- sum(out_cost^2)
        Lik.l <- ((2*pi*(exp(psi[2]))^2)^(-nobs/2))*exp(-SSR/(2*(exp(psi[2]))^2))
        minlogl <- -log(Lik.l)
        return(minlogl)
      }
      
      estim<-nlminb(as.numeric(log(prm.ML)),cost.res)
      return(estim$par)
    }
    
    
    # BOOTSTRAP
    
    out.res<-ode(y=c(L=0), times=time.sim, func = model1, parms =c(p=prm.ML,tau=tau), rtol = 1e-12,method = "lsode")

    res<- data.sim.noise-out.res[,2] # compute the residuals
    
    boot<-bootstrap(res,n.boot,theta.res)   # here we use n.boot bootstraps with the residual res and the function theta.res to estimate all boostrap step
    resu<-exp(boot$thetastar)  # this has the estimates for each bootstrap step
    
    # Get the confidence limits
    
    CI <- as.vector(quantile(resu, probs=c(0.025, 0.975)))
    ifelse(p > CI[1] & p < CI[2], coverage[1, x] <- 1, coverage[1, x] <- 0)
    
  }

  
  # PERFORMANCE MEASURES

  # Bias
  bias <- p - mean(estimates[1, ]) # mean of all parameter estimates stored in the estimates matrix 
  print(bias)
  
  # MSE
  MSE <- sum(estimates[1,]-p)^2*1/n.sim
  print(MSE)
  
  # RSE
  SE<- sd(estimates[1,])/sqrt(length(estimates[1,]))
  RSE<-SE/prm.ML[1]
  print(RSE)
  
  # Coverage
  cov <- mean(coverage[1, ])
  print(cov)
}



# Specify the model we use to generate the data
model1 <- function(t,state,pars){
  with(as.list(c(state,pars)),{
    dL <- ifelse(t<=tau,p -  p* L, -p*L)
    return(list(c(dL)))
  })
}

# Specify the timepoints 
time.sim <- c(0,7, 14, 28, 42, 63, 70, 84, 105, 126)

# Define tau
tau <- 63

# ALL SIMULATIONS WITH P = 0.01

# Define p
p <- 0.01

# No transformation
set.seed(987)
ML.function(n.sim = 100, time.sim = time.sim, model = model1, p = p, tau = tau, tf = "none", n.boot = 200)

# Arcsine square root transformation
set.seed(987)
ML.function(n.sim = 100, time.sim = time.sim, model = model1, p = p, tau = tau, tf = "arcsinesq", n.boot = 200)

# Square root transformation
set.seed(987)
ML.function(n.sim = 100, time.sim = time.sim, model = model1, p = p, tau = tau, tf = "sqroot", n.boot = 200)

# Log transformation
set.seed(987)
ML.function(n.sim = 100, time.sim = time.sim, model = model1, p = p, tau = tau, tf = "log", n.boot = 200)

# Quadratic root transformation
set.seed(987)
ML.function(n.sim = 100, time.sim = time.sim, model = model1, p = p, tau = tau, tf = "qdroot", n.boot = 200)

# Arcsine quadratic root transformation
set.seed(987)
ML.function(n.sim = 100, time.sim = time.sim, model = model1, p = p, tau = tau, tf = "arcsineqd", n.boot = 200)


# ALL SIMULATIONS WITH P = 0.0005

# Define p again
p <- 0.0005

# No transformation
set.seed(987)
ML.function(n.sim = 100, time.sim = time.sim, model = model1, p = p, tau = tau, tf = "none", n.boot = 200)

# Arcsine square root transformation
set.seed(987)
ML.function(n.sim = 100, time.sim = time.sim, model = model1, p = p, tau = tau, tf = "arcsinesq", n.boot = 200)

# Square root transformation
set.seed(987)
ML.function(n.sim = 100, time.sim = time.sim, model = model1, p = p, tau = tau, tf = "sqroot", n.boot = 200)

# Log transformation
set.seed(987)
ML.function(n.sim = 100, time.sim = time.sim, model = model1, p = p, tau = tau, tf = "log", n.boot = 200)

# Quadratic root transformation
set.seed(987)
ML.function(n.sim = 100, time.sim = time.sim, model = model1, p = p, tau = tau, tf = "qdroot", n.boot = 200)

# Arcsine quadratic root transformation
set.seed(987)
ML.function(n.sim = 100, time.sim = time.sim, model = model1, p = p, tau = tau, tf = "arcsineqd", n.boot = 200)

