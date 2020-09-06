require(pomp)

#################################
##########*** EAKF ***###########
#################################
eakf.internal <- function(object, params, Np, mifiter, rw.sd, cooling.fn,
                          tol = 0, max.fail = Inf, verbose,tolerance,
                          C,R,check.fail,
                          .indices = integer(0),
                          .gnsi = TRUE){
  #####* General initialization *#####
  
  ## tol is the tolerance used to test the filtering 
  ## failure which defined in Particle Filter
  tol <- as.numeric(tol)
  gnsi <- as.logical(.gnsi)
  verbose <- as.logical(verbose)
  mifiter <- as.integer(mifiter)
  Np <- as.integer(Np)
  
  if (length(tol) != 1 || !is.finite(tol) || tol < 0)
    pStop_(sQuote("tol")," should be a small nonnegative number.")
  
  do_ta <- length(.indices)>0L
  if (do_ta && length(.indices)!=Np[1L])
    pStop_(sQuote(".indices")," has improper length.")
  
  times <- time(object,t0=TRUE)
  ntimes <- length(times)-1
  
  ## Check the numer of filtering failure
  nfail <- 0
  
  #####* Kalman initialization *#####
  loglik <- numeric(ntimes)
  y <- obs(object) ## observation time series data
  nobs <- nrow(y) ## dim of observation 
  Y <- array(dim=c(nobs,Np),dimnames=list(variable=rownames(y),rep=NULL)) ## ensemble of measurements
  X <- rinit(object,params=coef(object),nsim=Np) ## ensemble of hidden states
  nvar.names <- rownames(rinit(object,params=coef(object),nsim=Np)) ## name of hidden state
  nvar <- nrow(X) ## dim of hidden state 
  
  ## names of parameters which do random walk
  npar.names <- names(cooling.fn(1,mifiter)$alpha*rw.sd[,1]) 
  
  ## dim of parameters which do random walk
  npar <- length(npar.names)
  
  ### Create extended observation matrix C.e
  C.p <- matrix(rep(0,npar), nobs, npar, byrow = TRUE)
  colnames(C.p) <- npar.names
  C.e <- cbind(C,C.p)
  
  #####* Start of the main loop *#####
  for (nt in seq_len(ntimes)) {
    
    ## decide parameter variance at nt
    pmag <- cooling.fn(nt,mifiter)$alpha*rw.sd[,nt]
    
    ## Forecast 1.0 : Add noise to parameters
    for (par.name in npar.names) {
      params[par.name,] <- rnorm(n = dim(params)[2], 
                                 mean = params[par.name,],
                                 sd = pmag[par.name])
    }
    
    ## Forecast 1.1 : Change parameters back to the original scale
    tparams <- pomp::partrans(object,params,dir="fromEst",.gnsi=gnsi)
    
    
    ## Forecast 1.2  Get initial states (if nt = 1) and add noise
    if (nt == 1L) {
      x <- pomp::rinit(object,params=tparams)
    }
    
    ## Forecast 2.1 Advance the state variables according to the process model
    X[,] <- pomp::rprocess(
      object,
      x0=x,
      t0=times[nt],
      times=times[nt+1],
      params=tparams,
      .gnsi=gnsi
    )
    
    ######* Check if failure happens, not necessary in algorithm*#####
    if(check.fail){
      weights <- pomp::dmeasure(
        object,
        y=object@data[,nt,drop=FALSE],
        x=array(X,dim=c(dim(X),1),dimnames=list(variable=rownames(X))),
        times=times[nt+1],
        params=tparams,
        log=FALSE,
        .gnsi=gnsi
      )
      
      ## Apply tolerance to decide which ensemble fail
      index.missing <- weights <= tol
      weights[index.missing] <- 0
      weights[is.na(weights)] <- 0
      loglik[nt] <- ifelse(sum(weights) == 0, 
                           #log(1e-17),
                           #log(.Machine$double.xmin), 
                           -Inf,
                           log(mean(weights)))
      
      
      if(  sum(weights) == 0){
        nfail <- nfail + 1
      }
    }

    #####* Compute weighted mean at last timestep if no faiure, otherwise take unweighted mean *#####
    if(nt == ntimes){
      weights <- pomp::dmeasure(
        object,
        y=object@data[,nt,drop=FALSE],
        x=array(X,dim=c(dim(X),1),dimnames=list(variable=rownames(X))),
        times=times[nt+1],
        params=tparams,
        log=FALSE,
        .gnsi=gnsi
      )
      
      ##*Apply tolerance
      index.missing <- weights <= tol
      weights[index.missing] <- 0
      
      ##*Normalize
      weights <- weights/sum(weights)
      
      ##*Remove NA, aovid all 0 since all 0 will produce NaN
      weights[is.na(weights)] <- 0
      
      ## Take weighted sum if no failure
      if(sum(weights) != 0){
        coef(object,transform=TRUE) <- apply(X = params,
                                             MARGIN = 1,
                                             FUN = weighted.mean,
                                             w = weights)
      }
      else{
      ## Otherwise take unweighted sum
        coef(object,transform=TRUE) <- apply(X = params,
                                             MARGIN = 1,
                                             FUN = mean)
      }
    }
    
    ## Forecast 3.0 Get the extended hidden state (hidden state + parameters which can change)
    x.e <- rbind(X,params[npar.names,,drop=F])
    
    #####* Updare steps of EAKF *#####
    
    ## Compute the prediction mean of the transformed parameters and hidden states 
    pm <- rowMeans(x.e)  ## prediction mean
    x.e <- x.e-pm
    
    ## Approximate the measurement noise by function R 
    ## and predicted mean parameters and hidden states 
    #R.num <- R(rowMeans(rbind(X,tparams[npar.names,])))
    R.num <- R(rowMeans(rbind(X,tparams[,])))
    if(any(diag(R.num) < 0) || any(is.na(R.num)) ){
      stop('Measurement error is uncomputable')
    }
    # sqrtR <- ifelse( any(R.num == 0),
    #                  0,t(chol(R.num)))
    
    
    if(length(R.num) > 1){
      if(any(diag(R.num) == 0)){
        index <- diag(R.num) == 0
        diag(R.num)[index] <- tolerance
        sqrtR <- t(chol(R.num))
        diag(sqrtR)[index] <- 0
      }else{
        sqrtR <- t(chol(R.num))
      }

    }else{
      if(R.num == 0){
        R.num <- tolerance
        sqrtR <- 0
      }else{
        sqrtR <- t(chol(R.num))
      }
    }
    
    
    
    
    

    
    pv <- tcrossprod(x.e)/(Np-1)      ## prediction variance
    svdV <- svd(pv,nv = 0)
    
    forecast <- C.e %*% pm        ## forecast (observables)
    resid <- y[,nt]-forecast     ## forecast error
    
    ww <- t(C.e%*%svdV$u)
    w <- crossprod(ww,svdV$d*ww)+R.num  ## forecast variance
    svdW <- svd(w,nv=0)
    
    ## compute the likelihood using the extended hidden state using matrix
    ##loglik[nt] <- sum(dnorm(x=crossprod(svdW$u,resid),mean=0,
                            ##sd=sqrt(svdW$d),log=TRUE))
    
    u <- sqrt(svdV$d)*ww
    u <- tcrossprod(u%*%sqrtR,u)
    svdU <- svd(u,nv=0)

    ## b is adjustment
    #### Since sqrt(svdV$d) may have 0 entry, divid by it will lead to Inf
    #### To avoid that, we give sqrt(svdV$d) a small value = tolerance
    if(any(svdV$d == 0) ){
      warning("SVD has numerical issue, adjustment is approximated in computation")
      b <- svdV$u%*%(sqrt(svdV$d)*svdU$u)%*%
        (1/sqrt(1+svdU$d)/(sqrt(svdV$d)+tolerance)*t(svdV$u))
    }else{
      b <- svdV$u%*%(sqrt(svdV$d)*svdU$u)%*%
        (1/sqrt(1+svdU$d)/sqrt(svdV$d)*t(svdV$u))
    }

    if( any(is.na(b)) ){
      stop('Adjustment is non-computable due to numberical reason, try to increase the number of ensembles to avoid that')
    }
    
    K <- tcrossprod(b%*%pv,b)%*%crossprod(C.e,sqrtR)   # Kalman gain
    
    if( any(is.na(K)) ){
      stop('Kalman Gain is non-computable due to numberical reason, try to increase the number of ensembles to avoid that')
    }
    
    fm <- pm+K%*%resid         ## filter mean
    
    x.e[,] <- b%*%x.e+as.vector(fm)  ## filter ensembles 
  
    
    ## break the extended hidden stateto changable parameters and hidden state
    x <- x.e[nvar.names,]
    params[npar.names,] <- x.e[npar.names,]
  }
  
  return(new(
    "pfilterd_pomp",
    as(object,"pomp"),
    paramMatrix=params,
    cond.loglik=loglik,
    indices=.indices,
    Np=Np,
    nfail = as.integer(nfail),
    tol=tol,
    loglik=sum(loglik)
  ))}



mif.eakf <- function(object, Nmif, rw.sd,
                     cooling.fraction.50,
                     C,R,
                     Np, max.fail = Inf,
                     ..., verbose = FALSE,check.fail = FALSE,
                     cooling.type = "geometric",tol = 0,
                     .ndone = 0L, .indices = integer(0), .paramMatrix = NULL,
                     .gnsi = TRUE,tolerance = 1e-200){
  
  ### Default initialization in POMP package
  verbose <- as.logical(verbose)
  check.fail <- as.logical(check.fail)
  object <- pomp(object,verbose=verbose)
  
  gnsi <- as.logical(.gnsi)
  
  if (length(Nmif) != 1 || !is.numeric(Nmif) || !is.finite(Nmif) || Nmif < 1)
    pStop_(sQuote("Nmif")," must be a positive integer.")
  Nmif <- as.integer(Nmif)
  
  if (is.null(.paramMatrix)) {
    start <- coef(object)
  } else {  
    start <- apply(.paramMatrix,1L,mean)
  }
  
  ntimes <- length(pomp::time(object))
  
  if (is.null(Np)) {
    pStop_(sQuote("Np")," must be specified.")
  } else if (is.function(Np)) {
    Np <- tryCatch(
      vapply(seq_len(ntimes),Np,numeric(1)),
      error = function (e) {
        pStop_("if ",sQuote("Np"),
               " is a function, it must return a single positive integer.")
      }
    )
  } else if (!is.numeric(Np)) {
    pStop_(sQuote("Np"),
           " must be a number, a vector of numbers, or a function.")
  }
  
  
  if (!all(is.finite(Np)) || any(Np <= 0))
    pStop_(sQuote("Np")," must be a positive integer.")
  
  Np <- as.integer(Np)
  
  if (missing(rw.sd))
    pStop_(sQuote("rw.sd")," must be specified!")
  
  ###### Create list of perturbed noise, use perturbn.kernel.sd function ######
  rw.sd <- perturbn.kernel.sd(rw.sd,time=time(object),paramnames=names(start))
  
  if (missing(cooling.fraction.50))
    pStop_(sQuote("cooling.fraction.50")," is a required argument.")
  if (length(cooling.fraction.50) != 1 || !is.numeric(cooling.fraction.50) ||
      !is.finite(cooling.fraction.50) || cooling.fraction.50 <= 0 ||
      cooling.fraction.50 > 1)
    pStop_(sQuote("cooling.fraction.50")," must be in (0,1].")
  cooling.fraction.50 <- as.numeric(cooling.fraction.50)
  
  ###### Choose cooling method and time, use mif2.cooling function ######
  cooling.fn <- mif2.cooling(
    type=cooling.type,
    fraction=cooling.fraction.50,
    ntimes=length(time(object))
  )
  
  if (is.null(.paramMatrix)) {
    paramMatrix <- array(data=start,dim=c(length(start),Np[1L]),
                         dimnames=list(variable=names(start),rep=NULL))
  } else {
    paramMatrix <- .paramMatrix
  }
  
  traces <- array(dim=c(Nmif+1,length(start)+1+
                          1),
                  dimnames=list(iteration=seq.int(.ndone,.ndone+Nmif),
                                variable=c("loglik","nfail",names(start))))
  traces[1L,] <- c(NA,NA,start)
  
  pomp::pompLoad(object,verbose=verbose)
  on.exit(pompUnload(object,verbose=verbose))
  
  ## change the parameter scale, change them back in the loop of enkf
  paramMatrix <- partrans(object,paramMatrix,dir="toEst",
                          .gnsi=gnsi)
  
  ###### iterate the filtering of EAKF #######
  
  ## iterate the filtering
  for (n in seq_len(Nmif)) {
    pfp <- eakf.internal(
      object=object,
      params=paramMatrix,
      Np= Np,
      R = R,
      C = C,
      mifiter=.ndone+n,
      cooling.fn=cooling.fn,
      rw.sd=rw.sd,
      tol=tol,
      check.fail = check.fail,
      max.fail=max.fail,
      verbose=verbose,
      .indices=.indices,
      .gnsi=gnsi,
      tolerance = tolerance
    )
    
    gnsi <- FALSE
    
    paramMatrix <- pfp@paramMatrix
    traces[n+1,-c(1L,2L)] <- coef(pfp)
    traces[n,1L] <- pfp@loglik
    traces[n,2L] <- pfp@nfail
    .indices <- pfp@indices
    
    if (verbose) cat("mif2 EnKF iteration",n,"of",Nmif,"completed\n")
    
  }
  pfp@paramMatrix <- partrans(object,paramMatrix,dir="fromEst",.gnsi=gnsi)
  new(
    "mif2d_pomp",
    pfp,
    Nmif=Nmif,
    rw.sd=rw.sd,
    cooling.type=cooling.type,
    cooling.fraction.50=cooling.fraction.50,
    traces=traces
  )}


#################################
##########*** EnKF ***###########
#################################
enkf.internal <- function(object, params, Np, mifiter, rw.sd, cooling.fn,
                            tol = 0, max.fail = Inf, verbose,
                            h,R,check.fail,
                            .indices = integer(0),
                            .gnsi = TRUE,tolerance){
  #####* General initialization *#####
  tol <- as.numeric(tol)
  gnsi <- as.logical(.gnsi)
  verbose <- as.logical(verbose)
  mifiter <- as.integer(mifiter)
  Np <- as.integer(Np)
  
  if (length(tol) != 1 || !is.finite(tol) || tol < 0)
    pStop_(sQuote("tol")," should be a small nonnegative number.")
  
  do_ta <- length(.indices)>0L
  if (do_ta && length(.indices)!=Np[1L])
    pStop_(sQuote(".indices")," has improper length.")
  
  times <- time(object,t0=TRUE)
  ntimes <- length(times)-1
  
  ## Check the numer of filtering failure
  nfail <- 0
  
  #####* Kalman initialization *#####
  loglik <- numeric(ntimes)
  y <- obs(object) ## observation time series data
  nobs <- nrow(y) ## dim of observation 
  Y <- array(dim=c(nobs,Np),dimnames=list(variable=rownames(y),rep=NULL)) ## ensemble of measurements
  X <- rinit(object,params=coef(object),nsim=Np) ## ensemble of hidden states
  nvar.names <- rownames(rinit(object,params=coef(object),nsim=Np)) ## name of hidden state
  nvar <- nrow(X) ## dim of hidden state 
  
  ## name of parameters which do random walk
  npar.names <- names(cooling.fn(1,mifiter)$alpha*rw.sd[,1]) 
  
  #####* Start of the main loop *#####
  for (nt in seq_len(ntimes)) {
    
    ## decide parameter variance at nt
    pmag <- cooling.fn(nt,mifiter)$alpha*rw.sd[,nt]
    
    ## Forecast 1.0 : Add noise to parameters
    for (par.name in npar.names) {
      params[par.name,] <- rnorm(n = dim(params)[2], 
                                 mean = params[par.name,],
                                 sd = pmag[par.name])
    }
    
    ## Forecast 1.1 : Change parameters back to original scale 
    tparams <- pomp::partrans(object,params,dir="fromEst",.gnsi=gnsi)
    
    
    ## Forecast 1.2 Get initial states (if nt = 1) and add noise
    if (nt == 1L) {
      x <- pomp::rinit(object,params=tparams)
    }
    
    ## Forecast 2.1 Advance the hidden states according to the process model
    X[,] <- pomp::rprocess(
      object,
      x0=x,
      t0=times[nt],
      times=times[nt+1],
      params=tparams,
      .gnsi=gnsi
    )
    

    ######* Check if failure happens, not necessary in algorithm*#####
    if(check.fail){
      weights <- pomp::dmeasure(
        object,
        y=object@data[,nt,drop=FALSE],
        x=array(X,dim=c(dim(X),1),dimnames=list(variable=rownames(X))),
        times=times[nt+1],
        params=tparams,
        log=FALSE,
        .gnsi=gnsi
      )
      
      ## Apply tolerance to decide which ensemble fail
      index.missing <- weights <= tol
      weights[index.missing] <- 0
      weights[is.na(weights)] <- 0
      loglik[nt] <- ifelse(sum(weights) == 0, 
                           #log(1e-17),
                           #log(.Machine$double.xmin), 
                           -Inf,
                           log(mean(weights)))
      
      
      if(  sum(weights) == 0){
        nfail <- nfail + 1
      }
    }
    
    
    #####* Compute weighted mean at last timestep if no faiure, otherwise take unweighted mean *#####
    if(nt == ntimes){
      weights <- pomp::dmeasure(
        object,
        y=object@data[,nt,drop=FALSE],
        x=array(X,dim=c(dim(X),1),dimnames=list(variable=rownames(X))),
        times=times[nt+1],
        params=tparams,
        log=FALSE,
        .gnsi=gnsi
      )
      
      ##*Apply tolerance
      index.missing <- weights <= tol
      weights[index.missing] <- 0
      
      ##*Normalize
      weights <- weights/sum(weights)
      
      ##*Remove NA, aovid all 0 since all 0 will produce NaN
      weights[is.na(weights)] <- 0
      
      if(sum(weights) != 0){
        coef(object,transform=TRUE) <- apply(X = params,
                                             MARGIN = 1,
                                             FUN = weighted.mean,
                                             w = weights)
      }
      else{
        coef(object,transform=TRUE) <- apply(X = params,
                                             MARGIN = 1,
                                             FUN = mean)
      }
    }
    
    ## Forecast 3.0 Get the extended hidden state (hidden state + parameters which can change)
    x.e <- rbind(X,params[npar.names,,drop=F])
    
    #####* Update states in EAKF
    pm <- rowMeans(x.e)  ## prediction mean of the hidden states
    Y[,] <- apply(rbind(X,tparams[npar.names,]),2,h) ## ensemble of forecasts
    

    ## Approximate the measurement noise by function R 
    ## and predicted mean parameters and hidden states 
    # R.num <- R(rowMeans(rbind(X,tparams[npar.names,])))
    R.num <- R(rowMeans(rbind(X,tparams[,])))
    if(any(diag(R.num) < 0) || any(is.na(R.num)) ){
      stop('Measurement error is uncomputable')
    }
    ## Here we assign R a small value since in Kalman algorithm, since we have:
    ##              fv <- tcrossprod(Y)/(Np-1)+R,    (Here Y can be 0)
    ##                svdS <- svd(fv,nv=0) 
    ##      Kt <- svdS$u%*%(crossprod(svdS$u,vyx)/svdS$d)
    ##            We need to avoid svdS$d = 0 to make /svdS$d work
    
    #R.num <- ifelse( any(diag(R.num) == 0),diag(tolerance),R.num)
    ## Replace the variance which is 0 with the tolerance
    if(length(R.num) > 1){
      if(any(diag(R.num) == 0)){
        index <- diag(R.num) == 0
        diag(R.num)[index] <- tolerance
      }

    }else{
      if(R.num == 0){
        R.num <- tolerance
      }
    }

    sqrtR <- t(chol(R.num))
    

    
    ym <- rowMeans(Y)                  ## forecast mean
    
    x.e <- x.e-pm
    Y <- Y-ym
    
    fv <- tcrossprod(Y)/(Np-1)+R.num    ## forecast variance
    
    vyx <- tcrossprod(Y,x.e)/(Np-1)   ## forecast/state covariance
    
    
    svdS <- svd(fv,nv=0)            ## singular value decomposition
    
    Kt <- svdS$u%*%(crossprod(svdS$u,vyx)/svdS$d) ## transpose of Kalman gain
    Ek <- sqrtR%*%matrix(rnorm(n=nobs*Np),nobs,Np) ## artificial noise
    resid <- y[,nt]-ym
    
    x.e <- x.e+pm+crossprod(Kt,resid-Y+Ek) ## updated hidden state
    
    ## break the extended hidden to changable parameters and hidden state
    x <- x.e[nvar.names,]
    params[npar.names,] <- x.e[npar.names,]
    
    ## compute the likelihood using the extended hidden state using matrix
    #loglik[nt] <- sum(dnorm(x=crossprod(svdS$u,resid),mean=0,sd=sqrt(svdS$d),log=TRUE))
    
  }
  
  return(new(
    "pfilterd_pomp",
    as(object,"pomp"),
    paramMatrix=params,
    cond.loglik=loglik,
    indices=.indices,
    Np=Np,
    nfail = as.integer(nfail),
    tol=tol,
    loglik=sum(loglik)
  ))}



mif.enkf <- function(object, Nmif, rw.sd,
                       cooling.fraction.50,
                       ## h,R is for Enkf
                       h,R,
                       Np, max.fail = Inf,
                       ..., verbose = FALSE,check.fail = FALSE,
                       cooling.type = "geometric",tol = 0,
                       .ndone = 0L, .indices = integer(0), .paramMatrix = NULL,
                       .gnsi = TRUE,tolerance = 1e-100){
  
  ### Default initialization in pomp package 
  verbose <- as.logical(verbose)
  check.fail <- as.logical(check.fail)
  object <- pomp(object,verbose=verbose)
  
  gnsi <- as.logical(.gnsi)
  
  if (length(Nmif) != 1 || !is.numeric(Nmif) || !is.finite(Nmif) || Nmif < 1)
    pStop_(sQuote("Nmif")," must be a positive integer.")
  Nmif <- as.integer(Nmif)
  
  if (is.null(.paramMatrix)) {
    start <- coef(object)
  } else {
    start <- apply(.paramMatrix,1L,mean)
  }
  
  ntimes <- length(pomp::time(object))
  
  if (is.null(Np)) {
    pStop_(sQuote("Np")," must be specified.")
  } else if (is.function(Np)) {
    Np <- tryCatch(
      vapply(seq_len(ntimes),Np,numeric(1)),
      error = function (e) {
        pStop_("if ",sQuote("Np"),
               " is a function, it must return a single positive integer.")
      }
    )
  } else if (!is.numeric(Np)) {
    pStop_(sQuote("Np"),
           " must be a number, a vector of numbers, or a function.")
  }
  
  if (!all(is.finite(Np)) || any(Np <= 0))
    pStop_(sQuote("Np")," must be a positive integer.")
  
  Np <- as.integer(Np)
  
  if (missing(rw.sd))
    pStop_(sQuote("rw.sd")," must be specified!")
  
  ###### Create the list of perturbed noise, use perturbn.kernel.sd function ######
  rw.sd <- perturbn.kernel.sd(rw.sd,time=time(object),paramnames=names(start))
  
  if (missing(cooling.fraction.50))
    pStop_(sQuote("cooling.fraction.50")," is a required argument.")
  if (length(cooling.fraction.50) != 1 || !is.numeric(cooling.fraction.50) ||
      !is.finite(cooling.fraction.50) || cooling.fraction.50 <= 0 ||
      cooling.fraction.50 > 1)
    pStop_(sQuote("cooling.fraction.50")," must be in (0,1].")
  cooling.fraction.50 <- as.numeric(cooling.fraction.50)
  
  ###### Choose cooling method and time, use mif2.cooling function ######
  cooling.fn <- mif2.cooling(
    type=cooling.type,
    fraction=cooling.fraction.50,
    ntimes=length(time(object))
  )
  
  if (is.null(.paramMatrix)) {
    paramMatrix <- array(data=start,dim=c(length(start),Np[1L]),
                         dimnames=list(variable=names(start),rep=NULL))
  } else {
    paramMatrix <- .paramMatrix
  }
  
  traces <- array(dim=c(Nmif+1,length(start)+1+
                          1),
                  dimnames=list(iteration=seq.int(.ndone,.ndone+Nmif),
                                variable=c("loglik","nfail",names(start))))
  traces[1L,] <- c(NA,NA,start)
  
  pomp::pompLoad(object,verbose=verbose)
  on.exit(pompUnload(object,verbose=verbose))
  
  ## change the parameter scale, change it back in the loop of enkf
  paramMatrix <- partrans(object,paramMatrix,dir="toEst",
                          .gnsi=gnsi)
  
  ###### iterate the filtering of EnKF #######
  
  ## iterate the filtering
  for (n in seq_len(Nmif)) {
    pfp <- enkf.internal(
      object=object,
      params=paramMatrix,
      Np= Np,
      R = R,
      h = h,
      ## n is the number in this loop
      mifiter=.ndone+n,
      cooling.fn=cooling.fn,
      rw.sd=rw.sd,
      check.fail = check.fail,
      tol=tol,
      max.fail=max.fail,
      verbose=verbose,
      .indices=.indices,
      .gnsi=gnsi,
      tolerance
    )
    
    gnsi <- FALSE
    
    paramMatrix <- pfp@paramMatrix
    traces[n+1,-c(1L,2L)] <- coef(pfp)
    traces[n,1L] <- pfp@loglik
    traces[n,2L] <- pfp@nfail
    .indices <- pfp@indices
    
    if (verbose) cat("mif2 EnKF iteration",n,"of",Nmif,"completed\n")
    
  }
  pfp@paramMatrix <- partrans(object,paramMatrix,dir="fromEst",.gnsi=gnsi)
  new(
    "mif2d_pomp",
    pfp,
    Nmif=Nmif,
    rw.sd=rw.sd,
    cooling.type=cooling.type,
    cooling.fraction.50=cooling.fraction.50,
    traces=traces
  )}








#################################
##########*** KF ***###########
#################################
kf.internal <- function(object, params, Np, mifiter, rw.sd, cooling.fn,
                        tol = 0, max.fail = Inf, verbose,tolerance,
                        A,Q,C,R,check.fail,
                        .indices = integer(0),
                        .gnsi = TRUE){
  #####* General initialization *#####
  
  ## tol is the tolerance used to test the filtering 
  ## failure which defined in Particle Filter
  tol <- as.numeric(tol)
  
  gnsi <- as.logical(.gnsi)
  verbose <- as.logical(verbose)
  mifiter <- as.integer(mifiter)
  Np <- as.integer(Np)
  
  if (length(tol) != 1 || !is.finite(tol) || tol < 0)
    pStop_(sQuote("tol")," should be a small nonnegative number.")
  
  do_ta <- length(.indices)>0L
  if (do_ta && length(.indices)!=Np[1L])
    pStop_(sQuote(".indices")," has improper length.")
  
  times <- time(object,t0=TRUE)
  ntimes <- length(times)-1
  
  nfail <- 0
  
  #####* Kalman initialization *#####
  
  loglik <- numeric(ntimes)
  y <- obs(object) ## observation data
  nobs <- nrow(y) ## dim of observation 
  Y <- array(dim=c(nobs,Np),dimnames=list(variable=rownames(y),rep=NULL)) ## ensemble of measurements
  X <- rinit(object,params=coef(object),nsim=Np) ## ensemble of hidden states
  nvar.names <- rownames(rinit(object,params=coef(object),nsim=Np)) ## name of hidden state
  nvar <- nrow(X) ## dim of hidden state 
  
  ### Create extended observation matrix C.e
  
  ## names of parameters which do random walk
  npar.names <- names(cooling.fn(1,mifiter)$alpha*rw.sd[,1]) 
  
  ## dim of parameters which do random walk
  npar <- length(npar.names)
  ## Create extended measurement matrix
  # C.p <- matrix(rep(0,npar), nobs, npar, byrow = TRUE)
  # colnames(C.p) <- npar.names
  # C.e <- cbind(C,C.p)
  C.e <- C
  ## Create extended process matrix
  A.e <- diag( rep(1,npar+nvar) )
  colnames(A.e) <- c(colnames(A),npar.names)
  rownames(A.e) <- c(colnames(A),npar.names)
  A.e[1:dim(A)[1],1:dim(A)[2]] <- A
  ## Create forecast variance
  fv <- matrix(0,nvar+npar,nvar+npar)
  colnames(fv) <- c(colnames(A),npar.names)
  rownames(fv) <- c(colnames(A),npar.names)
  ## Create process noise
  Q.e <- matrix(rep(0,npar+nvar), npar+nvar, npar+nvar, byrow = TRUE)
  colnames(Q.e) <- c(colnames(A),npar.names)
  rownames(Q.e) <- c(colnames(A),npar.names)
  Q.e[1:dim(Q)[1],1:dim(Q)[2]] <- Q

  ri <- solve(R)
  cric <- crossprod(C.e,ri)%*%C.e ## cric = t(C.e) %*% R^-1 %% C.e
  
  ## Number of filtering failure
  nfail <- 0
  
  
  #####* Start of the main loop *#####
  for (nt in seq_len(ntimes)) {
    
    ## decide parameter variance at nt
    pmag <- cooling.fn(nt,mifiter)$alpha*rw.sd[,nt]
    
    ## Forecast 1.0 : Add noise to parameters
    for (par.name in npar.names) {
      params[par.name,] <- rnorm(n = dim(params)[2], 
                                 mean = params[par.name,],
                                 sd = pmag[par.name])
    }
    Q.e[npar.names,npar.names] <- diag(pmag[npar.names])^2
    
    
    ## Forecast 1.1 : Change parameters back to the original scale
    tparams <- pomp::partrans(object,params,dir="fromEst",.gnsi=gnsi)
    
    
    ## Forecast 1.2  Get initial states (if nt = 1) 
    if (nt == 1L) {
      x <- pomp::rinit(object,params=tparams)
    }
    
    ## Forecast 2.1 Advance the state variables according to the process model
    X[,] <- pomp::rprocess(
      object,
      x0=x,
      t0=times[nt],
      times=times[nt+1],
      params=tparams,
      .gnsi=gnsi
    )
    
    #####* Check if filtering failure happens, not necessary in algorithm *#####
    if(check.fail){
      weights <- pomp::dmeasure(
        object,
        y=object@data[,nt,drop=FALSE],
        x=array(X,dim=c(dim(X),1),dimnames=list(variable=rownames(X))),
        times=times[nt+1],
        params=tparams,
        log=FALSE,
        .gnsi=gnsi
      )
      ## Apply tolerance to decide which ensemble fail
      index.missing <- weights <= tol
      weights[index.missing] <- 0
      loglik[nt] <- ifelse(mean(weights) == 0, 
                           #log(1e-17),
                           #log(.Machine$double.xmin), 
                           -Inf,
                           log(mean(weights)))
      
      
      if(  sum(weights) == 0  ){
        nfail <- nfail + 1
      }
    }
    
    #####* Compute weighted mean at last timestep *#####
    if(nt == ntimes){
      coef(object,transform=TRUE) <- params[,1]
    }
    
    ## Forecast 3.0 Get the extended hidden state (hidden state + parameters which can change)
    x.e <- rbind(X,t(t(params[names(pmag),])))
    
    #####* Updare steps of KF *#####
    ##* Prediction Step
    ## Compute the prediction mean of the transformed parameters and hidden states 
    predMeans <- pm <- x.e  ## prediction mean
    pv <- A.e%*%tcrossprod(fv,A.e)+Q.e       ## prediction variance pv = A.e %*% fv %*% t(A.e) + Q.e
    ##* Update Step
    svdV <- svd(pv,nv=0)
    
    resid <- y[,nt]-C.e%*%pm              # forecast error
    w <- tcrossprod(C.e%*%pv,C.e)+R        # forecast variance
    
    svdW <- svd(w,nv=0)
    
    # loglik[nt] <- sum(dnorm(x=crossprod(svdW$u,resid),mean=0,
    #                            sd=sqrt(svdW$d),log=TRUE))
    
    pvi <- svdV$u%*%(t(svdV$u)/svdV$d) # prediction precision pv^-1
    fvi <- pvi+cric                    # filter precision pv^-1 + t(C.e) %*% R^-1 %% C.e
    
    svdv <- svd(fvi,nv=0) 
    
    fv <- svdv$u%*%(t(svdv$u)/svdv$d)  # filter variance pv + 
    K <- fv%*%crossprod(C.e,ri)          # Kalman gain 
    # (pv^-1 + t(C.e) %*% R^-1 %% C.e)) %*% t(C.e) %*% R^-1
    #K <- tcrossprod(x.e,C.e%*%pm) %*% solve(C.e %*% pv %*% t(C.e) + R)
    
    fm <- pm+K%*%resid  # filter mean
    ###forecast <- C.e %*% pm
    
    ## break the extended hidden stateto changable parameters and hidden state
    x[,] <- fm[nvar.names,]
    params[npar.names,] <- fm[npar.names,]
    
  }
  return(new(
    "pfilterd_pomp",
    as(object,"pomp"),
    paramMatrix=params,
    cond.loglik=loglik,
    indices=.indices,
    Np=Np,
    nfail = as.integer(nfail),
    tol=tol,
    loglik=sum(loglik)
  ))}



mif.kf <- function(object, Nmif, rw.sd,
                   cooling.fraction.50,
                   A,Q,C,R,
                   Np = 1, max.fail = Inf,
                   ..., verbose = FALSE,check.fail = FALSE,
                   cooling.type = "geometric",tol = 0,
                   .ndone = 0L, .indices = integer(0), .paramMatrix = NULL,
                   .gnsi = TRUE,tolerance = 1e-200){
  
  ### Default initialization in POMP package
  verbose <- as.logical(verbose)
  check.fail <- as.logical(check.fail)
  object <- pomp(object,verbose=verbose)
  
  gnsi <- as.logical(.gnsi)
  
  if (length(Nmif) != 1 || !is.numeric(Nmif) || !is.finite(Nmif) || Nmif < 1)
    pStop_(sQuote("Nmif")," must be a positive integer.")
  Nmif <- as.integer(Nmif)
  
  if (is.null(.paramMatrix)) {
    start <- coef(object)
  } else {  
    start <- apply(.paramMatrix,1L,mean)
  }
  
  ntimes <- length(pomp::time(object))
  
  if (is.null(Np)) {
    pStop_(sQuote("Np")," must be specified.")
  } else if (is.function(Np)) {
    Np <- tryCatch(
      vapply(seq_len(ntimes),Np,numeric(1)),
      error = function (e) {
        pStop_("if ",sQuote("Np"),
               " is a function, it must return a single positive integer.")
      }
    )
  } else if (!is.numeric(Np)) {
    pStop_(sQuote("Np"),
           " must be a number, a vector of numbers, or a function.")
  }
  
  
  if (!all(is.finite(Np)) || any(Np <= 0))
    pStop_(sQuote("Np")," must be a positive integer.")
  
  Np <- as.integer(Np)
  
  if (missing(rw.sd))
    pStop_(sQuote("rw.sd")," must be specified!")
  
  ###### Create list of perturbed noise, use perturbn.kernel.sd function ######
  rw.sd <- perturbn.kernel.sd(rw.sd,time=time(object),paramnames=names(start))
  
  if (missing(cooling.fraction.50))
    pStop_(sQuote("cooling.fraction.50")," is a required argument.")
  if (length(cooling.fraction.50) != 1 || !is.numeric(cooling.fraction.50) ||
      !is.finite(cooling.fraction.50) || cooling.fraction.50 <= 0 ||
      cooling.fraction.50 > 1)
    pStop_(sQuote("cooling.fraction.50")," must be in (0,1].")
  cooling.fraction.50 <- as.numeric(cooling.fraction.50)
  
  ###### Choose cooling method and time, use mif2.cooling function ######
  cooling.fn <- mif2.cooling(
    type=cooling.type,
    fraction=cooling.fraction.50,
    ntimes=length(time(object))
  )
  
  if (is.null(.paramMatrix)) {
    paramMatrix <- array(data=start,dim=c(length(start),Np[1L]),
                         dimnames=list(variable=names(start),rep=NULL))
  } else {
    paramMatrix <- .paramMatrix
  }
  
  traces <- array(dim=c(Nmif+1,length(start)+1+
                          1),
                  dimnames=list(iteration=seq.int(.ndone,.ndone+Nmif),
                                variable=c("loglik","nfail",names(start))))
  traces[1L,] <- c(NA,NA,start)
  
  pomp::pompLoad(object,verbose=verbose)
  on.exit(pompUnload(object,verbose=verbose))
  
  ## change the parameter scale, change them back in the loop of enkf
  paramMatrix <- partrans(object,paramMatrix,dir="toEst",
                          .gnsi=gnsi)
  
  ###### iterate the filtering of EAKF #######
  ## iterate the filtering
  for (n in seq_len(Nmif)) {
    pfp <- kf.internal(
      object=object,
      params=paramMatrix,
      Np= Np,
      A = A,
      Q = Q,
      R = R,
      C = C,
      mifiter=.ndone+n,
      cooling.fn=cooling.fn,
      rw.sd=rw.sd,
      tol=tol,
      check.fail = check.fail,
      max.fail=max.fail,
      verbose=verbose,
      .indices=.indices,
      .gnsi=gnsi,
      tolerance = tolerance
    )
    
    gnsi <- FALSE
    
    paramMatrix <- pfp@paramMatrix
    traces[n+1,-c(1L,2L)] <- coef(pfp)
    traces[n,1L] <- pfp@loglik
    traces[n,2L] <- pfp@nfail
    .indices <- pfp@indices
    
    if (verbose) cat("mif2 EnKF iteration",n,"of",Nmif,"completed\n")
    
  }
  pfp@paramMatrix <- partrans(object,paramMatrix,dir="fromEst",.gnsi=gnsi)
  new(
    "mif2d_pomp",
    pfp,
    Nmif=Nmif,
    rw.sd=rw.sd,
    cooling.type=cooling.type,
    cooling.fraction.50=cooling.fraction.50,
    traces=traces
  )}














###### mif2.cooling function ######
mif2.cooling <- function (type, fraction, ntimes) {
  switch(
    type,
    geometric={
      factor <- fraction^(1/50)
      function (nt, m) {
        alpha <- factor^(nt/ntimes+m-1)
        list(alpha=alpha,gamma=alpha^2)
      }
    },
    hyperbolic={
      if (fraction < 1) {
        scal <- (50*ntimes*fraction-1)/(1-fraction)
        function (nt, m) {
          alpha <- (1+scal)/(scal+nt+ntimes*(m-1))
          list(alpha=alpha,gamma=alpha^2)
        }
      } else {
        function (nt, m) {
          list(alpha=1,gamma=1)
        }
      }
    }
  )
}

###### perturbn.kernel.sd function ######

perturbn.kernel.sd <- function(rw.sd, time, paramnames){
  
  if (is.matrix(rw.sd)) return(rw.sd)
  if (is(rw.sd,"safecall")) {
    enclos <- rw.sd@envir
    rw.sd <- as.list(rw.sd@call)[-1L]
  } else {
    pStop_(sQuote("rw.sd")," should be specified using the ",sQuote("rw.sd"),
           " function. See ",sQuote("?mif2"),".")
  }
  if (is.null(names(rw.sd)) | any(names(rw.sd)==""))
    pStop("rw.sd","parameters must be referenced by name.")
  if (!all(names(rw.sd) %in% paramnames)) {
    unrec <- names(rw.sd)[!names(rw.sd) %in% paramnames]
    pStop_("the following parameter(s), ",
           "given random walks in ",sQuote("rw.sd"),", are not present in ",
           sQuote("params"),": ",paste(sapply(unrec,sQuote),collapse=","),".")
  }
  ivp <- function (sd, lag = 1L) {
    sd*(seq_along(time)==lag)
  }
  sds <- lapply(rw.sd,eval,envir=list(time=time,ivp=ivp),enclos=enclos)
  for (n in names(sds)) {
    len <- length(sds[[n]])
    if (len==1) {
      sds[[n]] <- rep(sds[[n]],length(time))
    } else if (len!=length(time)) {
      pStop_(sQuote("rw.sd")," spec for parameter ",sQuote(n),
             " does not evaluate to a vector of the correct length (",
             sQuote("length(time(object))"),"=",length(time),").")
    }
  }
  do.call(rbind,sds)
}
