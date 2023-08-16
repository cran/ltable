
SetUpLogLin<-function (formula, data, contrasts=NULL )
{
  cal <- match.call()
  
  if (missing(data)) 
    data <- environment(formula)
  
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data"), names(mf), 0L)
  check.NA<-as.character(formula[2])
  stopifnot("\nThere are missing (NA) values in counts (dependent variable)\n"= !is.na(eval(parse(text=paste0("data","$",check.NA))))) 
  mf <- mf[c(1L, m)]
  mf$drop.unused.levels <- TRUE
  mf[[1L]] <- quote(stats::model.frame)
  mf <- eval(mf, parent.frame())
  mt <- attr(mf, "terms")
  Y <- model.response(mf, "any")
  if (length(dim(Y)) == 1L) {
    nm <- rownames(Y)
    dim(Y) <- NULL
    if (!is.null(nm)) 
      names(Y) <- nm
  }
  X <- if (!is.empty.model(mt)) 
    model.matrix(mt, mf, contrasts)
  else matrix(, NROW(Y), 0L)
  
  dim(Y)<-c(length(Y),1)
  list(Y=Y,X=X)
}



funcbeta<-function(beta, X, j, lambda, cphi){
-cphi*( sum(X[,j, drop=FALSE])*beta[j] + sum(lambda*exp(-X%*%beta)) ) 
                                            }

funcphi<-function(cphi, Xbeta, lambda, a, b){
  stopifnot("\nXbeta OR lambda is not column vector" = ncol(Xbeta) ==1 && ncol(lambda) ==1 ) 
  stopifnot("\nlength of lambda is different from # of counts" = nrow(Xbeta) == nrow(lambda)) 
  n<-nrow(Xbeta)
  -n*lgamma(cphi) + cphi*(n*log(cphi) - sum(Xbeta - log(lambda) + lambda*exp(-Xbeta))) +
    (a-1)*log(cphi) -b*cphi
}

zbetadraw<-function(betaold, X, j, lambda, cphi){
  zold <- funcbeta(betaold, X, j, lambda, cphi) 
  zold + log(runif(1))
}

zphidraw<-function(cphiold, Xbeta, lambda, a, b){
  zold <- funcphi(cphiold, Xbeta, lambda, a, b) 
  zold + log(runif(1))
}

slicebeta<-function(betaold, X, j, lambda, cphi, XLB, XUB){
  Z<-zbetadraw(betaold, X, j, lambda, cphi)
  
  betanew<-betaold
  u <- runif(1);
  W<-1;
  L <- max(XLB[j], betaold[j] - W*u ) ;
  R <- min(XUB[j], L+W);
  
  while (L > XLB[j]){
    betanew[j] <- L;
    fxl <- funcbeta (betanew, X, j, lambda, cphi) ;
    if (fxl <= Z) {break} ;
    L <- max(XLB[j], L - W) ;
  }
  
  while(R < XUB[j]){
    betanew[j] <- R;
    fxr <- funcbeta (betanew, X, j, lambda, cphi) ;
    if (fxr <= Z) {break} ;  
    R <- min(XUB[j], R + W);
  }
  
  fxsim<- Z-1;
  while (fxsim < Z){
    u <- runif(1);
    betanew[j] <- L + u*(R - L) ;
    fxsim <- funcbeta (betanew, X, j, lambda, cphi) ;
    if(betanew[j]> betaold[j]) {
      R <- betanew[j]
    } else {
      L <- betanew[j]
    }
  }
  return (betanew[j]);
}

slicephi<-function(phiold, Xbeta, lambda, a, b, LB, UB){
  Z<-zphidraw(phiold, Xbeta, lambda, a, b)
  u <- runif(1);
  W<-3;
  L <- max(LB, phiold - W*u ) ;
  R <- min(UB, L+W);
  
  while (L > LB){
    phisim <- L;
    fxl <- funcphi(phisim, Xbeta, lambda, a, b);
    if ( fxl <= Z) {break} ;
    L <- max(LB, L - W) ;
  }
  
  while(R < UB){
    phisim <- R;
    fxr <- funcphi(phisim, Xbeta, lambda, a, b) ;
    if ( fxr <= Z) {break} ;  
    R <- min(UB, R + W);
  }
  
  fxsim <- Z-1;
  while ( fxsim < Z){
    u <- runif(1);
    phisim <- L + u*(R - L) ;
    fxsim <- funcphi(phisim, Xbeta, lambda, a, b) ;
    if(phisim > phiold) {
      R <- phisim
    } else {
      L <- phisim
    }
  }
  return (phisim);
}



MCLogLin<-function(formula, data, offset, contrasts=NULL,  XLB=-100, XUB=100, a=0.1, b=0.1, DIC=FALSE, pcov=FALSE, draw=10000, burnin=3000 )
{
  
  stopifnot("\npars \"a\" and \"b\" should excced zero\n"= a>0 && b>0) 
  stopifnot("\npar \"draw\" should exceed par \"burnin\" at least by 3000\n"= draw-burnin > 500)
  stopifnot("\npar \"DIC\" is logical (TRUE or FALSE)\n"= is.logical(DIC))
  stopifnot("\npar \"pcov\" is logical (TRUE or FALSE)\n"= is.logical(pcov))
  
  
  ZData<-SetUpLogLin(formula, data, contrasts=NULL )
  lambda<-ZData$Y
  nlambda<-dim(lambda)[1]
  nbeta<-dim(ZData$X)[2]
  beta<-array(data=rep_len(0, nbeta), dim=c(nbeta,1))
  cphi<-1
  
  if (length(XLB) == 1) {XLB<-array(XLB, dim=c(nbeta,1))}
  else if (length(XLB) == nbeta) {dim(XLB)<-c(nbeta,1)}
  else {stop(simpleError("XLB is of wrong length. \nIt should be of length 1 or regr. effects number\n"))}
  if (length(XUB) == 1) {XUB<-array(XUB, dim=c(nbeta,1))}
  else if (length(XUB) == nbeta) {dim(XUB)<-c(nbeta,1)}
  else {stop(simpleError("XUB is of wrong length. \nIt should be of length 1 or regr. effects number\n"))}
  if (any (XLB >= XUB)) {stop(simpleError("Lower boundaries should be smaller than upper\n"))}
  else {}
  
  drawset<-matrix(data=NA, nrow=draw, ncol=nlambda+nbeta+1)
  colnames(drawset)<-c(paste0("lambda",1:nlambda),colnames(ZData$X),"phi")
  
  vrgamma<-Vectorize(rgamma, vectorize.args=c("shape","rate"))
  
  ##DRAWS BEGINNING
  
  if(!missing(offset)) 
  {
    Ylength<-length(ZData$Y)
    OFFSET=1
    ar<-eval(substitute(data$offset), envir=environment(MCLogLin))
    #exposure<-array(as.vector(ar, mode="double"), dim=c(Ylength,1))
    exposure<-array(ar, dim=c(Ylength,1))
    
  }
  else 
  {
    OFFSET=0 
    exposure<-0
  }
  
  for(i in 1:draw){
    Xbeta<-ZData$X%*%beta
    
    if(OFFSET==0) 
      lambda<- vrgamma(1, ZData$Y+cphi, 1+cphi*exp(-Xbeta))
    else if(OFFSET==1) 
      lambda<- vrgamma(1, ZData$Y+cphi, exposure+cphi*exp(-Xbeta))  
                        
    
    dim(lambda)<-c(nlambda,1)
    drawset[i,1:nlambda]<-lambda
    

      for(jj in 1:nbeta) {
        beta[jj]<-slicebeta(betaold=beta, X=ZData$X, j=jj, lambda=lambda, cphi=cphi, XLB=XLB, XUB=XUB)
                          }
    
    drawset[i, (nlambda+1):(nlambda+nbeta)]<-beta
    Xbeta<-ZData$X%*%beta
    cphi<-slicephi(phiold=cphi, Xbeta=Xbeta, lambda=lambda, a=a, b=b, LB=0, UB=100)
    drawset[i, nlambda+nbeta+1]<-cphi 
  }
  ##DRAWS END
  
  
  errors<-array(NA, c(nbeta+1,1))
  COVAR<-0
  lambda<-colMeans(drawset[burnin:draw,1:nlambda])
  cphi<-mean(drawset[burnin:draw, nlambda+nbeta+1])
  if(length(exposure)<=1) {
    for(i in 1:nlambda){
      k<-lambda[i]/(1+lambda[i]/cphi)
      COVAR<-COVAR + k*(t(ZData$X[i, , drop = FALSE]) %*% ZData$X[i, , drop = FALSE])
    }
    
  } else {
    for(i in 1:nlambda){
      k<-lambda[i]*exposure[i]/(1+lambda[i]*exposure[i]/cphi)
      COVAR<-COVAR + k*(t(ZData$X[i, , drop = FALSE]) %*% ZData$X[i, , drop = FALSE])
    }
  } 
  
  flag<-0
  COVAR<-tryCatch(solve(COVAR), error=function(c) {
    flag<<-1; 
    x<-drawset[,(nlambda+1):(nlambda+nbeta+1)]
    colnames(x)<-colnames(drawset[,(nlambda+1):(nlambda+nbeta+1)])
    cov(x,x)
  } 
  )
  
  errors<-sqrt(diag(COVAR))      
  
  sumcphi<-function(cphi, J, y) {
    zv<-1/cphi
    if ((y[J]-1) > 0)
    {
      for(i in 1:y[J]){
        zv<-zv+1/(i+cphi)
      } 
    } else { 
      zv<-0
    }
    zv
  }
  
  vsumcphi<-Vectorize(sumcphi, vectorize.args = "J")
  const1<-array(1:nlambda,c(nlambda,1))/(1/cphi)^4
  const2<-log(1+lambda/cphi) - array(vsumcphi(cphi,1:nlambda, ZData$Y), c(nlambda,1))
  const3<-lambda/((1/cphi)^2*(1+lambda/cphi))
  errors[nbeta+1]<-cphi*sqrt(1/sum(const1*const2^2+const3))*cphi
  
  
  cat("\nCall:\n")
  print (match.call())
  
  cat("\n\nCoefficients:\n")
  
  l<-40L
  cat(sprintf("%*s %*s %*s %*s" , l+2, "Estimate", 10, "Std.Error", 10, "|z-score|", 11,"Pr(>|z|)\n")) 
  
  pars<-colMeans(drawset[burnin:draw,(nlambda+1):(nlambda+nbeta+1)])
  vnames<-colnames(drawset[,(nlambda+1):(nlambda+nbeta+1)])
  z_value<-abs(pars/errors)
  
  for(j in 1:(nbeta+1)) {
    if (nchar(vnames[j]) >= l-5) {
      spl<-splitname(vnames[j],l-5)
      cat(sprintf("%s%*.3e%*.3e%*.3e%*.3e",spl$sname, l+2-nchar(spl$sname), pars[j], 11, errors[j],
                  11, z_value[j], 11, pnorm(z_value[j], lower.tail = FALSE)*2),"\n")
      cat(spl$adding,"\n")
    } else{
      
      cat(sprintf("%s%*.3e%*.3e%*.3e%*.3e",vnames[j], l+2-nchar(vnames[j]), pars[j], 11, errors[j],
                  11, z_value[j], 11, pnorm(z_value[j], lower.tail = FALSE)*2),"\n")
    }
  }
  
  
  m_lambda<-array(colMeans(drawset[(burnin+1):draw, 1:nlambda]), c(nlambda,1))
  m_cphi<-mean(drawset[(burnin+1):draw, nlambda+nbeta+1])

  Devmean<-output(ZData$Y, m_lambda, m_cphi, nlambda, nbeta, COVAR, exposure, Xbeta)
  
  
  if(DIC ==TRUE)
  {
    Davg<-array(0, c(draw-burnin,1))
    const<-lgamma(ZData$Y+1)
   if(length(exposure)<=1)
      {
    for(i in (burnin+1):draw) {
      Davg[i-burnin]<--2*sum(lgamma(ZData$Y+drawset[i,nlambda+nbeta+1]) - const - lgamma(drawset[i,nlambda+nbeta+1]) +
                               drawset[i,nlambda+nbeta+1]*log(drawset[i,nlambda+nbeta+1]/(array(exp(drawset[i,(nlambda+1):(nlambda+nbeta)]*ZData$X), c(nlambda,1))+drawset[i,nlambda+nbeta+1])) +
                               ZData$Y*log(array(exp(drawset[i,(nlambda+1):(nlambda+nbeta)]*ZData$X), c(nlambda,1))/(array(exp(drawset[i,(nlambda+1):(nlambda+nbeta)]*ZData$X), c(nlambda,1))+drawset[i,nlambda+nbeta+1])
                               )) 
                              }
      } else {
      for(i in (burnin+1):draw) {
        Davg[i-burnin]<--2*sum(lgamma(ZData$Y+drawset[i,nlambda+nbeta+1]) - const - lgamma(drawset[i,nlambda+nbeta+1]) +
                                 drawset[i,nlambda+nbeta+1]*log(drawset[i,nlambda+nbeta+1]/(array(exp(drawset[i,(nlambda+1):(nlambda+nbeta)]*ZData$X)*exposure, c(nlambda,1))+drawset[i,nlambda+nbeta+1])) +
                                 ZData$Y*log(array(exp(drawset[i,(nlambda+1):(nlambda+nbeta)]*ZData$X)*exposure, c(nlambda,1))/(array(exp(drawset[i,(nlambda+1):(nlambda+nbeta)]*ZData$X)*exposure, c(nlambda,1))+drawset[i,nlambda+nbeta+1])
                                 )) 
                                }
      
              }
    meanDavg<-mean(Davg, na.rm = TRUE)
    DIC<- -2*Devmean + 2*(meanDavg - (-2*Devmean))
    cat("DIC related components: \n")      
    cat("DIC = ", DIC, "\n")
    cat("pD = ", meanDavg - (-2*Devmean), "\n")
    cat("meanDev = ", meanDavg, "\n")
    cat("Devmean = ", -2*Devmean, "\n")
    cat("\n\n")
  }    
  
  
  if (pcov == TRUE) 
  {
    cat("Covariance matrix: \n")   
    print(COVAR, digits=4)
    cat("\n")
    COR<-matrix(0, nrow=nbeta, ncol=nbeta)
    tryCatch ({for(i in 1:nbeta){ 
      for(j in 1:nbeta) {COR[i,j]<-COVAR[i,j]/(sqrt(COVAR[i,i])*sqrt(COVAR[j,j]))}}
      cat("\nCorrelation matrix: \n")
      print(COR, digits=4)    
    }, error=function(c) {message("Can't calculate correlation matrix.")})
  }
  
  if (flag==1) warning ("MCMC based errors are used")
  drawset
  
} 
##END OF MCLogLin 


splitname<-function(sname, l){
  cnt<-0
  while (nchar(sname) >= l) 
  {
    cnt<-cnt+1
    s000<-simplify2array(strsplit(sname, ":"))
    sname<-character(0)
    len<-length(s000)
    for (i  in 1: (len-cnt)) {
      sname<-paste0(sname, ":", s000[i])
    }
    sname<-paste0(sname,":")
    sname<-substr(sname,2,nchar(sname))
  }
  adding<-""
  for(i in 1:cnt) adding<-paste0(adding,s000[length(s000)-cnt+i])
  
  return (list(sname=sname, adding=adding))
}


output<-function(Y, lambda, cphi, nlambda, nbeta, COVAR, exposure, LP)
{
  cat("\n\nModel fit:\n")
  cat("MCMC fitting\n")
  cat("Samplers : Gibbs for expected counts, Slice for regr. coeff. and inv.var.par. phi\n")
  cat("Language: R\n")
  cat("Jacobian reciprocal condition number = ", rcond(solve(chol(COVAR)), triangular=TRUE), "\n")
  if(length(exposure)<=1)  
  cat("chisq/n = ", sum((Y-lambda)^2/(lambda+lambda^2/cphi))/nlambda, "\n")
  else {
  cat("chisq/n = ", sum((Y-lambda*exposure)^2/(lambda*exposure+(lambda*exposure)^2/cphi))/nlambda, "\n")
  }
  m_expLP<-mean(exp(LP))
  
  if(length(exposure)<=1) 
    {
    for(i in 1:nlambda){
    dev<-0
    dev0<-0
    if (Y[i] > 0) {
      dev<-dev+Y[i]*log(Y[i]/exp(LP[i])) - (Y[i]+cphi)*log((Y[i]+cphi)/(exp(LP[i])+cphi))
      dev0<-dev0+Y[i]*log(Y[i]/m_expLP) - (Y[i]+cphi)*log((Y[i]+cphi)/(m_expLP+cphi))
    } else { 
      dev<-dev+cphi*log(1+exp(LP[i])/cphi)
      dev0<-dev0+cphi*log(1+m_expLP/cphi)
    }
                        }
    }
  else ##with offset
  {
    m_exposure<-mean(exposure)
    for(i in 1:nlambda){
      dev<-0
      dev0<-0
      if (Y[i] > 0) {
        dev<-dev+Y[i]*log(Y[i]/(exp(LP[i])*exposure[i])) - (Y[i]+cphi)*log((Y[i]+cphi)/(exp(LP[i])*exposure[i]+cphi))
        dev0<-dev0+Y[i]*log(Y[i]/(m_expLP*m_exposure)) - (Y[i]+cphi)*log((Y[i]+cphi)/(m_expLP*m_exposure+cphi))
      } else { 
        dev<-dev+cphi*log(1+exp(LP[i])*exposure[i]/cphi)
        dev0<-dev0+cphi*log(1+m_expLP*m_exposure/cphi)
      }
    }
  }
  
  dev<-2*dev
  dev0<-2*dev0
  cat("Deviance= ", dev, "\n")
  cat("NULL Deviance= ", dev0, "\n") 
  
  if(length(exposure)<=1) 
    {
  LL<-sum(lgamma(Y+cphi) - lgamma(Y+1) - lgamma(cphi) +
            cphi*log(cphi/(exp(LP)+cphi)) + Y*log(exp(LP)/(exp(LP)+cphi)))
    }
  else
    {
    LL<-sum(lgamma(Y+cphi) - lgamma(Y+1) - lgamma(cphi) +
              cphi*log(cphi/(exp(LP)*exposure+cphi)) + Y*log(exp(LP)*exposure/(exp(LP)*exposure+cphi)))
    }
  cat("Log.likelihood= ", LL, "\n")
  cat("AIC(1) = ", -2*(LL-nbeta), "\n")  
  cat("AIC(n) = ", -2*(LL-nbeta)/nlambda, "\n")    
  cat("BIC = ", -2*LL + nbeta*log(nlambda), "\n")

  if(length(exposure)<=1) {dif<-Y-lambda}
  else {Y<-Y/exposure; 
        dif<-Y-lambda}
  Pearson<-dif/sqrt(lambda+lambda^2/cphi)
  StdDeviance<-(3*cphi*((1+Y/cphi)^(2/3) - (1+lambda/cphi)^(2/3)) + 3*(Y^(2/3) - lambda^(2/3)))/(2*(lambda+lambda^2/cphi)^(1/6))
 
   if(length(exposure)<=1)
  {
    cat("\n\nResiduals report (Y denotes Counts):\n")
    cat(sprintf("%*s %*s %*s %*s %*s %*s" , 2, "Row", 4, "Ovserved Y", 6, "Predicted Y", 6,"Raw Residual", 6,"Pearson Residual", 6,"Anscombe Residual\n")) 
    
    for (i in 1:nlambda) 
    {
      cat(sprintf("%*d %*d %*.3f %*.3f %*.3f %*.3f", 2, i, 10, Y[i], 11, lambda[i], 12, dif[i], 16, Pearson[i], 16, StdDeviance[i]), "\n")
    }
    cat("\n\n")
  }
  else ## if(!missing(offset))
  {
    cat("\n\nResiduals report ( R is individual risk):\n")
    cat(sprintf("%*s %*s %*s %*s %*s %*s" , 2, "Row", 4, "Ovserved R", 6, "Predicted R", 6,"Raw Residual", 6,"Pearson Residual", 6,"Anscombe Residual\n")) 
    
    expectedrisk<-lambda
    for (i in 1:nlambda) 
    {
      cat(sprintf("%*d %*.5f %*.5f %*.5f %*.4f %*.4f", 2, i, 10, Y[i], 11, expectedrisk[i], 12, dif[i], 16, Pearson[i], 16, StdDeviance[i]), "\n")
    }
    cat("\n\n")
  }   
  
  LL
}

