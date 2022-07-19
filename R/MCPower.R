setClass("powerClass", slots = c(varnames="character", effectsname="character", 
                                  cal="call", Ntotal="numeric", 
                                  estim= "list", power1="list", power2="list",
                                  power3="list", power4="list", power5="list",
                                  power6="list", power7="list", power8="list",
                                  power9="list", power10="list",power11="list"))

setMethod(
  f= "[", 
  signature="powerClass", 
  definition=function(x,i,j, drop){
    if(i=="varnames"){return(x@varnames)} else {} 
    if(i=="effectsname"){return(x@effectsname)}else {}
    if(i=="cal"){return(x@cal)}else {}  
    if(i=="Ntotal"){return(x@Ntotal[j])}else {}
    if(i=="estim"){return(x@estim[[j]])}  else {} 
    if(i=="power1"){return(x@power1[[j]])} else {}
    if(i=="power2"){return(x@power2[[j]])} else {}
    if(i=="power3"){return(x@power3[[j]])} else {}
    if(i=="power4"){return(x@power4[[j]])} else {}
    if(i=="power5"){return(x@power5[[j]])} else {}
    if(i=="power6"){return(x@power6[[j]])} else {}
    if(i=="power7"){return(x@power7[[j]])} else {}
    if(i=="power8"){return(x@power8[[j]])} else {}
    if(i=="power9"){return(x@power9[[j]])} else {}
    if(i=="power10"){return(x@power10[[j]])}else {}
    if(i=="power11"){return(x@power11[[j]])}else {}
  }  
)


setReplaceMethod(
  f="[",
  signature="powerClass",
  definition=function(x,i,j, value){ 
    if(i=="varnames"){x@varnames<-value}else{} 
    if(i=="effectsname"){return(x@effectsname<-value)}else {}
    if(i=="cal"){return(x@cal<-value)}else {}   
    if(i=="Ntotal"){return(x@Ntotal[j]<-value)}else {}  
    if(i=="estim"){x@estim[[j]]<-value} else{} 
    if(i=="power1"){x@power1[[j]]<-value} else{}
    if(i=="power2"){x@power2[[j]]<-value} else{}
    if(i=="power3"){x@power3[[j]]<-value} else{}
    if(i=="power4"){x@power4[[j]]<-value} else{}
    if(i=="power5"){x@power5[[j]]<-value} else{}
    if(i=="power6"){x@power6[[j]]<-value} else{}
    if(i=="power7"){x@power7[[j]]<-value} else{}
    if(i=="power8"){x@power8[[j]]<-value} else{}
    if(i=="power9"){x@power9[[j]]<-value} else{}
    if(i=="power10"){x@power10[[j]]<-value} else{}
    if(i=="power11"){x@power11[[j]]<-value} else{}
    validObject(x)
    return (x)
  }
)



MCPower<-function (formula, data, contrasts=NULL, XLB=-100, XUB=100, a=0.1, b=0.1, scale_min=1, scale_max=5, effect, p_alpha=0.05, draw=10000, burnin=3000)
{
  stopifnot("\npars a and b should excced zero\n"= a>0 && b>0)
  stopifnot("\npar draw should exceed par burnin at least by 3000\n"= draw-burnin > 500)
  cal <- match.call()
  
  if (missing(data)) 
    data <- environment(formula)

  vars<-unlist(strsplit(effect,"[*]"))
  lengthvars<-length(vars)
  stopifnot("\nIn effect= use pattern varnames delimited with *,\ne.g., effect=\"vname1*vname3*vname4\"\n"= all(vars %in% colnames(data)) 
            && lengthvars>1)
  stopifnot("\nYou have used same varname several times in effect\n"=length(vars)==length(unique(vars)) )
  
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
  else matrix( , NROW(Y), 0L)
  
  nbeta<-dim(X)[2]
  if (length(XLB) == 1) {XLB<-array(XLB, dim=c(nbeta,1))}
  else if (length(XLB) == nbeta) {dim(XLB)<-c(nbeta,1)}
  else {stop(simpleError("XLB is of wrong length. \nIt should be of length 1 or regr. effects number\n"))}
  if (length(XUB) == 1) {XUB<-array(XUB, dim=c(nbeta,1))}
  else if (length(XUB) == nbeta) {dim(XUB)<-c(nbeta,1)}
  else {stop(simpleError("XUB is of wrong length. \nIt should be of length 1 or regr. effects number\n"))}
  if (any (XLB >= XUB)) {stop(simpleError("Lower boundaries should be smaller than upper\n"))}
  else {}
  
  #**pick up related to effect columns out of design mat X in var whichcols**#  
  varsp<-"["
  for(i in seq_along(vars) ) varsp<-paste0(varsp, vars[i])
  varsp<-paste0(varsp, "].+:") 
  patterneffect=""
  for(i in seq_along(vars)-1 ) patterneffect<-paste0(patterneffect, varsp)
  patterneffect<-substr(patterneffect,1,nchar(patterneffect)-1)
  whichcols<-which( sapply(colnames(X), function(x) grepl(patterneffect,x)) )
  purify=integer(0)
  c<-sapply(colnames(X)[whichcols], function(x) strsplit(x, "[:]"))
  for(i in seq_along(whichcols)) if (length(c[[i]])!=lengthvars) purify<-c(purify,i) 
  if(length(purify)>0) whichcols<-whichcols[-purify]
  
  scalevec=seq(scale_min,scale_max, by=(scale_max - scale_min)/10)


    powerClassObject<-new("powerClass", varnames=colnames(X), effectsname=colnames(X)[whichcols], cal=cal, Ntotal=c(sum(Y), scale_min, scale_max) )
 
  
  for (i in 1:11) {
    Z<-cbind(Y*scalevec[i],X)

    powerClassObject["estim",i]<-LogLinEst(formula=formula, data=Z, ii=i, contrasts=contrasts, XLB=XLB, XUB=XUB, a=a, b=b, draw=draw, burnin=burnin)
                  betas000initdraw<-powerClassObject["estim",i]$betas
                  if(i==1) {phi000initdraw<-powerClassObject["estim",i]$phi} else {}                  
                  }
 cl<-qnorm(p_alpha, lower.tail =FALSE)
 nsim<-length(whichcols)
    vbetas<-matrix(NA,nrow=11, ncol=nsim)
    verrors<-matrix(NA,nrow=11, ncol=nsim)
    for(i in 1:11) {
    vbetas[i,]<-powerClassObject["estim",i]$betas[whichcols]
    verrors[i,]<-powerClassObject["estim",i]$errors[whichcols]
    }
    results<-mapply(FUN = simpower, vbetas[1,], verrors[1,], cl)
     for(i in 1:nsim){
      powerClassObject["power1",i] <- list(betas=results[1,i]$betas, z= results[2,i]$z, power=results[3,i]$power)
      }
    results<-mapply(FUN = simpower, vbetas[2,], verrors[2,], cl)
    for(i in 1:nsim){
      powerClassObject["power2",i] <- list(betas=results[1,i]$betas, z= results[2,i]$z, power=results[3,i]$power)
    }
    results<-mapply(FUN = simpower, vbetas[3,], verrors[3,], cl)
    for(i in 1:nsim){
      powerClassObject["power3",i] <- list(betas=results[1,i]$betas, z= results[2,i]$z, power=results[3,i]$power)
    }
    results<-mapply(FUN = simpower, vbetas[4,], verrors[4,], cl)
    for(i in 1:nsim){
      powerClassObject["power4",i] <- list(betas=results[1,i]$betas, z= results[2,i]$z, power=results[3,i]$power)
    }
    results<-mapply(FUN = simpower, vbetas[5,], verrors[5,], cl)
    for(i in 1:nsim){
      powerClassObject["power5",i] <- list(betas=results[1,i]$betas, z= results[2,i]$z, power=results[3,i]$power)
    }
    results<-mapply(FUN = simpower, vbetas[6,], verrors[6,], cl)
    for(i in 1:nsim){
      powerClassObject["power6",i] <- list(betas=results[1,i]$betas, z= results[2,i]$z, power=results[3,i]$power)
    }
    results<-mapply(FUN = simpower, vbetas[7,], verrors[7,], cl)
    for(i in 1:nsim){
      powerClassObject["power7",i] <- list(betas=results[1,i]$betas, z= results[2,i]$z, power=results[3,i]$power)
    }
    results<-mapply(FUN = simpower, vbetas[8,], verrors[8,], cl)
    for(i in 1:nsim){
      powerClassObject["power8",i] <- list(betas=results[1,i]$betas, z= results[2,i]$z, power=results[3,i]$power)
    }
    results<-mapply(FUN = simpower, vbetas[9,], verrors[9,], cl)
    for(i in 1:nsim){
      powerClassObject["power9",i] <- list(betas=results[1,i]$betas, z= results[2,i]$z, power=results[3,i]$power)
    }
    results<-mapply(FUN = simpower, vbetas[10,], verrors[10,], cl)
    for(i in 1:nsim){
      powerClassObject["power10",i] <- list(betas=results[1,i]$betas, z= results[2,i]$z, power=results[3,i]$power)
    }
    results<-mapply(FUN = simpower, vbetas[11,], verrors[11,], cl)
    for(i in 1:nsim){
      powerClassObject["power11",i] <- list(betas=results[1,i]$betas, z= results[2,i]$z, power=results[3,i]$power)
    }

 rm(results)
    return(powerClassObject)
}
#END OF MCPower

LogLinEst<-function(formula, data, ii, contrasts, ibetas, iphi, XLB, XUB, a, b, draw, burnin )
{
  Y<-data[,1, drop=FALSE]
  X<-data[,2:dim(data)[2], drop=FALSE]
  lambda<-Y
  nlambda<-dim(lambda)[1]
  nbeta<-dim(X)[2]

  if (ii == 1) {
  beta<-array(data=rep_len(0, nbeta), dim=c(nbeta,1)); cphi<-1
              } else {
  beta<-get("betas000initdraw", envir=parent.frame(1))  
  cphi<-get("phi000initdraw", envir=parent.frame(1))  
                      }
 
  drawset<-matrix(data=NA, nrow=draw, ncol=nlambda+nbeta+1)
  colnames(drawset)<-c(paste0("lambda",1:nlambda),colnames(X),"phi")
  
 vrgamma<-Vectorize(rgamma, vectorize.args=c("shape","rate"))

  ##DRAWS BEGINNING
  for(i in 1:draw){
    Xbeta<-X%*%beta
    lambda<- vrgamma(1, Y+cphi, 1+cphi*exp(-Xbeta))
    dim(lambda)<-c(nlambda,1)
    drawset[i,1:nlambda]<-lambda
    
     for(jj in 1:nbeta) {
       beta[jj]<-slicebeta(betaold=beta, X=X, j=jj, lambda=lambda, cphi=cphi, XLB=XLB, XUB=XUB)
     }
     drawset[i, (nlambda+1):(nlambda+nbeta)]<-beta
     Xbeta<-X%*%beta
if (ii == 1) 
  {
  cphi<-slicephi(phiold=cphi, Xbeta=Xbeta, lambda=lambda, a=a, b=b, LB=0, UB=100)
  } else {}
     drawset[i, nlambda+nbeta+1]<-cphi 
                    }
  ##DRAWS END

  errors<-array(NA, c(nbeta+1,1))
  COVAR<-0
  lambda<-colMeans(drawset[burnin:draw,1:nlambda])
  cphi<-mean(drawset[burnin:draw, nlambda+nbeta+1])
  for(i in 1:nlambda){
    phi<-ifelse(cphi>1, cphi, 1)
    k<-lambda[i]/(1+lambda[i]/phi)
    COVAR<-COVAR + k*(t(X[i, , drop = FALSE]) %*% X[i, , drop = FALSE])
  }

  COVAR<-tryCatch(qr.solve(COVAR, tol=1e-7), error=function(c) {
    stop("Sorry, can't proceed with singular Hessian matrix\n", call. = FALSE)
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
 const2<-log(1+lambda/cphi) - array(vsumcphi(cphi,1:nlambda, Y), c(nlambda,1))
 const3<-lambda/((1/cphi)^2*(1+lambda/cphi))
 errors[nbeta+1]<-cphi*sqrt(1/sum(const1*const2^2+const3))*cphi

if( ii==1)
{
  m_lam<-mean(lambda)
  for(i in 1:nlambda){
    dev<-0
    dev0<-0
    if (Y[i] > 0) {
      dev<-dev+Y[i]*log(Y[i]/lambda[i]) - (Y[i]+cphi)*log((Y[i]+cphi)/(lambda[i]+cphi))
      dev0<-dev0+Y[i]*log(Y[i]/m_lam) - (Y[i]+cphi)*log((Y[i]+cphi)/(m_lam+cphi))
    } else { 
      dev<-dev+cphi*log(1+lambda[i]/cphi)
      dev0<-dev0+cphi*log(1+m_lam/cphi)
    }
  }
  dev<-2*dev
  dev0<-2*dev0
 
  LL<-sum(lgamma(Y+cphi) - lgamma(Y+1) - lgamma(cphi) +
            cphi*log(cphi/(lambda+cphi)) + Y*log(lambda/(lambda+cphi)))
  Jacobian_rcnumber <- rcond(solve(chol(COVAR)), triangular=TRUE)  
  chi_sq <- sum((Y-lambda)^2/(lambda+lambda^2/cphi))/nlambda 
} else {}
 
  if(nbeta > 1) 
    {
  betas<-colMeans(drawset[,(nlambda+1):(nlambda+nbeta)])
    } else {
  betas<-mean(drawset[,(nlambda+1)])    
            }
  phi<-mean(drawset[,nlambda+nbeta+1])
  if(ii == 1)
  {
    out <- list(betas=betas, phi=phi, errors=errors, dev=dev, dev0=dev0, 
              LL=LL, Jacobian_rcnumber=Jacobian_rcnumber, chi_sq=chi_sq)
  } else {  
    out <- list(betas=betas, phi=phi, errors=errors)
          }
    out
}

simpower<- function(beta, error, cl) {
  betav<-numeric(0)
  zv<-numeric(0)
  powerv=numeric(0)
  for (k in 1:100){
    simbetas<-rnorm(50, mean=beta, sd=error)
    zscore<-abs(simbetas/error)
    powerz<-ifelse(zscore>cl, 1, 0)
    powerm<-mean(powerz)
    i=round(runif(1,1,50))
    betav<-c(betav,simbetas[i])
    zv=c(zv,zscore[i])
    powerv=c(powerv, powerm)
  }
  powerlist<-list(betas=betav, z=zv, power=powerv)
  return (powerlist)
}
