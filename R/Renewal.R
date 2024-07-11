library(MASS)
library(nloptr)
renewal<-function(x, StatesNum=10, Astate=NULL, nForStart = 3000, epsilon1=0.05, epsilon2=0.05)
{
  stopifnot("\nChain x must be numeric vector \n"= is.numeric(x)) 
  stopifnot("\nChain length (x) must surpass nForStart arg (3000 by default)\n"= length(x)>nForStart) 
##  nForStart arg limits cprob to find astart
  if(epsilon1>0.1 || epsilon2 > 0.1) warning("arg epsilon is too large, recommended value is less than 0.05") 

  
  if(StatesNum>100) {
    warning("Number of states is limited to 100,\n 100 is used
            arg Astate also changed to median", call. = FALSE)
    StatesNum<-100}
  missval<-which(is.na(x))
  if (length(missval)>0) {x <- x[-missval]} else {}
  x<-outlierfree(x)
  n<-length(x)
  breaksX<-seq(from=min(x), to=max(x), length.out=StatesNum)
  statesX<-findInterval(x, vec=breaksX)
  
  if (is.null(Astate)) Astate<-round(median(statesX))
  stopifnot("\narg Astate should be within the range of \nnumber of states given by arg StatesNum\n"= (Astate>0) && (Astate<=StatesNum)) 
  
  
  if(sum(statesX==StatesNum)<4) 
    { statesX[statesX==StatesNum]<-StatesNum-1; 
      if(Astate==StatesNum) Astate<-Astate-1;
      StatesNum<-StatesNum-1
    } 

  statesX<-checkForScarseStates(statesX)
  statesX<-recodeState(statesX)
  NSRemodelled<-length(unique(statesX))
  if(NSRemodelled!=StatesNum) 
    { warning(paste("Number of states reduced to", NSRemodelled), call. = FALSE)
    StatesNum<-NSRemodelled;
    Astate<-round(median(statesX));
    warning(paste("arg Astate also changed to median that is to", Astate), call. = FALSE) 
    }
  transProb<-t(TransProbM(statesX, StatesNum))
  cprob<-condprob(transProb, Astate, nForStart)
##  cprob is conditional prob. vector for state Astate
  if (cprob[1]==0) stop("transition probability of chosen state=0")
  
pnormal<-fitdistr(statesX,"normal")
probAstate<-dnorm(Astate, mean=pnormal$estimate[[1]], sd=pnormal$estimate[[2]])
restraints<-deviations(cprob, probAstate)
restraintsr1<-restraints[[1]]/restraints[[2]]
r1<-optimr1(restraintsr1)$solution 
M1<-optimr2(restraints[[1]], r1)$solution
lengthcprob<-length(cprob)

pNotUntill<-function(vec, untill) 
{
v<-1-vec[1]
for(i in 2:untill) v<-v*(1-vec[i])
v
}

cprobNotUntill<-sapply(1:lengthcprob, pNotUntill, vec=cprob)
cprobNotUntillratio<-ifelse(cprobNotUntill[1]>0, cprobNotUntill/cprobNotUntill[1], NULL)
if (is.null(cprobNotUntillratio)) stop("transition probability of chosen state=0")
r2<-optimr1(cprobNotUntillratio)$solution
##  slight correction for lower bound of optimization routine
if(r2 <= 1.02) r2<-r2+0.1 
M2<-optimr2(cprobNotUntill, r2)$solution

pNotUntillN<-function(vec, untill) 
{
  if(untill==1) return (1-vec[1])
  nLastNot <- (untill-1)
  v <- 1-vec[1]
  for(i in 2:nLastNot) v<-v*(1-vec[i])
  v<-v*vec[untill]
  v
}

cprobNotUntillN<-sapply(1:lengthcprob, pNotUntillN, vec=cprob)
cprobNotUntillNratio<-ifelse(cprobNotUntillN[1]>0, cprobNotUntill/cprobNotUntill[1], NULL)
if (is.null(cprobNotUntillNratio)) stop("transition probability of chosen state=0")
r3<-optimr1(cprobNotUntillNratio)$solution
if(r3 <= 1.02) r3<-r3+0.2 
M3<-optimr2(cprobNotUntillN, r3)$solution


epsilon<-Epsilon(r1,r2,r3, M1,M2,M3, probAstate, lengthcprob)
astart<-2*which(epsilon<epsilon1)[1]

if(length(which(statesX[astart:lengthcprob]==Astate))<4) Astate<-round(median(statesX[astart:lengthcprob]))

span<-calcSpan(statesX, astart, Astate, epsilon2)
c(astart,span)
}

outlierfree<-function(vec)
{
  s<-summary(vec)
  IQ<-s[5]-s[2]
  OUTU<-s[5]+2*IQ
  OUTL<-s[2]-2*IQ
  outliers<-which(vec<=OUTL|vec>=OUTU)
  if (length(outliers)>0) vec<-vec[-outliers]
  vec
}

checkForScarseStates<-function(vec, ielem=1)
{
  if(sum(vec==ielem)==0) 
    stop("Some states occured 0 times\n  Decrease StatesNum or supply more data", call. = FALSE)
  else if(sum(vec==ielem)<=2) {
    warning("Some states occured less than 3 times. They are merged with neighbors", call. = FALSE)
    i<-0
    maxSt<-max(vec)
    while((sum(vec==ielem-i)<=2 && ielem!=1) || (sum(vec==ielem+i)<=2 && ielem==1))
    {
      
      if(ielem==maxSt) {if ((ielem-1-i) %in% vec) vec[vec==ielem]<-ielem-1-i else{}}
      else if(ielem==1){if ((ielem+1+i) %in% vec) vec[vec==ielem]<-ielem+1+i else{}}
      else             {if ((ielem-1-i) %in% vec) vec[vec==ielem]<-ielem-1-i  else{
                            if ((ielem+1+i) %in% vec) vec[vec==ielem]<-ielem+1+i else{}}}
      if (sum(vec==ielem) > 2) break 
      i<-i+1
    if(i>10) stop("Chain is with too many gaps between values.\n Try to reduce the number of states. ", call. = FALSE)
      }
    
    
  } 
  ## end of if(sum(vec==ielem)<=2)
  
  maxSt<-max(vec)
  ielem<-ielem+1  
  if(ielem<=maxSt) checkForScarseStates(vec, ielem)
  else  return (vec)
  }

recodeState<-function(vec)
{
  set<-unique(vec)
  set<-set[order(set)]
  l<-length(set)
  for(i in 1:l) 
  {vec[vec==set[i]]<-i}
  vec  
}

NumVals<-function(vecNext, StatesNum)
  {
  counts<-sapply(1:StatesNum, sumMatch, vec=vecNext)
  }

sumMatch<-function(vec, value)
  {
  sum(match(vec, table=value),  na.rm = TRUE)
  }



TransProbM<-function(vecX, Nstates)
  
{
  sapply(1:Nstates, getTransMatrow, x=vecX, Nst=Nstates)
}
  
getTransMatrow<-function(x, stateNumber, Nst)
  {
  nl<-length(x)  
  cc<-which(x==stateNumber)
  #exclude last instance of Astate i if it's end of chain
  if (cc[length(cc)]==nl) cc<-cc[-length(cc)]
  x[cc]
  transCountFrom<-NumVals(x[cc+1], Nst)
  transCountFrom[transCountFrom[]==0]<-1
  transCountFrom/sum(transCountFrom)
  
}

condprob<-function(probM, stateA, niterStart)
{
  p<-vector(mode="numeric", length=niterStart)
  S<-probM
  for(i in 1:niterStart){
    
    p[i]<-S[stateA,stateA]
    S<-S%*%probM
  }
  p
}

deviations<-function(cp, pA){
  
  depth<-length(cp)-1
  denom<-abs(cp[2]-pA)
  if(denom < 1e-4) denom<-1e-4 else{}
  num<-sapply(2:depth, function(x) abs(cp[x]-pA))
  list(num,denom)
}

optimr1<-function(restr){
  N<-length(restr)
  itern<-1:N
  eval_f0 <- function( x, restr, itern){ return( x )  }
  eval_g0 <- function( x, restr, itern){ return( restr - x^(-itern)) }
  nloptr( x0=c(1.5), eval_f=eval_f0,
                 lb = c(1.01),
                 ub = c(5),
                  eval_g_ineq = eval_g0,
                  opts = list("algorithm"="NLOPT_LN_COBYLA", "xtol_rel"=1.0e-8),
                  restr = restr,
                  itern=itern
                  )

  }

optimr2<-function(restr, ropt){
  N<-length(restr)
  itern<-1:N
  eval_f0 <- function( x, restr, itern){ return( x )  }
  eval_g0 <- function( x, restr, itern){ return( restr - x*ropt^(-itern)) }
  nloptr( x0=c(0.5), eval_f=eval_f0,
          lb = c(0.00001),
          ub = c(5),
          eval_g_ineq = eval_g0,
          opts = list("algorithm"="NLOPT_LN_COBYLA", "xtol_rel"=1.0e-8),
          restr = restr,
          itern=itern
  )
}

Epsilon<-function(r_1,r_2,r_3, M_1,M_2,M_3, pAstate, lcprob){
  depth<-1:lcprob
  rMin23<-ifelse(r_3>=r_2, r_2, r_3)
  k23<-ifelse(r_3>=r_2, r_3/r_2, r_2/r_3)
  rMin13<-ifelse(r_3>=r_1, r_1, r_3)
  k13<-ifelse(r_3>=r_1, r_3/r_1, r_1/r_3)  
  lratio23<-(1-k23^(-depth))/(k23-1)
  lratio13<-(1-k13^(-depth))/(k13-1)
  if(abs(r_3-1)>=0.0002) 
  add1<-2*M_3*r_3^(1-depth)/(r_3-1)
  else 
  add1<- 2*M_3*(1-depth)
  
  add2<-  M_2*M_3*r_3*lratio23*rMin23^(-depth-1)/(r_2-1)
  #add2<-  M_2*M_3*r_3*(r_3^(-depth)-r_2^(-depth))/((r_2-1)*(r_2-r_3))
  
  multypl<- M_1*M_2*M_3/(r_2-r_1)
  add3<-   multypl*r_1*r_3*lratio13*rMin13^(-depth-1)
  add4<-   multypl*r_2*r_3*lratio23*rMin23^(-depth-1)
  add1+add2+add3+add4
  
}

calcSpan<-function(vec, start, Ast, eps2)
{
 lvec<-length(vec)
 
  Astplaces<-which(vec[start:lvec]==Ast)
  m<-length(Astplaces)
  mint<-m-1 
  S<-vector(mode="numeric", length=mint)
  MeanVal<-mean(vec[Astplaces[1]:Astplaces[m]])
# MeanVal<-fitdistr(vec[Astplaces[1]:Astplaces[m]],"normal")$estimate[[1]]
  
  for(i in 1:mint){
    
  S[i]<-sum(vec[(Astplaces[i]+1):Astplaces[i+1]])

   TheorSum<-MeanVal*length((Astplaces[i]+1):Astplaces[i+1])
   S[i]<-S[i]-TheorSum
    
  }

  varSpan<-sum(S^2)
  round(sqrt(varSpan/eps2))
}
  

    

