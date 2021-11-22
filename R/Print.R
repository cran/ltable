setMethod(
  f= "print", 
  signature="powerClass",
  definition=function(x, choice="power",...){
    moreargs<-list(...)
    validObject(x)
    stopifnot("\nargument choice accepts values {\"power\", \"model\"} only"= (!missing(choice) && choice %in% c("power", "model")) || missing(choice))
    names<-x["effectsname"]
    vbetasize<-length(names)
    vnames<-x["varnames"]    
    vparsize<-length(vnames)
    samplesize<-seq(from=x["Ntotal",2], to=x["Ntotal",3], by=(x["Ntotal",3] - x["Ntotal",2])/10)*x["Ntotal",1]
 
    
    if (choice=="power") {
    
      for(i in 1:vbetasize) {
       
    cat(paste0("\n\nEffect:  ", names[i],"\n\n"))
  
        z0.5<-numeric(0)
        z0.4<-numeric(0)
        z0.3<-numeric(0)
        z0.2<-numeric(0)
        z0.1<-numeric(0)
        z0.05<-numeric(0)
        z0.025<-numeric(0)
        
        p0.5<-numeric(0)
        p0.4<-numeric(0)
        p0.3<-numeric(0)
        p0.2<-numeric(0)
        p0.1<-numeric(0)
        p0.05<-numeric(0)
        p0.025<-numeric(0)
        
        cat("Test statistic Z:                    Quantiles\n")
        cat("Sample size:",  "Q0.025", "Q0.05", "Q0.1", "Q0.2", "Q0.3", "Q0.4", "Q0.5", "\n", sep="\t")
        
        for(j in 1:11) {
          zq<-quantile(x[paste0("power",j),i]$z, probs=c(0.025,0.05,0.1,0.2,0.3,0.4,0.5))
          z0.025<-c(z0.025,zq[1])
          z0.05<-c(z0.05,zq[2])
          z0.1<-c(z0.1,zq[3])
          z0.2<-c(z0.2,zq[4])
          z0.3<-c(z0.3,zq[5])
          z0.4<-c(z0.4,zq[6])
          z0.5<-c(z0.5,zq[7])
          cat(sprintf("%s%*.3f%*.3f%*.3f%*.3f%*.3f%*.3f%*.3f", round(samplesize[j]), 18, z0.025[j], 8, z0.05[j], 7, z0.1[j], 8, z0.2[j], 8, z0.3[j], 8, z0.4[j], 8, z0.5[j]), "\n")
          
      }
        
      
      
      cat("\nPower:                               Quantiles\n")
      
        
      cat("Sample size:",  "Q0.025", "Q0.05", "Q0.1", "Q0.2", "Q0.3", "Q0.4", "Q0.5", "\n", sep="\t")
        
        for(j in 1:11) {
          pq<-quantile(x[paste0("power",j),i]$power, probs=c(0.025,0.05,0.1,0.2,0.3,0.4,0.5))
          p0.025<-c(p0.025,pq[1])
          p0.05<-c(p0.05,pq[2])
          p0.1<-c(p0.1,pq[3])
          p0.2<-c(p0.2,pq[4])
          p0.3<-c(p0.3,pq[5])
          p0.4<-c(p0.4,pq[6])
          p0.5<-c(p0.5,pq[7])
          cat(sprintf("%s%*.2f%*.2f%*.2f%*.2f%*.2f%*.2f%*.2f", round(samplesize[j]), 18, p0.025[j], 8, p0.05[j], 7, p0.1[j], 8, p0.2[j], 8, p0.3[j], 8, p0.4[j], 8, p0.5[j]), "\n")
        }
        
        
      }
    
    } else if (choice=="model"){
      cat("\nCall:\n")
      print (x["cal"])
      
      cat("\n\nCoefficients:\n")
 
      l<-40L
      cat(sprintf("%*s %*s %*s %*s" , l+2, "Estimate", 8, "Std.Error", 8, "|z-score|", 11,"Pr(>|z|)\n")) 
    
      for(j in 1:vparsize) {
        z_value<-abs(x["estim",1]$betas[j]/x["estim",1]$errors[j])  
        
        if (nchar(vnames[j]) >= l-5) {
        spl<-splitname(vnames[j],l-5)
        cat(sprintf("%s%*.2e%*.2e%*.2e%*.1e",spl$sname, l+2-nchar(spl$sname), x["estim",1]$betas[j], 10, x["estim",1]$errors[j],
              10, z_value, 10, pnorm(z_value, lower.tail = FALSE)*2),"\n")
        cat(spl$adding,"\n")
      
        } else{
          
        cat(sprintf("%s%*.2e%*.2e%*.2e%*.1e",vnames[j], l+2-nchar(vnames[j]), x["estim",1]$betas[j], 10, x["estim",1]$errors[j],
                        10, z_value, 10, pnorm(z_value, lower.tail = FALSE)*2),"\n")
          
    }
      }
      cat("\n\nModel fit:\n")
      cat("Weighted nonlinear least-squares fitting\n")
      cat("Method|Solver: Levenberg-Marquardt\n")
      cat("Algorithm: trust region\n")
      cat("initial |f(x)| = ", x["estim",1]$`initial |f(x)|`, "\n")
      cat("final |f(x)| = ", x["estim",1]$`final |f(x)|`, "\n")
      cat("Jacobian reciprocal condition number = ", x["estim",1]$`Jacobian reciprocal condition number`, "\n")
      cat("number of iterations = ", x["estim",1]$`number of iterations`, "\n")
      cat("reason for stopping: ", x["estim",1]$`reason for stopping`, "\n")
      cat("chisq/dof = ", x["estim",1]$`chisq/dof`, "\n")
      cat("status: ", x["estim",1]$status, "\n")
   
      message("\nWarning:\nJacobian reciprocal condition number measures the inverse sensitivity of the solution to small perturbations in the input data. It tends to zero as J tends to singularity indicating solution instability.")
      message("\nThe value of ch-squared per degree of freedom (chisq/dof) approximately 1 indicates a good fit.")
      message("If chisq/dof >> 1  the error estimates obtained from the covariance matrix will be too small and should be multiplied by square root of chisq/dof .")
      message("Poor fit will result from the use of an inappropriate model, and the scaled error estimates may then be outside the range of validity for Gaussian errors.")
      message("BEWARE:\nPoor fit jeopardizes the validity of power analysis.")
    }
      else{}
  }
)
 
splitname<-function(sname, l){
cnt<-0
  while (nchar(sname) >= l) {
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