setMethod(
  f= "plot",
  signature="powerClass",
  definition=function(x,stencil,...){

    moreargs<-list(...)
    validObject(x)

    stopifnot("\nargument stencil accepts values {1,2,3} only"= (!missing(stencil) && stencil %in% 1:3) || missing(stencil))
    names<-x["effectsname"]
    vbetasize<-length(names)
#    default<-par(mar=c(3,3,1,1), cex=0.7, cex.main=0.8, cex.axis=0.8)
    default_par<-par(no.readonly = TRUE)
    on.exit(par(default_par), add = TRUE)
    samplesizeidx<-seq(from=x["Ntotal",2], to=x["Ntotal",3], by=(x["Ntotal",3] - x["Ntotal",2])/10)*x["Ntotal",1]
    samplesize<-rep(samplesizeidx, each=100)



    if (missing(stencil)) {
      plotz(x,vbetasize, samplesize, samplesizeidx)
      plotp(x,vbetasize, samplesize, samplesizeidx)
    }


     else if(stencil==1) {plotz(x,vbetasize, samplesize, samplesizeidx)}
     else if(stencil==2) {plotp(x,vbetasize, samplesize, samplesizeidx)}
     else if(stencil==3) {

       for(i in 1:vbetasize){
       plotzp(x,i, samplesize, samplesizeidx)
       }
     }

        else{}


 #   par(default)
  }

      )

plotz<-function(x,vbetasize, samplesize, samplesizeidx) {
  par(mfrow=c(1,1))

  for(i in 1:vbetasize){
    z<-numeric(0)
    z0.5<-numeric(0)
    z0.1<-numeric(0)
    z0.05<-numeric(0)

    for(j in 1:11) {
      z<-c(z, x[paste0("power",j),i]$z)
      zq<-quantile(x[paste0("power",j),i]$z, probs=c(0.05,0.1,0.5))
      z0.5<-c(z0.5,zq[3])
      z0.1<-c(z0.1,zq[2])
      z0.05<-c(z0.05,zq[1])
    }
      plot(samplesize,z, main=paste0("Effect: ",x["effectsname"][i]),cex=0.2, cex.main=0.5, cex.axis=0.6, col="grey", ylab="", xlab="")
      mtext("sample size", side = 1, line = 2, cex=0.7)
      mtext("z-score", side = 2, line = 2, cex=0.7)
      lines(samplesizeidx, z0.5, lty="solid", col="steelblue3", lwd=1)
      lines(samplesizeidx, z0.1, lty="dashed", col="steelblue2", lwd=1)
      lines(samplesizeidx, z0.05, lty="dotted", col="steelblue1", lwd=1)
      usr <- par("usr")
      lx<-usr[1]+(usr[2]-usr[1])*0.68
      ly<-usr[3]+(usr[4]-usr[3])*0.2
      legend(lx, ly, bty="n", c("median","0.1 quantile", "0.05 quantile"), lty=c("solid", "dashed", "dotted"),
             col=c("steelblue3", "steelblue2", "steelblue1"), lwd=1, cex=0.6)
    }

}

plotp<-function(x, vbetasize, samplesize, samplesizeidx) {
  par(mfrow=c(1,1))

  for(i in 1:vbetasize){
    p<-numeric(0)
    p0.5<-numeric(0)
    p0.1<-numeric(0)
    p0.05<-numeric(0)

    for(j in 1:11) {
      p<-c(p, x[paste0("power",j),i]$power)
      pq<-quantile(x[paste0("power",j),i]$power, probs=c(0.05,0.1,0.5))
      p0.5<-c(p0.5,pq[3])
      p0.1<-c(p0.1,pq[2])
      p0.05<-c(p0.05,pq[1])
    }
      plot(samplesize,p, main=paste0("Effect: ",x["effectsname"][i]), cex.main=0.5, cex=0.2, cex.axis=0.6, col="grey", ylab="", xlab="")
      mtext("sample size", side = 1, line = 2, cex=0.7)
      mtext("power", side = 2, line = 2, cex=0.7)
      lines(samplesizeidx, p0.5, lty="solid", col="indianred4", lwd=1)
      lines(samplesizeidx, p0.1, lty="dashed", col="indianred3", lwd=1)
      lines(samplesizeidx, p0.05, lty="dotted", col="indianred2", lwd=1)
      usr <- par("usr")
      lx<-usr[1]+(usr[2]-usr[1])*0.68
      ly<-usr[3]+(usr[4]-usr[3])*0.2
      legend(lx, ly, bty="n", c("median","0.2 quantile", "0.05 quantile"), lty=c("solid", "dashed", "dotted"),
             col=c("indianred4", "indianred3", "indianred2"), lwd=1, cex=0.6)

  }

  }

plotzp<-function(x, number, samplesize, samplesizeidx) {
  par(mfrow=c(1,2))
  par(oma=c(0,0,1,0)) ##leave upper outer margin for Effects name
  par(mgp=c(2,0,0)) ##place labels tight on axes
  par(mar=c(2,1,1,0))

  i<-number
  z<-numeric(0)
  p<-numeric(0)
  z0.5<-numeric(0)
  z0.1<-numeric(0)
  z0.05<-numeric(0)
  p0.5<-numeric(0)
  p0.1<-numeric(0)
  p0.05<-numeric(0)

  for(j in 1:11) {
    z<-c(z, x[paste0("power",j),i]$z)
    p<-c(p, x[paste0("power",j),i]$power)
    zq<-quantile(x[paste0("power",j),i]$z, probs=c(0.05,0.1,0.5))
    pq<-quantile(x[paste0("power",j),i]$power, probs=c(0.05,0.1,0.5))
    z0.5<-c(z0.5,zq[3])
    z0.1<-c(z0.1,zq[2])
    z0.05<-c(z0.05,zq[1])
    p0.5<-c(p0.5,pq[3])
    p0.1<-c(p0.1,pq[2])
    p0.05<-c(p0.05,pq[1])
  }
  plot(samplesize,z, main="z-score",cex=0.2, cex.main=0.6, cex.axis=0.6, col="grey", ylab="", xlab="", tck=0)
  mtext("sample size", side = 1, line = 1, cex=0.6)
  lines(samplesizeidx, z0.5, lty="solid", col="steelblue3", lwd=1)
  lines(samplesizeidx, z0.1, lty="dashed", col="steelblue2", lwd=1)
  lines(samplesizeidx, z0.05, lty="dotted", col="steelblue1", lwd=1)
  usr <- par("usr")
  lx<-usr[1]+(usr[2]-usr[1])*0.5
  ly<-usr[3]+(usr[4]-usr[3])*0.2
  legend(lx, ly, bty="n", c("median","0.1 quantile", "0.05 quantile"), lty=c("solid", "dashed", "dotted"),
         col=c("steelblue3", "steelblue2", "steelblue1"), lwd=1, cex=0.5)

  plot(samplesize,p, main="power", cex.main=0.6,cex=0.2, cex.axis=0.6, col="grey", ylab="", xlab="", tck=0)
  mtext("sample size", side = 1, line = 1, cex=0.6)
  lines(samplesizeidx, p0.5, lty="solid", col="indianred4", lwd=1)
  lines(samplesizeidx, p0.1, lty="dashed", col="indianred3", lwd=1)
  lines(samplesizeidx, p0.05, lty="dotted", col="indianred2", lwd=1)
  usr <- par("usr")
  lx<-usr[1]+(usr[2]-usr[1])*0.5
  ly<-usr[3]+(usr[4]-usr[3])*0.2
  legend(lx, ly, bty="n", c("median","0.2 quantile", "0.05 quantile"), lty=c("solid", "dashed", "dotted"),
         col=c("indianred4", "indianred3", "indianred2"), lwd=1, cex=0.5)

 mtext(paste0("Effect: ",x["effectsname"][i]), side=3, outer = TRUE, cex=0.6)

 par(mgp=c(3,1,0))
 par(oma=c(0,0,0,0))
 par(mfrow=c(1,1))
 par(mar=c(3,3,1,1))
  }
