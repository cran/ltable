tableToData<-function(tname, numerictype="", orderedtype="") {
  stopifnot("\nfirst argument is not a table, conversion aborted"=is.data.frame(tname))
  
  x<-dim(tname)
  
  if(colnames(tname)[x[2]]=="Total,p") 
    stop("\nThis table contains relative frequencies.\nPlease use type=1 argument in table_f() to build\ntable appropriate for conversion to data")
  stopifnot("\ntable has only one row, conversion aborted"= x[1]>1)
  stopifnot("\ntable is not produced by table_f"=colnames(tname)[x[2]] %in% c("Total,N", "Total, N") )
  stopifnot("numerictype and orderedtype inputs should be character strings of varnames delimited with comma"
            =is.character(numerictype) && is.character(orderedtype))
  lname<-tname[-x[1], -x[2], drop=FALSE]
  indx<-which(grepl(".+:.+", colnames(lname)))
  checknames<-colnames(lname)[indx]
  checked<-sapply(checknames,strsplit, ":")
  lchecked<-length(checked)
  checkstr<-vector(mode = "character", length = lchecked)
  for(i in 1:lchecked) checkstr[i]<-checked[[i]][1]
  stopifnot ("\nIn order to proceed please remove \":\" symbol in all but last variable names" = all(diff(indx)==1) && length(unique(checkstr)) == 1)
  lname<-reshape(lname, varying=indx, direction="long", sep=":")
  x<-dim(lname)
  lname<-lname[,-(x[2])]
  tmpname<-colnames(lname)[x[2]-1]
  colnames(lname)[x[2]-1]<-"Counts"
  colnames(lname)[x[2]-2]<-tmpname
  row.names(lname)<-1:x[1]
  numvar<-unlist(simplify2array(strsplit(numerictype,",")))
  lnumerictype<-length(unique(numvar))
  if(lnumerictype>0) {
    mvar<-match(numvar, colnames(lname))
    mvar<-mvar[!is.na(mvar)]
    stopifnot("\nVariable names in numerictype input don't correspond to that of the table"
              =length(unique(mvar))==lnumerictype)
   tryCatch(
     for(i in mvar) lname[,i]<-as.numeric(lname[,i]),
    warning=function(c){
    c$message<-paste0(c$message," in some of the variables.\nNumeric transformation of variables ", 
                      numerictype, " is aborted.\nTry again with revised variables list")
    message(c)
          })
    }
  else {}
  
  numvar<-unlist(simplify2array(strsplit(orderedtype,",")))
  stopifnot("\nYou can't order variable Counts"=isFALSE("Counts" %in% numvar))
  lnumerictype<-length(unique(numvar))
  if(lnumerictype>0) {
    mvar<-match(numvar, colnames(lname))
    mvar<-mvar[!is.na(mvar)]
    stopifnot("\nVariable names in orderedtype input don't correspond to that of the table"
              =length(unique(mvar))==lnumerictype)
    tryCatch(
      for(i in mvar) lname[,i]<-as.ordered(lname[,i]),
      warning=function(c){
        c$message<-paste0(c$message," in some of the variables.\nOrdering transformation of variables ", 
                          orderedtype, " is aborted.\nTry again with revised variables list")
        message(c)
      })
  }
  else {} 
  
  lname$Counts<-as.numeric(lname$Counts)
  
  return (lname)
  }
  