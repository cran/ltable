table_f<-function (data, datavars, type = 1, digits = 2, extended = FALSE, 
          MV = FALSE, cb = FALSE) 
{
  flds <- names(data)
  svars <- unlist(strsplit(datavars, ","))
  pointer <- as.vector(unlist(sapply(svars, function(x) which(flds == x))))
  lengthPointer <- length(pointer)
  switch(type, ifelse(MV == TRUE, t <- ftable(data[, pointer], 
        exclude = NaN), t <- ftable(data[, pointer])), t <- round(prop.table(ftable(data[,
        pointer]), 1), digits = digits), t <- round(prop.table(ftable(data[, 
        pointer]), 2), digits = digits), t <- round(sweep(ftable(data[, 
        pointer]), 1, sum(as.matrix(ftable(data[, pointer]))), 
        "/"), digits = digits)) 
  t[is.na(t)] <- 0
  ifelse(lengthPointer > 1, tm <- addmargins(t, margin = seq_along(dim(t)), 
  FUN = sum, quiet = TRUE), tm <- addmargins(t, margin = 2, 
  FUN = sum, quiet = TRUE))
  dtm <- dim(tm)
  if (type == 2){ 
    tm[dtm[1], ] <- c(round(apply(t(tm[dtm[1], -dtm[2]]), 
    1, function(x) x/sum(t(x))), digits = digits), 1)
  } else if (type == 3) { 
      if (lengthPointer == 1) 
      stop("Please select at least two variables for type=3")
      else {
    tsub <- addmargins(ftable(data[, pointer]), margin = 2, 
                       FUN = sum, quiet = TRUE)
    tm[, dim(tm)[2]] <- c(round(apply(t(tsub[, dim(tsub)[2]]), 
                        1, function(x) x/sum(t(x))), digits = digits), 1)
      }
  }
  str <- ifelse(type %in% 2:4, "Total, p", "Total, N")
  if (lengthPointer > 1) {
    v = as.matrix(unlist(lapply(rownames(as.matrix(t)), function(x) strsplit(x, "_"))))
    dim(v) <- c(lengthPointer - 1, dim(t)[1])
    xx <- rbind(t(v), rep(str, lengthPointer - 1))
  }
  ll <- levels(as.factor(data[, pointer[lengthPointer]]))
  suppressWarnings(ifelse(lengthPointer == 1, newdata <- data.frame(tm), 
                          newdata <- data.frame(cbind(xx, tm))))
  suppressWarnings(if (lengthPointer == 1) 
    names(newdata) <- c(paste(svars[length(svars)], ":", ll, sep = ""), str)
    else names(newdata) <- c(svars[-length(svars)], paste(svars[length(svars)], ":", ll, sep = ""), str))
  if (MV == TRUE & lengthPointer > 1) {
    emptyrec = sapply(newdata[, 1:lengthPointer - 1], function(x) which(x == 
                      "NA" & newdata[, dim(newdata)[2]] == 0))
    delrec = unique(as.array(unlist(emptyrec)))
    suppressWarnings(if (dim(delrec) > 0) 
      newdata <- newdata[-delrec, ])
  }
  if (MV == TRUE & type == 1 & is.na(names(newdata)[dim(newdata)[2]])) {
    names(newdata)[dim(newdata)[2] - 1] <- "NA"
    names(newdata)[dim(newdata)[2]] <- "Total, N"
  }
  if (extended == TRUE & type %in% 2:4) {
    t <- addmargins(ftable(data[, pointer]), margin = seq_along(dim(ftable(data[, pointer]))), FUN = sum, quiet = TRUE)
    N <- sum(ftable(data[, pointer]))
    rb <- rbind(as.matrix(newdata), c(rep("Total, N", lengthPointer -1), t[dim(t)[1], ]))
    suppressWarnings(newdata <- as.data.frame(cbind(rb, t[, dim(t)[2]])))
    newdata[dim(newdata)[1], dim(newdata)[2]] <- N
    names(newdata) <- c(svars[-length(svars)], paste(svars[length(svars)], 
                                ":", ll, sep = ""), "Total, p", "Total, N")
  }
  
  if (extended==TRUE & type %in% 2:4) z<-dim(newdata)[1]-2 else z<-dim(newdata)[1]-1
  suppressWarnings(row.names(newdata)[1:z]<-1:z)
  
  if (cb == TRUE) {
    message("\nWarning:\nYou are about to use the clipboard. All previous data in clipboard will be lost. Do You want to proceed?")
    message("\nType 1 to proceed, Enter to refuse:\n")
    
    suppressMessages(i<-scan(n=1, what="integer"))
    
  if(length(i)>0 && i=="1") {clipr::write_clip(newdata, object_type = "table")} else {message("No clipping was made")} 
  
  }
    cat("\n")
  newdata
}
