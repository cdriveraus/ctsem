listOfMatrices <- function(df,matnames=NA){
  mlist <- list()
  for(mi in unique(df$matrix)){
    sel <- df$matrix %in% mi
    rc <- cbind(df$row[sel], df$col[sel])
    value <- df$value[sel]
    fixed <- !is.na(value) #a fixed cell holds its value, a free one its parameter
    m <- matrix(NA, max(rc[,1]), max(rc[,2]))
    # Values, then parameters: a matrix with a free cell is character, its
    # values as as.character() writes them.
    if(any(fixed)) m[rc[fixed,,drop=FALSE]] <- value[fixed]
    if(any(!fixed)) m[rc[!fixed,,drop=FALSE]] <- df$param[sel][!fixed]
    mlist[[mi]] <- m
  }

  if(!is.na(matnames[1])){
    for(mi in matnames){
      if(is.null(mlist[[mi]])) mlist[[mi]] <- matrix(NA,0,0)
    }
  }
  return(mlist)
}
