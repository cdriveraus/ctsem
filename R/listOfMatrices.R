listOfMatrices <- function(df,matnames=NA){
  mats <- unique(df$matrix)
  mlist <- list()
  for(mi in mats){
    mlist[[mi]] <- matrix(NA, max(df$row[df$matrix %in% mi]), max(df$col[df$matrix %in% mi]))
    for(ci in 1:max(df$col[df$matrix %in% mi])){
      for(ri in 1:max(df$row[df$matrix %in% mi])){
        if(is.na(df$value[df$matrix %in% mi & df$row == ri & df$col==ci])) {
          mlist[[mi]][ri,ci] <- df$param[df$matrix %in% mi & df$row == ri & df$col==ci]
        } else {
          mlist[[mi]][ri,ci] <- df$value[df$matrix %in% mi & df$row == ri & df$col==ci]
        }
      }
    }
  }
  
  if(!is.na(matnames[1])){
    for(mi in matnames){
      if(is.null(mlist[[mi]])) mlist[[mi]] <- matrix(NA,0,0)
    }
  }
  return(mlist)
}
