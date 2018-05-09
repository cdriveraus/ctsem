tformshapes <- function(){
paste0('  if(transform==0) out = param * meanscale * multiplier + offset;
  if(transform==1) out = log(1+(exp(param * meanscale))) * multiplier + offset ;
  if(transform==2) out = exp(param * meanscale) * multiplier + offset;
  if(transform==3) out = inv_logit(param*meanscale) * multiplier + offset;
')
  }

tform <- function(param, transform, multiplier, meanscale, offset){
  if(!is.na(suppressWarnings(as.integer(transform)))) eval(parse(text=paste0(tformshapes())))
  if(is.na(suppressWarnings(as.integer(transform)))) out <- eval(parse(text=transform))  
  return(out)
}

