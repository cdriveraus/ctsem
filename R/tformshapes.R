tformshapes <- function(singletext=FALSE,transform=NA,jacobian=FALSE,driftdiag=FALSE, parname='param'){
  out = c('(param * meanscale * multiplier +inneroffset + offset)',
    '(log1p(exp(param * meanscale+inneroffset)) * multiplier + offset)',
    '(exp(param * meanscale+inneroffset) * multiplier + offset)',
    '(exp(param*meanscale+inneroffset)/(1+exp(param)) * multiplier + offset)',
    '(((param*meanscale+inneroffset)^3)*multiplier + offset)')
  
  out=gsub('param',parname,out,fixed=TRUE)
  
  
  if(driftdiag && jacobian) out = paste0(out,' * param')
  out = sapply(out,Simplify)
  names(out)=paste0('fn',1:length(out))
  if(jacobian) out = jacobianSymb(out,variables='param')
  if(!is.na(transform)) out = out[transform+1] else transform = 0:4
  if(!singletext)   out = paste0('if(transform==',transform+ifelse(jacobian,ifelse(driftdiag,60,50),0),') out = ',out,';\n',collapse='')
  out=gsub('  ','',out,fixed=TRUE)
  return(out)
}

tform <- function(param, transform, multiplier, meanscale, offset, inneroffset, extratforms='',singletext=FALSE,jacobian=FALSE,driftdiag=FALSE){
  
  if(!is.na(suppressWarnings(as.integer(transform)))) {
    out <- tformshapes(singletext=singletext,transform=as.integer(transform),jacobian=jacobian)
    if(!singletext) paste0(out,extratforms)
    if(singletext) {
      for(i in c('param','multiplier', 'meanscale',  'inneroffset','offset')){
        irep = get(i)
        out <- gsub(pattern = i,replacement = irep,out)
      }
    }
  }
  if(jacobian) transform <- transform + ifelse(driftdiag,60,50)
  if(is.na(suppressWarnings(as.integer(transform)))) out <- transform
  if(!singletext) out <- eval(parse(text=out))
  return(out)
}

Jtformshapes <- function(){
  fn=sapply(tformshapes(singletext = TRUE),function(x) Simplify(x))
names(fn)=paste0('fn',1:length(fn))
jacobianSymb(fn,variables = c('param'))
}
