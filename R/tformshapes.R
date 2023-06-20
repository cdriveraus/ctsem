tformshapes <- function(singletext=FALSE,transform=NA,jacobian=FALSE,driftdiag=FALSE, parname='param',stan=FALSE){
 out = c('param',
    '(log1p_exp(param))',
    '(exp(param))',
    '(1/(1+exp(-param)))',
    '((param)^3)',
    'log1p(param)', #why is this here? results in NA's / warnings. 
    'meanscale',
    '1/(1+exp(-param))',
    'exp(param)',
    '1/(1+exp(-param))-(exp(param)^2)/(1+exp(param))^2',
    '3*param^2',
    '1/(1+param)')
 
 tfvec=c(0:5,50:55)
  
  out=gsub('param',parname,out,fixed=TRUE)
  
  # if(driftdiag && jacobian) out = paste0(out,' * param')
  # out = sapply(out,Simplify)
  # names(out)=paste0('fn',1:length(out))
  # if(jacobian) out = jacobianSymb(out,variables='param')
  
  if(!is.na(transform)&&transform!=0) out = out[tfvec == transform] #ifelse(jacobian,0,1):(length(out)-ifelse(jacobian,1,0))
  if(!singletext) {
    out = paste0('if(transform==', tfvec,') param = ',out,';\n',collapse='')
  
  if(!stan) out <- paste0('param = parin * meanscale + inneroffset; \n ',out,'
  param=param*multiplier;
    if(transform < 49) param = param+offset;')
  if(stan) out <- paste0('if(meanscale!=1.0) param *= meanscale; 
  if(inneroffset != 0.0) param += inneroffset; \n',out,'
  if(multiplier != 1.0) param *=multiplier;
  if(transform < 49 && offset != 0.0) param+=offset;')
  }
  if(singletext) out <- paste0('offset + multiplier*',gsub('param','(param*meanscale+inneroffset)',out))
  
  out=gsub('  ','',out,fixed=TRUE)
  return(out)
}

tform <- function(parin, transform, multiplier, meanscale, offset, inneroffset, extratforms='',singletext=FALSE,jacobian=FALSE,driftdiag=FALSE){
  param=parin
  if(!is.na(suppressWarnings(as.integer(transform)))) {
    out <- tformshapes(singletext=singletext,transform=as.integer(transform))#,jacobian=jacobian)
    if(!singletext) paste0(out,extratforms)
    if(singletext) {
      for(i in c('param','multiplier', 'meanscale',  'inneroffset','offset')){
        irep = get(i)
        out <- gsub(pattern = i,replacement = irep,out)
      }
    }
  }
  # if(jacobian) transform <- transform + ifelse(driftdiag,60,50)
  if(is.na(suppressWarnings(as.integer(transform)))) out <- transform
  if(!singletext) out <- eval(parse(text=out))
  return(out)
}

# Jtformshapes <- function(){
#   fn=sapply(tformshapes(singletext = TRUE),function(x) Simplify(x))
#   names(fn)=paste0('fn',1:length(fn))
#   jacobianSymb(fn,variables = c('param'))
# }
