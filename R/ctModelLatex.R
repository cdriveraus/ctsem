ctModelBuildPopCov <- function(ctm,linearise){ #for latex
  ctm <- T0VARredundancies(ctm)
  pars <- unique(ctm$pars$param[ctm$pars$indvarying])
  d=length(pars)
  m <- matrix(paste0(ifelse(linearise,'','raw'),'PCov_',rep(1:d,d),'_',rep(1:d,each=d)),d,d,dimnames = list(pars,pars))
  m[upper.tri(m)]=t(m)[upper.tri(m)]
  return(m)
}

# getPopEffectsFromFit <- function(x,linearise=TRUE,digits=3){
#   # browser()
#   ms=ctMatsetupFreePars(x$setup$matsetup)
#   e=x$stanfit$transformedparsfull#ctExtract(x)
#   if(x$standata$ntipred > 0){
#     if(linearise) timat <- e$linearTIPREDEFFECT
#     if(!linearise) timat <- e$TIPREDEFFECT
#     dimnames(timat) <- list(iter=1:(dim(timat)[1]),param= ms$parname,
#       TIpred= x$ctstanmodel$TIpredNames)
#     timat <- timat[,apply(x$standata$TIPREDEFFECTsetup,1,function(x) any(x!=0)),,drop=FALSE]
#   } else timat <- diag(0,0)
#   
#   if(!is.null(e$rawpopc) && !is.null(e$popcov)){ 
#     if(!linearise) popcov <- matrix(e$rawpopc[1,4,,],dim(e$rawpopc)[3])
#     if(linearise) {
#       popcov <- matrix(e$popcov[1,,],dim(e$popcov)[3])
#       if(x$standata$intoverpop==1){
#         t0index <- ms$indvarying[ms$param > 0 & ms$row <= x$standata$nlatent & ms$matrix %in% 1 & ms$indvarying > 0]
#         popcov[t0index,t0index] <- e$pop_T0VAR[1,
#           t0index,t0index] #is this correct...?
#       }
#     }
#     dimnames(popcov) <- list(#iter=1:(dim(popcov)[1]),
#       ms$parname[as.logical(ms$indvarying)],
#       ms$parname[as.logical(ms$indvarying)] )
#     
#   }
#   return(list(popcov=popcov,tieffects=timat))
# }

ctModelBuildTIeffects <- function(ctm){ #for latex
  ctm$pars <- ctStanModelCleanctspec(ctm$pars)
  tieffects <- unique(colnames(ctm$pars)[grep('_effect',colnames(ctm$pars),fixed=TRUE)])
  pars <- unique(ctm$pars$param[apply(ctm$pars[,tieffects,drop=FALSE],1,any)])
  timat <- matrix(0,length(pars),length(tieffects),dimnames = list(pars,gsub('_effect','',tieffects)))
  if(length(tieffects)){
    for(p in 1:length(pars)){
      timat[p,] <- unlist(ctm$pars[match(x = pars[p],ctm$pars$param),tieffects,drop=FALSE])
    }
  }
  for(i in 1:nrow(timat)){
    for(j in 1:ncol(timat)){
      if(timat[i,j] !=0) timat[i,j] = paste0(pars[i],'_',gsub('_effect','',tieffects[j],fixed=TRUE))
    }}
  return(timat)
}

ctMatsetupFreePars <- function(m,intoverpop){
  m=m[m$when %in% c(0,-1) & m$param > 0,,drop=FALSE]
  m=m[match(unique(m$param),m$param),,drop=FALSE]
  m = m[order(m$param),,drop=FALSE]
}

ctMatvalueFreePars <- function(ms,mv){
  mv=mv[ms$when %in% c(0,-1) & ms$param > 0,,drop=FALSE]
  ms=ms[ms$when %in% c(0,-1) & ms$param > 0,,drop=FALSE]
  mv=mv[match(unique(ms$param),ms$param),,drop=FALSE]
  mv = mv[order(ms$param),,drop=FALSE]
}

# ctCovTransform <- function(rawpopcov, rawpopmeans, ms, mv){
#   mv <- ctMatvalueFreePars(ms,mv)
#   ms <- ctMatsetupFreePars(ms)
#   rawpopmeans=rawpopmeans[ms$param[ms$indvarying>0]]
#   
#   d=nrow(rawpopcov)
#   n=10000
#   mc <- t(chol(rawpopcov))
#   
#   x <- matrix(rnorm(n*d),n,d)
#   x <- t(apply(x,1,function(y)  mc %*% (y) + rawpopmeans ))
#   tx <- x
#   for(i in 1:nrow(mc)){
#       tx[,i] <- ctsem:::tform(x[,i],ms$transform[i],mv$multiplier[i], mv$meanscale[i], mv$offset[i],mv$inneroffset[i])
#   }
#   return(cov(tx))
# }

texPrep <- function(x){ #replaces certain characters with tex safe versions
  for(i in 1:length(x)){
    x[i]=gsub('_', '\\_',x[i],fixed=TRUE)
    x[i]=gsub('^', '\\textasciicircum',x[i],fixed=TRUE)
  }
  return(x)
}

bmatrix = function(x, digits=NULL,nottext=FALSE, ...) {
  if(!is.null(x)){
    if(!nottext){
      for(i in 1:length(x)){
        if(is.na(suppressWarnings(as.numeric(x[i]))) & #if x[i] cannot be numeric and 
            grepl('\\',x[i],fixed=TRUE) == FALSE) {
          x[i] = texPrep(x[i])
          x[i] = paste0('\\text{',x[i],'}')
        }
      }
    }
    x=as.matrix(x)
    
    out=c()
    for(i in 1:nrow(x)){
      for(j in 1:ncol(x)){
        out=c(out,x[i,j])
        if(j!=ncol(x)) out=c(out,paste0(' & '))
        if(j==ncol(x) & i!=nrow(x)) out=c(out,paste0('\\\\ \n'))
        if(j==ncol(x) & i==nrow(x)) out=c(out, '\n')
      }
    }
    
    out = paste0("\\begin{bmatrix}\n",
      paste0(out,collapse=''),
      "\\end{bmatrix}",collapse='')
  } else out=""
  return(out)
}


#' Generate and optionally compile latex equation of subject level ctsem model.
#'
#' @param x ctsem model object or ctStanFit object.
#' @param matrixnames Logical. If TRUE, includes ctsem matrix names such as DRIFT and DIFFUSION under the matrices.
#' @param digits Precision of decimals for numeric values.
#' @param linearise Logical. Show the linearised normal approximation for subject parameters and 
#' covariate effects, or the raw parameters?
#' @param textsize Standard latex text sizes -- 
#' tiny scriptsize footnotesize small normalsize large Large LARGE huge Huge. 
#' Useful if output overflows page. 
#' @param filename filename, without suffix, to output .tex and .pdf files too.
#' @param folder Character string specifying folder to save to, defaults to temporary directory, use "./" for working directory.
#' @param tex Save .tex file? Otherwise latex is simply returned within R as a string.
#' @param equationonly Logical. If TRUE, output is only the latex relevant to the equation, not a compileable document.
#' @param minimal if TRUE, outputs reduced form version displaying matrix dimensions and equation structure only.
#' @param compile Compile to .pdf? (Depends on \code{tex = TRUE}) 
#' @param open Open after compiling? (Depends on \code{compile = TRUE})
#' @param includeNote Include text describing matrix transformations and subject notation?
#' triangular matrices (which results in a covariance or Cholesky matrix) is shown -- 
#' the latter is a more direct representation of the model, while the former is often simpler to convey.
#'
#' @return character string of latex code. Side effects include saving a .tex, .pdf, and displaying the pdf. 
#' @export
#' @importFrom tools texi2pdf
#'
#' @examples
#' ctmodel <- ctModel(type='stanct', 
#' n.latent=2, n.manifest=1, 
#' manifestNames='sunspots', 
#' latentNames=c('ss_level', 'ss_velocity'),
#' LAMBDA=matrix(c( 1, 'ma1' ), nrow=1, ncol=2),
#' DRIFT=matrix(c(0, 1,   'a21', 'a22'), nrow=2, ncol=2, byrow=TRUE),
#' MANIFESTMEANS=matrix(c('m1'), nrow=1, ncol=1),
#' CINT=matrix(c(0, 0), nrow=2, ncol=1),
#' DIFFUSION=matrix(c(
#'   0, 0,
#'   0, "diffusion"), ncol=2, nrow=2, byrow=TRUE))
#'   
#' l=ctModelLatex(ctmodel,compile=FALSE, open=FALSE)
#' cat(l)
ctModelLatex<- function(x,matrixnames=TRUE,digits=3,linearise=class(x) %in% 'ctStanFit',textsize='normalsize',folder=tempdir(),
  filename=paste0('ctsemTex',as.numeric(Sys.time())),tex=TRUE, equationonly=FALSE, compile=TRUE, open=TRUE, includeNote=TRUE,
  minimal=FALSE){
  #library(ctsem)
  dopopcov <- FALSE
  
  if('ctStanFit' %in% class(x)){
    ms=ctMatsetupFreePars(x$setup$matsetup)
    e=ctExtract(x)
    if(x$standata$ntipred > 0){
      if(linearise) timat <- round(ctCollapse(e$linearTIPREDEFFECT,1,mean),digits)
      if(!linearise) timat <- round(ctCollapse(e$TIPREDEFFECT,1,mean),digits)
      rownames(timat) <- ms$parname
      colnames(timat) <- x$ctstanmodel$TIpredNames
      timat <- timat[apply(x$standata$TIPREDEFFECTsetup,1,function(x) any(x!=0)),,drop=FALSE]
    } else timat <- diag(0,0)
    
    if(!is.null(e$rawpopcov)){ 
      
      if(!linearise) popcov <- round(ctCollapse(e$rawpopcov,1,mean),digits)
      if(linearise) {
        popcov <- stan_constrainsamples(x$stanmodel,x$standata,matrix(x$stanfit$rawest,nrow=1),
          cores=1,pcovn =1000,dokalman=FALSE,savesubjectmatrices = FALSE)$popcov
        popcov <- round(ctCollapse(e$popcov,1,mean),digits=digits)
        if(x$standata$intoverpop==1){
          t0index <- ms$indvarying[ms$param > 0 & ms$row <= x$standata$nlatent & ms$matrix %in% 1 & ms$indvarying > 0]
          popcov[t0index,t0index] <- round(ctCollapse(e$pop_T0cov,1,mean),digits=digits)[
            t0index,t0index] 
        }
        rownames(popcov) <- ms$parname[as.logical(ms$indvarying)]
        colnames(popcov) <- ms$parname[as.logical(ms$indvarying)]
      }
    } else popcov <- diag(0,0)
    
    if(!linearise) popmeans <- round(ctCollapse(e$rawpopmeans,1,mean),digits)[
      as.logical(ms$indvarying + ms$tipred),drop=FALSE]
    if(linearise) {
      popmeans <- round(ctCollapse(e$popmeans,1,mean),digits)[
        as.logical(ms$indvarying + ms$tipred),drop=FALSE]
      if(x$standata$intoverpop==1){
        popmeans[t0index]<- round(ctCollapse(e$pop_T0MEANS,1,mean),digits=digits)[
          t0index,1]
      }
    }
    
    # parmats <- summary(x,residualcov=FALSE,priorcheck=FALSE,digits=digits)
    # parmats <- data.frame(parmats$parmatrices,matrix=rownames(parmats$parmatrices))
    ctmodelmats <- listOfMatrices((x$ctstanmodelbase$pars))
    ctmodel <- x$ctstanmodelbase
    ####################################################################
    for(mi in names(ctmodelmats)){
      mimean <- ctCollapse(e[[paste0('pop_',mi)]],1,mean)
      for(i in 1:nrow(ctmodelmats[[mi]])){
        for(j in 1:ncol(ctmodelmats[[mi]])){
          ctmodelmats[[mi]][i,j] <- round(mimean[i,j],digits)
          # if(ctmodel$pars$matrix[i] %in% parmats$matrix){
          #   try(ctmodel$pars$value[i] <- parmats[parmats$matrix %in% ctmodel$pars$matrix[i] & 
          #       ctmodel$pars$row[i] == parmats$Row & ctmodel$pars$col[i]==parmats$Col,'Mean'],silent=TRUE)
        }
      }
    }
    ctmodel <- c(ctmodel,ctmodelmats)
    class(ctmodel) <- 'ctStanModel'
  } else ctmodel <- x
  
  if('ctStanModel' %in% class(ctmodel)) {
    
    if(!'ctStanFit' %in% class(x)){ #construct pop effects
      popcov <- ctModelBuildPopCov(ctmodel,linearise=linearise)
      if(ctmodel$n.TIpred > 0) timat <- ctModelBuildTIeffects(ctmodel) else timat <- diag(0,0)
      if(!linearise) timat[,] <- paste0('raw_',timat)
      popmeans<-paste0(ifelse(linearise,'','raw_'),unique(c(rownames(popcov),rownames(timat))))
      ctmodel<-T0VARredundancies(ctmodel)
      ctmodel <- c(ctmodel,listOfMatrices(ctmodel$pars)) 
    }
    
    dopopcov <- as.logical(nrow(popcov))
    doti <- as.logical(nrow(timat))
    
    if(doti){
      # both <- rownames(timat) %in% rownames(popcov)
      pars <- unique(c(rownames(popcov),rownames(timat)))
      newpopcov <- matrix(0,length(pars),length(pars),dimnames=list(pars,pars))
      newtimat <- matrix(0,length(pars),ncol(timat),dimnames=list(pars,colnames(timat)))
      newtimat[na.omit(match(rownames(timat),rownames(newtimat))),] <- timat
      newpopcov[na.omit(match(rownames(popcov),rownames(newpopcov))), 
        na.omit(match(rownames(popcov),rownames(newpopcov)))] <- popcov
      popcov <- newpopcov
      timat <- newtimat
    }
    
    dopop <- doti||dopopcov
    
    ### this section replaced t0var fixed values with params, seemed broken...
    # if(!'ctStanFit' %in% class(x)){ #if a model
    #   
    #   t0index <- unique(ctmodel$pars$row[ctmodel$pars$matrix %in% 'T0MEANS' & 
    #       ctmodel$pars$indvarying & is.na(ctmodel$pars$value)]) #which t0means are indvarying
    #   if(length(t0index)){
    #     t0varpopcov <- matrix(
    #       paste0('Pcorsqrt_',t0index,'_',rep(t0index,each=length(t0index))),
    #       nrow=length(t0index),ncol=length(t0index))
    #     t0varpopcov[upper.tri(t0varpopcov)] <- 0
    #     # ms$indvarying[ms$param > 0 & ms$row <= x$standata$nlatent & ms$matrix %in% 1 & ms$indvarying > 0]
    #     ctmodel$T0VAR[t0index,t0index] <- t0varpopcov
    #   }
    # }
    
    continuoustime <- ctmodel$continuoustime
  } else {
    dopop <- FALSE
    if(! 'ctsemInit' %in% class(ctmodel)) stop('not a ctsem model!')
    continuoustime <- TRUE
  }
  
  if(equationonly) compile <- FALSE
  

  
  
  
  
  W <- diag(1,1)
  if(continuoustime) diag(W) <- 't-u'
  
  #out = 'Hello' 
  
  
  out <- ifelse(equationonly,"",paste0("
\\documentclass[a4paper]{article}
\\usepackage{geometry}
\\geometry{paperwidth=\\maxdimen,paperheight=\\maxdimen,margin=1cm}

\\usepackage[fleqn]{amsmath} %for multiple line equations
\\usepackage[active,tightpage,displaymath]{preview}
\\usepackage{bm}
\\newcommand{\\vect}[1]{\\boldsymbol{\\mathbf{#1}}}

\\begin{document}
\\pagenumbering{gobble}
\\begin{",textsize,"}
"))
  
  if (minimal){
    dict = list('A' = 'DRIFT','b'='CINT','M'='TDPREDEFFECT','G'='DIFFUSION','tau'='MANIFESTMEANS')
    
    for (name in names(dict)) {
      chmat = ctmodel[[dict[[name]]]]
      if (!is.numeric(chmat)){
        dict[[name]] = TRUE
      } else {
        #print('Recognized as numeric')
        if (max(abs(chmat)) < 1e-3) {
          dict[[name]] = FALSE
        } else dict[[name]] = TRUE
      }
    }
    
    nu = ctmodel$n.latent
    c = ctmodel$n.manifest
    l = ctmodel$n.TDpred
    
    tablestring <- '\\begin{center}
\\begin{tabular}'#{c|c|c|c} missing
    equationstring <- '\\begin{align*} \n'
    
    tabledim = '{c'
    tablecont1 = '$\\eta(t)$'
    tablecont2 = paste0('$',nu,'$')
    noisestring = ''
    equationcont = 'd\\eta(t) &= '
    
    if (!dict[['A']] & !dict[['b']] & !dict[['M']] ){ #if all of these are not in the equation, we only have noise.
      if (!dict[['G']]) equationcont = paste0(equationcont,'\\mathbf{0}') 
    } else {
      equationcont = paste0(equationcont, '\\left(')
      if (dict[['A']]){
        tabledim = paste0(tabledim,'|c')
        tablecont1 = paste0(tablecont1,'& $\\mathbf{A}$')
        tablecont2 = paste0(tablecont2,'& $',nu,'\\times', nu,'$')
        equationcont = paste0(equationcont,'\\mathbf{A} \\eta(t)')
      }
      if (dict[['b']]){
        tabledim = paste0(tabledim,'|c')
        tablecont1 = paste0(tablecont1,'& $\\mathbf{b}$')
        tablecont2 = paste0(tablecont2,'& $',nu,'$')
        if(dict[['A']]) equationcont = paste0(equationcont,'+')
        equationcont = paste0(equationcont,'\\mathbf{b}')
      }
      if (dict[['M']] && !l==0){
        tabledim = paste0(tabledim,'|c|c')
        tablecont1 = paste0(tablecont1,'& $\\mathbf{M}$ & $\\chi(t)$')
        tablecont2 = paste0(tablecont2,'& $',nu,'\\times', l,'$','& $',l,'$')
        if(dict[['A']]||dict[['b']]) equationcont = paste0(equationcont,'+')
        equationcont = paste0(equationcont,'\\mathbf{M} \\chi(t)')
      }
    }
    
    if (dict[['A']] || dict[['b']] || dict[['M']] ) equationcont = paste0(equationcont,'\\right) dt')
    
    
    if (dict[['G']]){
      tabledim = paste0(tabledim,'|c|c')
      tablecont1 = paste0(tablecont1,'& $\\mathbf{G}$ & $d\\mathbf{W}(t) $')
      tablecont2 = paste0(tablecont2,'& $',nu,'\\times', nu,'$','& $',nu,'$')
      noisestring = paste0(noisestring,'\\mathbf{W}(t+\\Delta t)-\\mathbf{W}(t) &\\sim N(\\mathbf{0},\\mathrm{diag}(\\Delta t)) \\\\','\n' )
      if (dict[['A']] || dict[['b']] || dict[['M']] ) equationcont = paste0(equationcont,'+')
      equationcont = paste0(equationcont,'\\mathbf{G} d\\mathbf{W}(t)')
    }
    
    tabledim = tabledim = paste0(tabledim,'|c|c')
    tablecont1 = paste0(tablecont1,'& $\\mathbf{y}$ & $\\Lambda $')
    tablecont2 = paste0(tablecont2,'& $',c,'$','& $',nu,'\\times', c,'$')
    equationcont = paste0(equationcont,'\\\\','\n', '\\mathbf{y}(t) &= \\Lambda \\eta(t)')
    
    if (dict[['tau']]){
      tabledim = paste0(tabledim,'|c')
      tablecont1 = paste0(tablecont1,'& $\\tau$')
      tablecont2 = paste0(tablecont2,'& $',c,'$')
      equationcont = paste0(equationcont,'+ \\tau')
    }
    
    tabledim = tabledim = paste0(tabledim,'|c|c}')
    tablecont1 = paste0(tablecont1,'& $\\epsilon(t)$ & $\\Theta$ \\\\','\n', '\\hline')
    tablecont2 = paste0(tablecont2,'& $',c,'$','& $',c,'\\times', c,'$')
    noisestring = paste0(noisestring,'\\epsilon(t) &\\sim N(\\mathbf{0},\\Theta) \\\\','\n' )
    equationcont = paste0(equationcont,'+ \\epsilon(t)')
    
    
    tablestring = paste0(tablestring,tabledim,'\n',tablecont1,'\n',tablecont2,'\n','\\end{tabular}','\n','\\end{center}')
    
    equationstring = paste0(equationstring,noisestring,'\\\\',equationcont,'\n','\\end{align*}')
    
    out= paste0(out,tablestring,'\n',equationstring)
    
    
  } else { #end minimal
    showd <- ifelse(continuoustime,"\\mathrm{d}","") #for continuous or discrete system
    
    # if(covMatrices){
    #   if('ctStanFit' %in% class(m)){
    #     cp <- ctStanContinuousPars(m)
    #     ctmodel$T0VAR <- cp$T0COV
    #     ctmodel$DIFFUSION <- cp$DIFFUSIONcov
    #     } else {
    #       ctmodel$T0VAR[upper.tri(ctmodel$T0VAR)] <- t(ctmodel$T0VAR)[upper.tri(ctmodel$T0VAR)]
    #       ctmodel$DIFFUSION[upper.tri(ctmodel$DIFFUSION)] <- t(ctmodel$DIFFUSION)[upper.tri(ctmodel$DIFFUSION)]
    #     }
    # }
          
        
        
    
    out <- paste0(out, "
 \\setcounter{MaxMatrixCols}{200}
 \\begin{flalign*}
  &\\begin{aligned}
  ",if(dopop) paste0("\\parbox{10em}{\\centering{Subject\\linebreak parameter\\linebreak distribution:}}
             &\\underbrace{",bmatrix(matrix(paste0('\\text{',
               texPrep(colnames(popcov)),'}_i')),nottext=TRUE)," 
            }_{\\vect{\\phi}(i)} ",ifelse(linearise,"\\approx","\\sim"),
    ifelse(linearise,"","\\textrm{tform}\\left\\{"),
    " \\mathrm{N} \\left(
              ",bmatrix(popmeans),", ", bmatrix(popcov)," \\right) ",
    if(doti) paste0(" + \\underbrace{",bmatrix(timat),"}_{\\vect{",ifelse(linearise,"\\hat",""),"\\beta}}","
  \\underbrace{
    ",bmatrix(matrix(colnames(timat))),"}_{\\vect{z}}"),
    ifelse(linearise,"","\\right\\}")," \\\\"),
      "\\parbox{10em}{\\centering{Initial\\linebreak latent\\linebreak state:}}
  &\\underbrace{",bmatrix(matrix(paste0(ctmodel$latentNames)))," 
    \\big{(}t_0\\big{)}}_{\\vect{\\eta} (t_0)}	\\sim \\mathrm{N} \\left(
              \\underbrace{
        ",bmatrix(ctmodel$T0MEANS),"
      ",ifelse(!matrixnames,"}_{{", "}_{\\underbrace{"),"\\vect{}}",ifelse(!matrixnames,"}","_\\textrm{T0MEANS}}"),",
      \\underbrace{UcorSDtoCov \\left\\{","
        ",bmatrix(ctmodel$T0VAR),"\\right\\}","
      ",ifelse(!matrixnames,"}_{{", "}_{\\underbrace{"),"\\vect{Q^{*}}_{t0}}",ifelse(!matrixnames,"}","_\\textrm{T0VAR}}"),"
      \\right) \\\\
      \\parbox{10em}{\\centering{Deterministic\\linebreak change:}}
  &\\underbrace{",showd,"
    ",bmatrix(matrix(paste0(ctmodel$latentNames)))," 
    \\big{(}t\\big{)}}_{",showd," \\vect{\\eta} (t)}	=  \\left(
      \\underbrace{
        ",bmatrix(ctmodel$DRIFT),"
      ",ifelse(!matrixnames,"}_{{", "}_{\\underbrace{"),"\\vect{A}}",ifelse(!matrixnames,"}","_\\textrm{DRIFT}}")," \\underbrace{
        ",bmatrix(matrix(paste0(ctmodel$latentNames)))," 
        \\big{(}t\\big{)}
      }_{\\vect{\\eta} (t",ifelse(continuoustime,"","-1"),")}	+ \\underbrace{
        ",bmatrix(ctmodel$CINT),"
      ",ifelse(!matrixnames,"}_{{", "}_{\\underbrace{"),"\\vect{b}}",ifelse(!matrixnames,"}","_\\textrm{CINT}}"),
      if(ctmodel$n.TDpred > 0) paste0( "+ \\underbrace{
        ",bmatrix(ctmodel$TDPREDEFFECT),"
      ",ifelse(!matrixnames,"}_{{", "}_{\\underbrace{"),"\\vect{M}}",ifelse(!matrixnames,"}","_\\textrm{TDPREDEFFECT}}"),"
      \\underbrace{
        ",bmatrix(matrix(ctmodel$TDpredNames))," 
      }_{\\vect{\\chi} (t)}"),
      "\\right) ",ifelse(continuoustime,"\\mathrm{d}t","")," \\quad + \\nonumber \\\\ \\\\
    \\parbox{10em}{\\centering{Random\\linebreak change:}}
    & \\qquad \\qquad \\quad \\underbrace{UcorSDtoChol\\left\\{
      ",bmatrix(ctmodel$DIFFUSION),"\\right\\}
    ",ifelse(!matrixnames,"}_{{", "}_{\\underbrace{"),"\\vect{G}}",ifelse(!matrixnames,"}","_\\textrm{DIFFUSION}}"),"
    \\underbrace{",showd,"
      ",bmatrix(matrix(paste0('W_{',1:ctmodel$n.latent,'}')),nottext=TRUE)," 
      (t)}_{",showd," \\vect{W}(t)} \\\\ \\\\
              \\parbox{10em}{\\centering{Observations:}}
&\\underbrace{
      ",bmatrix(matrix(ctmodel$manifestNames),nottext=FALSE),"  
      (t)}_{\\vect{Y}(t)} = 
        \\underbrace{
          ",bmatrix(ctmodel$LAMBDA)," 
        ",ifelse(!matrixnames,"}_{{", "}_{\\underbrace{"),"\\vect{\\Lambda}}",ifelse(!matrixnames,"}","_\\textrm{LAMBDA}}")," \\underbrace{
          ",bmatrix(matrix(ctmodel$latentNames))," 
          (t)}_{\\vect{\\eta}(t)} +
        \\underbrace{
          ",bmatrix(ctmodel$MANIFESTMEANS)," 
        ",ifelse(!matrixnames,"}_{{", "}_{\\underbrace{"),"\\vect{\\tau}}",ifelse(!matrixnames,"}","_\\textrm{MANIFESTMEANS}}")," + \\nonumber \\\\ \\\\
    \\parbox{10em}{\\centering{Observation\\linebreak noise:}}
    & \\qquad \\qquad \\quad  \\underbrace{
                ",bmatrix(ctmodel$MANIFESTVAR),"  
              ",ifelse(!matrixnames,"}_{{", "}_{\\underbrace{"),"\\vect{\\Theta}}",ifelse(!matrixnames,"}","_\\textrm{MANIFESTVAR}}"),"
              \\underbrace{
          ",bmatrix(matrix(paste0('\\epsilon_{',1:ctmodel$n.manifest,'}')))," 
          (t)}_{\\vect{\\epsilon}(t)} \\\\ \\\\
                \\parbox{10em}{\\centering{System noise\\linebreak distribution per time step:}}
          &",ifelse(continuoustime,'\\Delta ',''),"\\big[W_{j \\in [1,",ctmodel$n.latent,"]}\\big](t",
      ifelse(continuoustime,'-u',''),")   \\sim  \\mathrm{N}(0,",W,") \\quad
              \\parbox{10em}{\\centering{Observation noise\\linebreak distribution:}}
            ",bmatrix(matrix(paste0('\\epsilon_{j \\in [1,',ctmodel$n.latent,']}')))," 
            (t) \\sim  \\mathrm{N}(0,1) \\\\ \\\\
      \\end{aligned} \\\\",
      if(includeNote) paste0("&\\textrm{Note: } UcorSDtoChol\\textrm{ converts lower tri matrix of standard deviations and unconstrained correlations to Cholesky factor,} \\\\
&UcorSDtoCov =\\textrm{ transposed cross product of UcorSDtoChol, to give covariance, See Driver \\& Voelkle (2018) p11.} \\\\",
if(dopop) paste0("&\\textrm{Individual specific notation (subscript i) only shown for subject parameter distribution -- pop. means shown elsewhere.} \\\\
",if(linearise) "&\\textrm{Linearised approximation of subject parameter distribution shown.} \\\\")),
"\\end{flalign*}
      ")
  }
  
  if(!equationonly) out <- paste0(out, 
    "  \\end{",textsize,"}
\\end{document}")
  
  
  if(tex) {
    oldwd <- getwd()
    setwd(dir = folder)
    on.exit(setwd(oldwd))
    write(x = out,file = paste0(filename,'.tex'))
    if(compile){
      hastex <- !Sys.which('pdflatex') %in% ''
      a=try(tools::texi2pdf(file=paste0(filename,'.tex'),quiet=FALSE, clean=TRUE))
      if('try-error' %in% class(a)) {
        
        if(!grepl('SunOS',Sys.info()['sysname']) && requireNamespace('tinytex',quietly=TRUE)){
          a=try(tinytex::pdflatex(file=paste0(filename,'.tex'), clean=TRUE))
          if('try-error' %in% class(a)) 'Error - Perhaps tinytex needs to be installed via: tinytex::install_tinytex()' 
        } else {
          open <- FALSE
          message('Tex compiler not found -- you could install the tinytex package using:\ninstall.packages("tinytex")\ntinytex::install_tinytex()')
        }
      }
      
      if(!'try-error' %in% class(a) && interactive() && open) try(openPDF(paste0(filename,'.pdf')))
    }
    
  }
  return(invisible(out))
}

