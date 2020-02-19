ctModelBuildPopCov <- function(ctm){ #for latex
  ctm <- T0VARredundancies(ctm)
  pars <- ctm$pars$param[ctm$pars$indvarying]
  d=length(pars)
  m <- matrix(paste0('PCov_',rep(1:d,d),'_',rep(1:d,each=d)),d,d,dimnames = list(pars,pars))
  m[upper.tri(m)]=t(m[lower.tri(m)])
  return(m)
}

ctModelBuildTIeffects <- function(ctm){ #for latex
  ctm$pars <- ctStanModelCleanctspec(ctm$pars)
  tieffects <- colnames(ctm$pars)[grep('_effect',colnames(ctm$pars),fixed=TRUE)]
  pars <- ctm$pars$param[apply(ctm$pars[,tieffects,drop=FALSE],1,any)]
  timat <- matrix(0,length(pars),length(tieffects),dimnames = list(pars,gsub('_effect','',tieffects)))
  if(length(tieffects)){
    for(p in 1:length(pars)){
      timat[p,] <- unlist(ctm$pars[match(x = pars[p],ctm$pars$param),tieffects,drop=FALSE])
    }
  }
  for(i in 1:nrow(timat)){
    for(j in 1:ncol(timat)){
      if(timat[i,j] !=0) timat[i,j] = paste0('tip_',pars[i],'_',gsub('_effect','',tieffects[j],fixed=TRUE))
    }}
  return(timat)
}

ctMatsetupFreePars <- function(m){
  m=m[m$when==0 & m$param > 0,]
  m=m[match(unique(m$param),m$param),]
}



#' Generate and optionally compile latex equation of subject level ctsem model.
#'
#' @param x ctsem model object or ctStanFit object.
#' @param matrixnames Logical. If TRUE, includes ctsem matrix names such as DRIFT and DIFFUSION under the matrices.
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
ctModelLatex<- function(x,matrixnames=TRUE,textsize='normalsize',folder=tempdir(),
  filename=paste0('ctsemTex',as.numeric(Sys.time())),tex=TRUE, equationonly=FALSE, compile=TRUE, open=TRUE,
  minimal=FALSE){
  #library(ctsem)
  
  dopopcov <- FALSE
  
  if('ctStanFit' %in% class(x)){
    ms=ctMatsetupFreePars(x$setup$matsetup)
    e=ctExtract(x)
    
    if(x$standata$ntipred > 0){
    timat <- ctCollapse(e$TIPREDEFFECT,1,mean)
    rownames(timat) <- ms$parname
    colnames(timat) <- x$ctstanmodel$TIpredNames
    timat <- timat[apply(timat,1,function(x) any(x!=0)),,drop=FALSE]
    } else timat <- diag(0,0)
    
    if(!is.null(e$rawpopcov)){
      browser()
    } else popCov <- diag(0,0)
    
    parmats <- summary(x,residualcov=FALSE,priorcheck=FALSE)
    parmats <- data.frame(parmats$parmatrices,matrix=rownames(parmats$parmatrices))
    ctmodel <- x$ctstanmodelbase
    for(i in 1:nrow(ctmodel$pars)){
      if(is.na(ctmodel$pars$value[i])){
        ctmodel$pars$value[i] <- parmats[parmats$matrix %in% ctmodel$pars$matrix[i] & 
            ctmodel$pars$row[i] == parmats$Row & ctmodel$pars$col[i]==parmats$Col,'Mean']
      }
    }
    
    
  } else ctmodel <- x
  
  if('ctStanModel' %in% class(ctmodel)) {
    
    if(!'ctStanFit' %in% class(x)){ #construct pop effects
    popCov <- ctModelBuildPopCov(ctmodel)
    if(ctmodel$n.TIpred > 0) timat <- ctModelBuildTIeffects(ctmodel) else timat <- diag(0,0)
    ctmodel<-T0VARredundancies(ctmodel)
    }
    
    dopopcov <- as.logical(nrow(popCov))
    doti <- as.logical(nrow(timat))
    
    if(doti){
      # 
      # both <- rownames(timat) %in% rownames(popCov)
      # browser()
      pars <- unique(c(rownames(popCov),rownames(timat)))
      newpopcov <- matrix(0,length(pars),length(pars),dimnames=list(pars,pars))
      newtimat <- matrix(0,length(pars),ncol(timat),dimnames=list(pars,colnames(timat)))
      newtimat[na.omit(match(rownames(timat),rownames(newtimat))),] <- timat
      newpopcov[na.omit(match(rownames(newpopcov),rownames(popCov))), 
        na.omit(match(rownames(newpopcov),rownames(popCov)))] <- popCov
      popCov <- newpopcov
      timat <- newtimat
    }
    
    dopop <- doti||dopopcov
    ctmodel <- c(ctmodel,listOfMatrices(ctmodel$pars)) 
    continuoustime <- ctmodel$continuoustime
  } else {
    if(class(ctmodel) != 'ctsemInit') stop('not a ctsem model!')
    continuoustime <- TRUE
  }
  
  if(equationonly) compile <- FALSE
  
  bmatrix = function(x, digits=NULL,nottext=FALSE, ...) {
    if(!is.null(x)){
      if(!nottext){
        for(i in 1:length(x)){
          if(is.na(suppressWarnings(as.numeric(x[i]))) & #if x[i] cannot be numeric and 
              grepl('\\',x[i],fixed=TRUE) == FALSE) {
            x[i]=gsub('_', '\\_',x[i],fixed=TRUE)
            x[i]=gsub('^', '\\textasciicircum',x[i],fixed=TRUE)
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
  
  
  W <- diag(1,1)
  if(continuoustime) diag(W) <- 't-u'

#out = 'Hello' 

  
out <- ifelse(equationonly,"","
\\documentclass[a4paper]{article}
\\usepackage{geometry}
\\geometry{paperwidth=\\maxdimen,paperheight=\\maxdimen,margin=1cm}

\\usepackage[fleqn]{amsmath} %for multiple line equations
\\usepackage[active,tightpage,displaymath]{preview}
\\usepackage{bm}
\\newcommand{\\vect}[1]{\\boldsymbol{\\mathbf{#1}}}

\\begin{document}
\\pagenumbering{gobble}")

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
} else {
  showd <- ifelse(continuoustime,"\\mathrm{d}","") #for continuous or discrete system
out <- paste0(out, "
 \\begin{",textsize,"}
 \\setcounter{MaxMatrixCols}{200}
  \\begin{align*}
  ",if(dopop) paste0("\\textrm{Subject specific parameters: }
             &\\underbrace{",bmatrix(matrix(paste0('\\text{',
             gsub('_','\\_',colnames(popCov),fixed=TRUE),'}_i')),nottext=T)," 
            }_{\\vect{\\phi}(i)} \\sim  \\textrm{tform}\\left\\{ \\mathrm{N} \\left(
              ",bmatrix(paste0(colnames(popCov))),"
              ,
                ",bmatrix(popCov)," \\right) ",
      if(doti) paste0(" + \\underbrace{",bmatrix(timat),"}_{\\vect{\\beta}}","
  \\underbrace{
    ",bmatrix(matrix(colnames(timat))),"}_{\\vect{z}}"),
    " \\right\\} \\\\"),
  "\\textrm{Initial latent states: }
  &\\underbrace{",bmatrix(matrix(paste0(ctmodel$latentNames)))," 
    \\big{(}t_0\\big{)}}_{\\vect{\\eta} (t_0)}	\\sim \\mathrm{N} \\left(
              \\underbrace{
        ",bmatrix(ctmodel$T0MEANS),"
      ",ifelse(!matrixnames,"}_{{", "}_{\\underbrace{"),"\\vect{}}",ifelse(!matrixnames,"}","_\\textrm{T0MEANS}}"),",
      \\underbrace{covsdcor \\left\\{
        ",bmatrix(ctmodel$T0VAR),"\\right\\}
      ",ifelse(!matrixnames,"}_{{", "}_{\\underbrace{"),"\\vect{Q^{*}}_{t0}}",ifelse(!matrixnames,"}","_\\textrm{T0VAR}}"),"
      \\right) \\\\
      \\textrm{Deterministic change:}
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
        ",bmatrix(matrix(paste0('\\chi_{',1:ncol(ctmodel$TDPREDEFFECT),'}')))," 
      }_{\\vect{\\chi} (t)}"),
    "\\right) ",ifelse(continuoustime,"\\mathrm{d}t","")," \\quad + \\nonumber \\\\ \\\\
    \\textrm{Random change: }
    & \\qquad \\qquad \\quad \\underbrace{cholsdcor\\left\\{
      ",bmatrix(ctmodel$DIFFUSION),"\\right\\}
    ",ifelse(!matrixnames,"}_{{", "}_{\\underbrace{"),"\\vect{G}}",ifelse(!matrixnames,"}","_\\textrm{DIFFUSION}}"),"
    \\underbrace{",showd,"
      ",bmatrix(matrix(paste0('W_{',1:ctmodel$n.latent,'}')),nottext=TRUE)," 
      (t)}_{",showd," \\vect{W}(t)} \\\\ \\\\
              \\textrm{Observations: }
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
        ",ifelse(!matrixnames,"}_{{", "}_{\\underbrace{"),"\\vect{\\tau}}",ifelse(!matrixnames,"}","_\\textrm{MANIFESTMEANS}}")," + 
              \\underbrace{
                ",bmatrix(ctmodel$MANIFESTVAR),"  
              ",ifelse(!matrixnames,"}_{{", "}_{\\underbrace{"),"\\vect{\\Theta}}",ifelse(!matrixnames,"}","_\\textrm{MANIFESTVAR}}"),"
              \\underbrace{
          ",bmatrix(matrix(paste0('\\epsilon_{',1:ctmodel$n.manifest,'}')))," 
          (t)}_{\\vect{\\epsilon}(t)} \\\\ \\\\
                \\textrm{Latent noise per time step : }
          &",ifelse(continuoustime,'\\Delta ',''),"\\big[W_{j \\in [1,",ctmodel$n.latent,"]}\\big](t",
          ifelse(continuoustime,'-u',''),")   \\sim  \\mathrm{N}(0,",W,") \\quad
              \\textrm{Observation noise: }
            ",bmatrix(matrix(paste0('\\epsilon_{j \\in [1,',ctmodel$n.latent,']}')))," 
            (t) \\sim  \\mathrm{N}(0,1) \\\\ \\\\
            &\\textrm{Indivividual specific notation (subscript i) not shown for system dynamics and observation model.} \\\\
&cholsdcor\\textrm{ converts lower tri matrix of std dev and unconstrained correlation to Cholesky factor covariance.} \\\\
&covsdcor =\\textrm{ transposed cross product of cholsdcor, to give covariance.} \\\\
&\\textrm{See Driver \\& Voelkle (2018) p11.}
      \\end{align*}
      \\end{",textsize,"}
      ")
}
  
  if(!equationonly) out <- paste0(out, "\\end{document}")
  
  
  if(tex) {
    oldwd <- getwd()
    setwd(dir = folder)
    on.exit(setwd(oldwd))
    write(x = out,file = paste0(filename,'.tex'))
    if(compile){
      if(requireNamespace('tinytex',quietly=TRUE)){
        tinytex::pdflatex(file=paste0(filename,'.tex'), clean=TRUE)
        
      } else{
      hastex <- !Sys.which('pdflatex') %in% ''
      a=try(tools::texi2pdf(file=paste0(filename,'.tex'),quiet=FALSE, clean=TRUE))
      if('try-error' %in% class(a) && !hastex) {
        open <- FALSE
        message('Tex compiler not found -- you could install the tinytex package using:\ninstall.packages("tinytex")\ntinytex::install_tinytex()')
        # dotiny <- readline('Y/N?')
        # if(dotiny %in% c('Y','y')){
        #   utils::install.packages('tinytex')
        #   if(requireNamespace('tinytex',quietly=TRUE)) tinytex::install_tinytex()
        # }
      }
      } #end else
      
      if(interactive() && open) try(openPDF(paste0(filename,'.pdf')))
    }
    
  }
  return(invisible(out))
}

