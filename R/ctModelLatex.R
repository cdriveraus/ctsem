#' Generate and optionally compile latex equation of subject level ctsem model.
#'
#' @param ctmodel ctsem model object
#' @param matrixnames Logical. If TRUE, includes ctsem matrix names such as DRIFT and DIFFUSION under the matrices.
#' @param textsize Standard latex text sizes -- 
#' tiny scriptsize footnotesize small normalsize large Large LARGE huge Huge. 
#' Useful if output overflows page. 
#' @param filename filename, without suffix, to output .tex and .pdf files too.
#' @param folder Character string specifying folder to save to, defaults to temporary directory, use "./" for working directory.
#' @param tex Save .tex file? Otherwise latex is simply returned within R as a string.
#' @param equationonly Logical. If TRUE, output is only the latex relevant to the equation, not a compileable document.
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
ctModelLatex<- function(ctmodel,matrixnames=TRUE,textsize='normalsize',folder=tempdir(),
  filename=paste0('ctsemTex',as.numeric(Sys.time())),tex=TRUE, equationonly=FALSE, compile=TRUE, open=TRUE){
  
  if(class(ctmodel) == 'ctStanModel') {
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
  
  
  W <- diag(1,ctmodel$n.latent)
  if(continuoustime) diag(W) <- 'u-t'
  
  out <- ifelse(equationonly,"","
\\documentclass[a4paper,landscape]{report}
\\usepackage[margin=1cm]{geometry}
\\usepackage{amsmath} %for multiple line equations
\\usepackage{bm}
\\newcommand{\\vect}[1]{\\boldsymbol{\\mathbf{#1}}}


\\begin{document}
\\thispagestyle{empty}")
  out <- paste0(out, "
\\setcounter{MaxMatrixCols}{200}
 \\begin{",textsize,"}
  \\begin{align*}
  &\\underbrace{",ifelse(continuoustime,"\\mathrm{d}",""),"
    ",bmatrix(matrix(paste0(ctmodel$latentNames)))," 
    \\big{(}t\\big{)}}_{",ifelse(continuoustime,"\\mathrm{d}","")," \\vect{\\eta} (t)}	=  \\left(
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
    & \\qquad \\qquad \\quad \\underbrace{
      ",bmatrix(ctmodel$DIFFUSION),"
    ",ifelse(!matrixnames,"}_{{", "}_{\\underbrace{"),"\\vect{G}}",ifelse(!matrixnames,"}","_\\textrm{DIFFUSION}}"),"
    \\underbrace{",ifelse(continuoustime,"\\mathrm{d}",""),"
      ",bmatrix(matrix(paste0('W_{',1:ctmodel$n.latent,'}')),nottext=TRUE)," 
      (t)}_{",ifelse(continuoustime,"\\mathrm{d}","")," \\vect{W}(t)} \\\\ \\\\
          &",if(continuoustime) paste0("\\underbrace{
            ",bmatrix(matrix(paste0('W_{',1:ctmodel$n.latent,'}')),nottext=TRUE),"  
            (t+u)}_{\\vect{W}(t+u)} - "),"  \\underbrace{",bmatrix(matrix(paste0('W_{',1:ctmodel$n.latent,'}')),nottext=TRUE),"  
            (t)}_{\\vect{W}(t)} \\sim  \\mathrm{N} \\left(
              ",bmatrix(matrix(0,ctmodel$n.latent,1)),", ",bmatrix(W)," \\right) \\\\ \\\\
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
          &\\underbrace{
            ",bmatrix(matrix(paste0('\\epsilon_{',1:ctmodel$n.manifest,'}')))," 
            (t)}_{\\vect{\\epsilon}(t)} \\sim  \\mathrm{N} \\left(
              ",bmatrix(matrix(0,ctmodel$n.manifest,1)),"
              ,
                ",bmatrix(diag(1,ctmodel$n.manifest))," \\right) 
      \\end{align*}
      \\end{",textsize,"}
      ")
  
  if(!equationonly) out <- paste0(out, "\\end{document}")
  
  
  
  if(tex) {
    oldwd <- getwd()
    setwd(dir = folder)
    on.exit(setwd(oldwd))
    write(x = out,file = paste0(filename,'.tex'))
    if(compile){
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
      if(interactive() && open) try(openPDF(paste0(filename,'.pdf')))
    }
    
  }
  return(out)
}


