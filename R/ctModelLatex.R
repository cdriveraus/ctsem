#' Generate and optionally compile latex equation of subject level ctsem model.
#'
#' @param ctmodel ctsem model object
#' @param textsize Standard latex text sizes -- 
#' tiny scriptsize footnotesize small normalsize large Large LARGE huge Huge. 
#' Useful if output overflows page. 
#' @param filename filename, without suffix, to output .tex and .pdf files too.
#' @param folder Character string specifying folder to save to, defaults to working directory.
#' @param tex Save .tex file? Otherwise latex is simply returned within R as a string.
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
#' l=ctModelLatex(ctmodel,folder=tempdir(),filename="ctsemTex")
#' cat(l)
ctModelLatex<- function(ctmodel,textsize='normalsize',folder='./',filename='ctsemTex',tex=TRUE, compile=TRUE, open=TRUE){
  
  if(class(ctmodel) == 'ctStanModel') ctmodel <- c(ctmodel,listOfMatrices(ctmodel$pars)) else if(class(ctmodel) != 'ctsemInit') stop('not a ctsem model!')
  
  bmatrix = function(x, digits=NULL,nottext=FALSE, ...) {
    if(!is.null(x)){
      if(!nottext){
        for(i in 1:length(x)){
          if(is.na(suppressWarnings(as.numeric(x[i]))) & #if x[i] cannot be numeric and 
              grepl('\\',x[i],fixed=TRUE) == FALSE) {
            x[i]=gsub('_', '\\_',x[i],fixed=TRUE)
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
      
      return(paste0("\\begin{bmatrix}\n",
        paste0(out,collapse=''),
        "\\end{bmatrix}",collapse=''))
    } else ""
  }
  
  W <- diag(0,ctmodel$n.latent)
  diag(W) <- 'u'
  
  out <- paste0(" 
\\documentclass[a4paper]{report}
\\usepackage[margin=1cm]{geometry}
\\usepackage{amsmath} %for multiple line equations
\\usepackage{bm}
\\newcommand{\\vect}[1]{\\boldsymbol{\\mathbf{#1}}}

\\begin{document}
 \\begin{",textsize,"}
  \\begin{align*}
  &\\underbrace{\\mathrm{d}
    ",bmatrix(matrix(paste0(ctmodel$latentNames)))," 
    \\big{(}t\\big{)}}_{\\mathrm{d} \\vect{\\eta} (t)}	=  \\left(
      \\underbrace{
        ",bmatrix(ctmodel$DRIFT),"
      }_{\\underbrace{\\vect{A}}_\\textrm{DRIFT}} \\underbrace{
        ",bmatrix(matrix(paste0(ctmodel$latentNames)))," 
        \\big{(}t\\big{)}
      }_{\\vect{\\eta} (t)}	+ \\underbrace{
        ",bmatrix(ctmodel$CINT),"
      }_{\\underbrace{\\vect{b}}_\\textrm{CINT}} ",
    if(ctmodel$n.TDpred > 0) paste0( "+ \\underbrace{
        ",bmatrix(ctmodel$TDPREDEFFECT),"
      }_{\\underbrace{\\vect{M}}_\\textrm{TDPREDEFFECT}}
      \\underbrace{
        ",bmatrix(matrix(paste0('\\chi_',1:ncol(ctmodel$TDPREDEFFECT))))," 
      }_{\\vect{\\chi} (t)}"),
    "\\right) \\mathrm{d}t \\quad + \\nonumber \\\\ \\\\
    & \\qquad \\qquad \\quad \\underbrace{
      ",bmatrix(ctmodel$DIFFUSION),"
    }_{\\underbrace{\\vect{G}}_\\textrm{DIFFUSION}}
    \\underbrace{\\mathrm{d}
      ",bmatrix(matrix(paste0('W_',1:ncol(ctmodel$DIFFUSION))),nottext=TRUE)," 
      (t)}_{\\mathrm{d} \\vect{W}(t)} \\\\ \\\\
          &\\underbrace{
            ",bmatrix(matrix(paste0('W_',1:ctmodel$n.latent)),nottext=TRUE),"  
            (t)}_{\\vect{W}(t+u)} -  \\underbrace{",bmatrix(matrix(paste0('W_',1:ctmodel$n.latent)),nottext=TRUE),"  
            (t)}_{\\vect{W}(t)} \\sim  \\mathrm{N} \\left(
              ",bmatrix(matrix(0,ctmodel$n.latent,1)),", ",bmatrix(W)," \\right) 
\\end{align*}
\\begin{align*}
    &\\underbrace{
      ",bmatrix(matrix(ctmodel$manifestNames),nottext=FALSE),"  
      (t)}_{\\vect{Y}(t)} = 
        \\underbrace{
          ",bmatrix(ctmodel$LAMBDA)," 
        }_{\\underbrace{\\vect{\\Lambda}}_\\textrm{LAMBDA}} \\underbrace{
          ",bmatrix(matrix(paste0('\\eta_',1:nrow(ctmodel$DRIFT))))," 
          (t)}_{\\vect{\\eta}(t)} +
        \\underbrace{
          ",bmatrix(ctmodel$MANIFESTMEANS)," 
        }_{\\underbrace{\\vect{\\tau}}_\\textrm{MANIFESTMEANS}} + 
        \\underbrace{
          ",bmatrix(matrix(paste0('\\epsilon_',1:ctmodel$n.manifest))),"  
          (t)}_{\\vect{\\epsilon}(t)} \\\\ \\\\
          &\\underbrace{
            ",bmatrix(matrix(paste0('\\epsilon_',1:ctmodel$n.manifest))),"  
            (t)}_{\\vect{\\epsilon}(t)} \\sim  \\mathrm{N} \\left(
              ",bmatrix(matrix(0,ctmodel$n.manifest,1)),"
              ,
              \\underbrace{
                ",bmatrix(ctmodel$MANIFESTVAR),"  
              }_{\\underbrace{\\vect{\\Theta}}_\\textrm{MANIFESTVAR}} \\right) 
      \\end{align*}
      \\end{",textsize,"}
   \\end{document}")
  
  openPDF <- function(f) {
    os <- .Platform$OS.type
    if (os=="windows")
      shell.exec(normalizePath(f))
    else {
      pdf <- getOption("pdfviewer", default='')
      if (nchar(pdf)==0)
        stop("The 'pdfviewer' option is not set. Use options(pdfviewer=...)")
      system2(pdf, args=c(f))
    }
  }
  
  if(tex) {
    oldwd <- getwd()
    setwd(dir = folder)
    write(x = out,file = paste0(filename,'.tex'))
    if(compile){
    try(tools::texi2pdf(file=paste0(filename,'.tex'),clean=TRUE))
    if(open) try(openPDF(paste0(filename,'.pdf')))
    }
    setwd(oldwd)
  }
  return(out)
}


