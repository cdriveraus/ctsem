#' Runs stan, and plots sampling information while sampling.
#'
#' @param iter Number of iterations
#' @param chains Number of chains
#' @param ... All the other regular arguments to stan()
#' @export
#' @examples
#' \dontrun{
#' #### example 1 
#' scode <- "
#' parameters {
#'   real y[2]; 
#' } 
#' model {
#'   y[1] ~ normal(0, .5);
#'   y[2] ~ double_exponential(0, 2);
#' } 
#' "
#' fit1 <- stanWplot(iter = 100000,chains=2,cores=1, model_code = scode)
#' }

stanWplot <- function(iter=2000,chains=4,...){


tmpdir=tempdir()
tmpdir=gsub('\\','/',tmpdir,fixed=TRUE)

windows= Sys.info()[1]=='Windows'

stanplot<-function(chains,seed){
  wd<-  paste0("setwd('",tmpdir,"')")
  
  if(1==99) shiny::runApp(appDir = getwd(), {})
  
  writeLines(text=paste0(wd,'
    seed<-',seed,';
    chains<-',chains,';
    iter<-',iter,';
    
    notyet<-TRUE
    while(notyet==TRUE){
    Sys.sleep(1);
    samps<-try(read.csv(file=paste0(seed,"samples_1.csv"),comment.char="#"),silent=TRUE)
    if(class(samps) != "try-error") notyet<-FALSE
    }
    varnames<-colnames(samps);
    # require(shiny); 
    shiny::runApp(appDir=list(server=function(input, output,session) {
    
    output$chainPlot <- renderPlot({
    parameter<-input$parameter
    begin<-input$begin
    refresh <- input$refresh
    colimport<-rep("NULL",length(varnames))
    colimport[which(varnames %in% parameter)]<-NA
    begin<-input$begin
    samps<-list()
    for(chaini in 1:chains) {
    samps[[chaini]]<-try(read.csv(file=paste0(seed,"samples_",chaini,".csv"),comment.char="#",colClasses = colimport),silent=TRUE)
    if(class(samps[[chaini]])=="try-error") samps[[chaini]]=samps[[1]][1,,drop=FALSE]
}
    
    mini<-min(unlist(lapply(1:chains,function(chaini) samps[[chaini]][-1:-begin,parameter])),na.rm=T)
    maxi<-max(unlist(lapply(1:chains,function(chaini) samps[[chaini]][-1:-begin,parameter])),na.rm=T)
    lengthi<-max(unlist(lapply(1:chains,function(chaini) length(samps[[chaini]][,parameter]))),na.rm=T) #-1:-begin
    
    plot((begin):(lengthi+begin-1),
    c(samps[[1]][-1:-begin,parameter],rep(NA,lengthi-length(samps[[1]][-1:-begin,parameter]))),
    type="l",xlab="",ylab="",main=parameter,
    log=ifelse(parameter %in% c("stepsize__"),"y",""),
    xlim=c(begin,lengthi),
    ylim=c(mini,maxi)
    )
    
    if(chains > 1) for(chaini in 2:chains){
    points(begin:(lengthi+begin-1),c(samps[[chaini]][-1:-begin,parameter],rep(NA,lengthi-length(samps[[chaini]][-1:-begin,parameter]))),type="l",xlab="",ylab="",main=parameter,col=chaini)
    }
    grid()
    
    })
    },ui=shiny::fluidPage(
    # Application title
    shiny::titlePanel("stan mid-sampling plots..."),
    shiny::sidebarLayout(
    # Sidebar with a slider input for number of observations
    shiny::sidebarPanel(
    shiny::sliderInput("begin", "Start of range:", min = 1,max=iter,value = 1,step=1), 
    shiny::selectInput("parameter", "Choose a parameter:", choices = varnames),
    shiny::actionButton("refresh", "Refresh sample data")
    ),
    
    # Show a plot of the generated distribution
    shiny::mainPanel(
    shiny::plotOutput("chainPlot")
    )
    ))),
    launch.browser=TRUE)
    quit(save="no")'),con=paste0(tmpdir,"/stanplottemp.R"))
  
  if(windows) system(paste0("Rscript --slave --no-restore -e source(\'",tmpdir,"/stanplottemp.R\')"),wait=FALSE) else
    system(paste0(R.home(component = "home"),"/Rscript --slave --no-restore -e source\\(\\\'",tmpdir,"\\/stanplottemp.R\\\'\\)"),wait=FALSE)
  
}

stanseed<-floor(as.numeric(Sys.time()))

  sample_file<-paste0(tmpdir,'/',stanseed,'samples', ifelse(chains==1,'_1',''),'.csv')

  stanplot(chains=chains,seed=stanseed)
  
  out=sampling(iter=iter,chains=chains,sample_file=sample_file,...)

  for(chaini in 1:chains) system(paste0("rm ",tmpdir,'/',stanseed,"samples_",chaini,".csv"))
  system(paste0('rm ',tmpdir,'/stanplottemp.R'))
  return(out)
}
