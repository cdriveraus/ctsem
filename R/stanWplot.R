#' Runs stan, and plots sampling information while sampling.
#'
#' @param object stan model object
#' @param iter Number of iterations
#' @param chains Number of chains
#' @param ... All the other regular arguments to stan()
#' @export
#' @details On windows, requires Rtools installed and able to be found by pkgbuild::rtools_path()
#' @examples
#' library(rstan)
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
#' #Uncomment the following lines -- launches rscript not compatible with cran check.
#' #sm <- stan_model(model_code = scode)
#' #fit1 <- stanWplot(object = sm,iter = 100000,chains=2,cores=1)

stanWplot <- function(object,iter=2000,chains=4,...){
  
    tmpdir=tempdir()
    tmpdir=gsub('\\','/',tmpdir,fixed=TRUE)
    windows= Sys.info()[1]=='Windows'
  
    stanplot<-function(chains,seed,iter){
      wd<-  paste0("setwd('",tmpdir,"')")
      
      writeLines(text=paste0(wd,'
    seed<-',seed,';
    chains<-',chains,';
    iter<-',iter,';
    
    notyet<-TRUE
     while(any(notyet==TRUE)){
      Sys.sleep(1);
       cmd = paste0(\'findstr "^[^#]" \', seed,"samples_1.csv") #needed because of comments after warmup
      samps<-try(data.table::fread(skip="lp__",
        cmd=cmd,
      # file=paste0(seed,"samples_1.csv"),
        ),silent=TRUE)
      if(!"try-error" %in% class(samps) && length(samps) > 0) notyet<-FALSE
     }
    
    varnames<-colnames(samps);
    
    
    
    server=function(input, output,session) {
      # 
      session$onSessionEnded(function() {
        stopApp()
      })

      # observe({ # stop shiny
      #   if (all(input$nsamps %in% c(0,input$iter))) close.window()#stopApp()
      #   if(all(unlist(lapply(inputsamps,class))=="try-error")) window.close()#stopApp()
      # })
    
    output$chainPlot <- renderPlot({
    parameter<-input$parameter
    refresh <- input$refresh
    begin<-input$begin
    samps<-list()
    
    colclasses <- rep("NULL",length(varnames))
colclasses[which(varnames %in% parameter)] <- NA
        
        nsamps <- rep(NA,chains)
    for(chaini in 1:chains) {
      cmd = paste0(\'findstr "^[^#]" \', seed,"samples_",chaini,".csv") #needed because of comments after warmup
      samps[[chaini]]<-try(
        # as.matrix(read.table(
        #   file=paste0(seed,"samples_",chaini,".csv"),sep=",",comment.char="#",
        #   header = TRUE,colClasses = colclasses)
        
        data.table::fread(cmd = cmd,
          colClasses = colclasses,
          select = parameter,
          skip="lp__")
      ,silent=TRUE)
      if("try-error" %in% class(samps[[chaini]]) || nrow(samps[[chaini]]) ==0) {
       samps[[chaini]]=samps[[1]][1,,drop=FALSE]
       nsamps[chaini] <- NA
      } else  nsamps[chaini] <- nrow(samps[[chaini]])
    }
    
      # Check if all chains have reached the iter value, and if so, close the browser window
  observe({
    if (all(nsamps == iter)) {
      for(chaini in 1:chains){
        system(paste0("rm ",tmpdir,"/",stanseed,"samples_",chaini,".csv"))
      }
      runjs("window.close();")
      
    }
  })
 
    
    mini<-min(unlist(lapply(1:chains,function(chaini) samps[[chaini]][begin:nrow(samps[[chaini]]),parameter,with=FALSE])),na.rm=T)
    maxi<-max(unlist(lapply(1:chains,function(chaini) samps[[chaini]][begin:nrow(samps[[chaini]]),parameter,with=FALSE])),na.rm=T)
    lengthi<-max(unlist(lapply(1:chains,function(chaini) nrow(samps[[chaini]][,parameter,with=FALSE]))),na.rm=TRUE) #-1:-begin
    
    plot(begin:nrow(samps[[1]]),
      unlist(samps[[1]][begin:nrow(samps[[1]]),parameter,with=FALSE]),
      type="l",xlab="",ylab="",main=parameter,
      log=ifelse(parameter %in% c("stepsize__"),"y",""),
      xlim=c(begin,lengthi),
      ylim=c(mini,maxi)
    )
    
    if(chains > 1) for(chaini in 2:chains){
      points(begin:nrow(samps[[chaini]]),
      unlist(samps[[chaini]][begin:nrow(samps[[chaini]]),parameter,with=FALSE]), type="l",xlab="",ylab="",main=parameter,col=chaini)
    }
    grid()
    
    })
    } #end server function
    
    
    
    
    ui=shiny::fluidPage(
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
    )) #end ui function
    
    
    shiny::runApp(appDir=list(server=server,ui=ui), launch.browser=TRUE)
    quit(save="no")'),con=paste0(tmpdir,"/stanplottemp.R"))
      
      if(windows) system(paste0("Rscript --slave --no-restore -e source(\'",tmpdir,"/stanplottemp.R\')"),wait=FALSE) else
        system(paste0(R.home(component = "home"),"/Rscript --slave --no-restore -e source\\(\\\'",tmpdir,"\\/stanplottemp.R\\\'\\)"),wait=FALSE)
      
    }
    
    stanseed<-floor(as.numeric(Sys.time()))
    
    on.exit({
      for(chaini in 1:chains) system(paste0("rm ",tmpdir,'/',stanseed,"samples_",chaini,".csv"))
      system(paste0('rm ',tmpdir,'/stanplottemp.R'))
    })
    
    
    sample_file<-paste0(tmpdir,'/',stanseed,'samples', ifelse(chains==1,'_1',''),'.csv')
    
    if(requireNamespace(c('shiny')))  stanplot(chains=chains,seed=stanseed,iter=iter)
    
    # out=rstan::sampling(object=object,iter=iter,chains=chains,sample_file=sample_file)
    out=rstan::sampling(object=object,iter=iter,chains=chains,sample_file=sample_file,...)
    
    return(out)
  
}
