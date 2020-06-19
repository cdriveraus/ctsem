ctStanSubjectTIpred <- function(fit){
  if(class(fit$stanfit) !='stanfit' && class(fit$stanfit$stanfit) %in% 'stanfit'){
    sub <- ctStanKalman(fit,collapsefunc = mean,subjectpars = TRUE,tformsubjectpars=FALSE)
    sub <- data.table(matrix(sub,nrow=dim(sub)[2],ncol=dim(sub)[3],
      dimnames = dimnames(sub)[-1]))

    dat <- data.frame(sub, fit$standata$tipredsdata)
  }
}

ctStanIndPars <- function(fit){
  if(class(fit$stanfit) !='stanfit' && class(fit$stanfit$stanfit) %in% 'stanfit'){
    base <- 
    sub <- ctStanKalman(fit,collapsefunc = mean,subjectpars = TRUE)
    sub <- data.table(matrix(sub,nrow=dim(sub)[2],ncol=dim(sub)[3],
      dimnames = dimnames(sub)[-1]))
    
    dat <- data.frame(sub, fit$standata$tipredsdata)
  }
}



#' Plot marginal relationships between covariates and parameters for a ctStanFit object.
#'
#' @param fit ctStanFit object.
#' @param tipred character vector representing which tipreds to use.
#' @param pars Character vector of either parameter names, or if matrices=TRUE,
#' system matrices with optional indices -- e.g, 'DRIFT' or 'DIFFUSION[2,1]'.
#' @param matrices Logical, see pars argument.
#' @param plot Logical, whether to plot.
#'
#' @return If \code{plot=TRUE}, nothing, otherwise an array that can be used with ctPlotArray.
#' @export
#'
#' @examples
#' if(w32chk()){#'
#' ctStanTIpredMarginal(ctstantestfit(),pars=c('CINT'),
#'   matrices=TRUE, tipred=c('TI1'))
#' }
ctStanTIpredMarginal<-function(fit,tipred,pars,matrices=TRUE, plot=TRUE){
  
  if(class(fit$stanfit) =='stanfit'){ #very slow!
    e<-ctExtract(fit,subjectMatrices = TRUE)
    dt <- data.table(Parameter = '', TIpred = '', y=0,Covariate=0)
    for(p in pars){
      for(ti in tipred){
        tin <- match(ti,fit$ctstanmodelbase$TIpredNames)
        for(i in 1:dim(e[[p]])[3]){
          for(j in 1:dim(e[[p]])[4]){
            dts <-data.frame(Parameter = paste0(p, '[',i,',',j,']'),
              TIpred = ti, ParValue=c(aaply(e[[p]],c(2,3,4),mean,.drop=FALSE)[,i,j,]), Covariate=c(fit$standata$tipredsdata[,tin]))
            
            dt <- rbind(dt,dts)
          }
        }
      }
    }
    dt=dt[-1,]
  }
  
  
  
  
  if(1==99) ..=NULL
  if(class(fit$stanfit) !='stanfit' && class(fit$stanfit$stanfit) %in% 'stanfit'){
    sub <- ctStanKalman(fit,collapsefunc = mean,subjectpars = TRUE)
    sub <- data.table(matrix(sub,nrow=dim(sub)[2],ncol=dim(sub)[3],
      dimnames = dimnames(sub)[-1]))
    
    if(matrices) subcols=vgrep(pars,colnames(sub))
    
    if(!matrices){ #find and rename columns to represent parameter names
      subcols = findmatrixslots(pars,listOfMatrices(fit$ctstanmodelbase$pars))
      for(pi in 1:length(subcols)){
        colnames(sub)[colnames(sub) %in% subcols[[pi]]] <- names(subcols)[pi]
      }
      subcols <- colnames(sub) %in% pars
      if(any(!pars %in% colnames(sub))) stop('Invalid par names!')
    }
    
    if(FALSE) ..subcols <- 1
    sub <- sub[,..subcols]

      tidat <- fit$standata$tipredsdata[,vgrep(tipred,colnames(fit$standata$tipredsdata)),drop=FALSE]
      tidat <- melt(data.table(id=fit$setup$idmap[,1],tidat),
        measure.vars=colnames(tidat),id.vars = 'id',variable.name = 'TIpred',value.name = 'Covariate')
      sub2 <- melt(data.table(id=fit$setup$idmap[,1],sub),
        measure.vars=colnames(sub),idvars='id',variable.name = 'Parameter',value.name = 'ParValue')
      
      dt <- merge(tidat,sub2,all = TRUE,allow.cartesian = TRUE)
      # browser()
    }
    
    
    # colours = suppressWarnings(RColorBrewer::brewer.pal(length(tipred),'Set1'))[1:length(tipred)]
    
    if(1==99) ParValue <- Covariate <- Parameter <- TIpred <- NULL
    
    g<- ggplot(data = dt,aes(y=ParValue,x=Covariate,colour=TIpred,fill=TIpred)) +
      theme_minimal() +
      # stat_density_2d(aes(alpha=..nlevel..),linetype='dotted',show.legend = FALSE, contour = TRUE) +
      # stat_bin2d(aes(alpha=..density..,fill=TIpred,colour=TIpred),geom='tile',linetype=0,contour=FALSE,show.legend = FALSE) +
      # layer(geom = 'raster',stat=StatDensity2d,params=list(contour=FALSE,linetype=0,alpha=0.2)
      # ,mapping=aes(alpha=(..ndensity..)),position='identity') +
      # geom_quantile(method = "rqss",aes(alpha=.5-abs(.5-(..quantile..))),quantiles = seq(.01,.99,.01))+
      # stat_summary(geom="ribbon", 
      #   fun.ymin = function(x) stat_quantile(aes(x=x), 0.05), 
      #   fun.ymax = function(x) stat_quantile(aes(x=x), 0.95))
      geom_point(alpha=.5,
        # data=dt[sample(1:nrow(dt),min(nrow(dt),300),replace=FALSE),],
        show.legend = TRUE) +
      # geom_point(data = subset(pd, Source==grouplab[1]),show.legend = TRUE) +
      # scale_fill_gradient (low = "white", high = colours[1,guide='none')+
      # scale_fill_manual(values=setNames(colours,tipred)) + #values=c("a"="#FF0000", "b"="#00FF00")) +
      # labs(x=xlab,y=ylab,title = title, color  = "Source", shape = "Source")+ 
      theme_minimal() +
      scale_alpha(guide = 'none') +
      # scale_colour_manual(values=setNames(colours, grouplab))+
      # scale_shape_manual(values=setNames(grouppch,grouplab)) +
      theme(legend.title = element_blank()
        ,panel.background=element_rect(fill="transparent",colour=NA),
        plot.background=element_rect(fill="transparent",colour=NA),
        legend.key = element_rect(fill = "transparent", colour = "transparent")
        # ,panel.background = element_rect(fill = "white",colour = 'white'), # or theme_blank()
        # panel.grid.minor = theme_blank(), 
        # panel.grid.major = theme_blank(),
        # plot.background = element_rect(fill = "white",colour = 'white')
      ) +
      facet_wrap(vars(Parameter),scales='free')
    
    # g<-g+geom_quantile(quantiles=qseq,method='rqss')
    
    # for(i in 1:(length(qrnames)/2)){
    #   g <- g+ geom_ribbon(mapping=aes_string(ymin=qrnames[i],ymax=qrnames[i+(length(qrnames)/2)]),linetype=0,alpha=.01)
    # }
    
    
    if(plot) print(g)
    
  }
  
