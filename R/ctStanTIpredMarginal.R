#' Plot marginal relationships between covariates and parameters for a ctStanFit object.
#'
#' @param fit ctStanFit object.
#' @param tipred character vector representing which tipreds to use.
#' @param pars Subject level matrices from the ctStanFit output -- e.g, 'DRIFT' or 'DIFFUSION'.
#' @param plot Logical, whether to plot.
#'
#' @return If \code{plot=TRUE}, nothing, otherwise an array that can be used with ctPlotArray.
#' @export
#'
#' @examples
#' if(w32chk()){
#'
#' ctStanTIpredMarginal(ctstantestfit(),pars=c('CINT'),tipred=c('TI1'))
#' }
ctStanTIpredMarginal<-function(fit,tipred,pars, plot=TRUE){
  e<-ctExtract(fit,subjectMatrices = TRUE)
  # qseq <- seq(.01,.99,.01)
  
  dt <- data.table(Parameter = '', TIpred = '', y=0,Covariate=0)

  # if(requireNamespace('quantreg')){
  # qr <- data.table(matrix(0,ncol=length(qseq)))
  # qrnames <- paste0('q',qseq)
  # names(qr) <- qrnames
  # dt <- cbind(dt,qr)
  # } else stop('To use this function, run:  install.packages("quantreg")')
  for(p in pars){
    for(ti in tipred){
      tin <- match(ti,fit$ctstanmodelbase$TIpredNames)
      for(i in 1:dim(e[[p]])[3]){
        for(j in 1:dim(e[[p]])[4]){
          dts <-data.frame(Parameter = paste0(p, '[',i,',',j,']'),
            TIpred = ti, y=c(aaply(e[[p]],c(2,3,4),mean,.drop=FALSE)[,i,j,]), Covariate=c(e$tipreds[1,,tin]))
          
          # if(requireNamespace('quantreg')){
          #   # browser()
          #   qr=data.table(sapply(qseq,function(d){
          #   quantreg::predict.rqss(quantreg::rqss(formula = as.formula('y ~ qss(Covariate)'),data=dts,tau=d),dts)
          # }))
          # names(qr) <- qrnames
          # dts <- cbind(dts,qr)
          # }

          dt <- rbind(dt,dts)
        }
      }
    }
  }
  dt=dt[-1,]
  
  # colours = suppressWarnings(RColorBrewer::brewer.pal(length(tipred),'Set1'))[1:length(tipred)]
  
  if(1==99) y <- Covariate <- Parameter <- TIpred <- NULL
  
 g<- ggplot(data = dt,aes(y=y,x=Covariate,colour=TIpred,fill=TIpred)) +
    theme_minimal() +
    # stat_density_2d(aes(alpha=..nlevel..),linetype='dotted',show.legend = FALSE, contour = TRUE) +
    # stat_bin2d(aes(alpha=..density..,fill=TIpred,colour=TIpred),geom='tile',linetype=0,contour=FALSE,show.legend = FALSE) +
    # layer(geom = 'raster',stat=StatDensity2d,params=list(contour=FALSE,linetype=0,alpha=0.2)
      # ,mapping=aes(alpha=(..ndensity..)),position='identity') +
    # geom_quantile(method = "rqss",aes(alpha=.5-abs(.5-(..quantile..))),quantiles = seq(.01,.99,.01))+
    # stat_summary(geom="ribbon", 
    #   fun.ymin = function(x) stat_quantile(aes(x=x), 0.05), 
    #   fun.ymax = function(x) stat_quantile(aes(x=x), 0.95))
    geom_point(data=dt[sample(1:nrow(dt),min(nrow(dt),300),replace=FALSE),],show.legend = TRUE) +
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
