plotdensity2dby2  <- function(x1,y1,x2=NA,y2=NA,xlab,ylab,title,
  grouplab=c('Observed','Model'),
  grouppch=c(20,19),
  colours=c('blue','red'),
  group1samples = 300,
  group2samples=50000,resolution=100){
  
  if(1==99) Source <- ..density.. <- ..level.. <- ..b <- y <- x <- ..ndensity.. <- ..nlevel.. <-  NULL
  #remove missings
  
  nomiss <- intersect(which(!is.na(x1)),which(!is.na(y1)))
  y1 <- y1[nomiss]
  x1 <- x1[nomiss]
  
  
  
  if(length(y1) > group1samples) {
    datasample <- sample(1:length(y1),group1samples) 
    y1 <- y1[datasample]
    x1 <- x1[datasample]
  } 
  
  pd <- data.table(Source=grouplab[1],
    x=x1,
    y=y1)
  
  if(any(!is.na(x2))){
    nomiss <- intersect(which(!is.na(x2)),which(!is.na(y2)))
    y2 <- y2[nomiss]
    x2 <- x2[nomiss]
    
    if(length(y2) > group2samples) {
      datasample <- sample(1:length(y2),group2samples) 
      y2 <- y2[datasample]
      x2 <- x2[datasample]
    }
    
    
    pd <- rbind(pd,data.table(Source=grouplab[2], 
      x=x2,
      y=y2))
  }
  
  lims <- lapply(c('x','y'),function(b) sapply(grouplab,function (a) quantile(pd[Source==a,..b],c(.05,.95),na.rm=TRUE)))
  names(lims) <- c('x','y')
  lims <- lapply(lims, function(x) c(min(x[1,]),max(x[2,])))
  lims <- lapply(lims,function(x) x + c(-1,1) *sd(x)/10)
  
  
  
  # pd$xd <- as.numeric(cut_width(pd$x, diff(range(pd$x))/30))

  g<-ggplot(data=pd,aes(y=y,x=x,shape=Source,colour=Source)) 
  
  if(all(!is.na(x2))) g <- g+ stat_density_2d(data = subset(pd, Source==grouplab[2]),
    geom="raster",interpolate=TRUE,
    aes(alpha=..ndensity..,fill = ..ndensity..),
    show.legend = FALSE, contour = FALSE,n = resolution, 
    h=c(MASS::bandwidth.nrd(subset(pd, Source==grouplab[2])$x)*1.5,
      MASS::bandwidth.nrd(subset(pd, Source==grouplab[2])$y)*1.5)
  ) +
    scale_fill_gradient (low = "white", high = colours[2],guide='none')
    
    g <- g+ 
    stat_density_2d(data = subset(pd, Source==grouplab[1]),
      aes(alpha=..nlevel..),linetype='dotted',show.legend = FALSE, contour = TRUE,
      h=c(MASS::bandwidth.nrd(subset(pd, Source==grouplab[1])$x)*1.5,
        MASS::bandwidth.nrd(subset(pd, Source==grouplab[1])$y)*1.5)) +
    coord_cartesian(xlim = lims$x, ylim = lims$y) +
    geom_point(data = subset(pd, Source==grouplab[1]),show.legend = TRUE) +
    labs(x=xlab,y=ylab,title = title, color  = "Source", shape = "Source")+ 
    theme_minimal() + 
    scale_alpha(guide = 'none') +
    scale_colour_manual(values=setNames(colours, grouplab))+
    scale_shape_manual(values=setNames(grouppch,grouplab)) +
    theme(legend.title = element_blank())
  
}
