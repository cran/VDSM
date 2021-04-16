#' @title VDSM-Scatterplot.
#' @description Plotting the VDSM-Scatterplot.
#' @param X A m*p matrix which contains m different p-dimensional models. All the elements are either 0 or 1.
#' @param f A vector with m elements which represent each model's frequency in X.
#' @param p The number of variate in the model.
#' @param Anchor.estimate An estimation for the anchor model.
#' @param Anchor.model A vector containing p elements with either 1 or 0 value and must be found in X. Default is the model with the highest frequency.
#' @param xlim A vector with two elements which determine the range of x-axis in the plot.
#' @param ylim A vector with two elements which determine the range of y-axis in the plot.
#' @param circlesize customize the size of the circle in the plot, default is 10.
#' @param fontsize Customize the size of the font in the circles, default is 1.5.
#' @returns A list with components
#' \item{Scatterplot.info}{The table includes all the information about each group, i.e., the total possible number of models in the group and the actual existing number of model in the group.}
#' \item{Hplus.histogram}{The frequency of Hamming distance plus.}
#' \item{Hminus.weighted.histogram}{The frequency of Hamming distance minus-weighted.}
#' @export
#' @examples
#' data(exampleX)
#' X=exampleX
#' data(examplef)
#' f=examplef
#' p=8
#' Anchor.estimate=c(3,2.5,2,1.5,1,0,0,0)
#' Scatter_example1 = VDSM_scatterplot(X,f,p,Anchor.estimate)
#' Scatter_example2 = VDSM_scatterplot(X,f,p,Anchor.estimate,xlim=c(0,5),
#' ylim=c(0,8),circlesize=15,fontsize=2)


VDSM_scatterplot <- function(X,f,p,Anchor.estimate,xlim=NULL,ylim=NULL,Anchor.model=NULL,circlesize=NULL,fontsize=NULL){
  Rts=VDSM_scatter_heat(X,f,p,Anchor.estimate,Anchor.model)
  AT=Rts[[1]]
  BT=Rts[[2]]
  CT=Rts[[3]]
  true.x.loc=Rts[[4]]
  x.comb=Rts[[5]]
  combination=Rts[[6]]
  hat.beta.true=Rts[[7]]


  if(is.null(circlesize)){
    circlesize = 10
  }
  if(is.null(fontsize)){
    fontsize = 1.5
  }
  AT2=AT[which(AT[,4]>0),]
  p1<-
    ggplot(AT2, aes(x=xloc, y=yloc, colour=frequency)) +
    geom_point(shape=16, size=circlesize) +
    geom_text(aes(x = xloc, y = yloc, label = round(frequency,3)), size=fontsize, color="white", na.rm=TRUE) +
    xlab(expression(paste(H[w]^{"-"}))) +
    ylab("") +
    theme_bw() +
    scale_color_viridis("frequency", option = "D", limits = c(0, 1)) +
    theme(text=element_text(size=8)
          ,axis.text.y=element_blank()
          ,axis.text.x = element_text(size = 6, vjust = 0.5, hjust = 0.5, angle = 90)
          ,plot.margin=unit(c(0.2,0.5,0.5,0), "cm"))


  p2<-
    ggplot(BT, aes(x=yloc, y=xloc)) +
    geom_bar(stat = "identity") +
    xlab(expression(paste(H^{"+"}))) +
    ylab("frequency") +
    theme_bw() +
    coord_flip()+scale_y_reverse(limits=c(1,0), breaks=seq(1, 0, by=-0.5)) +
    theme(text=element_text(size=8)
          ,axis.text.x = element_text(size = 6)
          ,legend.position="none"
          ,plot.margin=unit(c(0.5,0.1,0.5,0.5), "cm")
    )
  ##########################################
  if(is.null(xlim)){
    limitLeft=0
    limitRight=min(CT[floor(nrow(CT)/3),1],max(AT2[,1]))
    xlim=c(limitLeft,limitRight)
  }else{
    limitLeft=xlim[1]
    limitRight=xlim[2]
  }
  if(is.null(ylim)){
    limitLow=0
    limitHigh=min(10,max(AT2[,2]))
    ylim=c(limitLow,limitHigh)
  }else{
    limitLow=ylim[1]
    limitHigh=ylim[2]
  }
  CT=CT[which(CT[,1]>=limitLeft),]
  CT=CT[which(CT[,1]<=limitRight),]

  p3<-ggplot(CT, aes(x=xloc, y=yloc)) +
    geom_line() +
    geom_point(data = CT[which(CT[,2]>0),], shape=1, size=2, color="black", na.rm=TRUE) +
    xlab("") +
    ylab("") +
    scale_y_continuous(limits=c(0,1), breaks=seq(0, 1, by=0.5)) +
    theme_bw() +
    theme(text=element_text(size=8)
          ,axis.text.x=element_blank()
          ,legend.position="none"
          ,axis.text.y=element_text(size = 6)
          ,plot.margin=unit(c(0.5,0.5,-0.2,0.5), "cm")
    )

  limits <- c(limitLow-0.5, limitHigh+2+0.5)
  breaks <- seq(limits[1]+0.5, limits[2]-0.5, by=1)

  x4.max=max(AT2[,1])
  limits1 <- c(limitLeft, limitRight)
  breaks1 <- x.comb[unique(AT2[,3])]

  # assign common axis to both plots
  targ=unique(AT2[,3])[-1]
  if (length(targ)==0){
    sd.xname="0"
  }else{
    sd.xname=c("0", rep("",length(targ)))
    for(i in 1:length(targ)){
      if(length(combination[[targ[i]]])==1){
        sd.xname[i+1]=paste("b",combination[[targ[i]]],sep = "")
      }else{
        sd.xname[i+1]=paste("b",combination[[targ[i]]],sep = "",collapse="+")
      }
    }
  }

  anno<-data.frame(x = breaks1, y = rep(limitHigh+1,length(breaks1)), label = sd.xname)
  anno=anno[match(unique(anno[,1]),anno[,1]),]
  p1.common.y <- p1 + scale_y_continuous(limits=limits, breaks=breaks) +
    scale_x_continuous(limits=limits1, breaks=breaks1, labels = round(breaks1,2)) +
    geom_text(data=anno, aes(x = x, y = y, label = label), size=3, angle=90, color="black")
  p2.common.y <- p2 + scale_x_continuous(limits=limits, breaks=breaks)
  p3.common.x <- p3 + scale_x_continuous(limits=limits1, breaks=breaks1)

  # At this point, they have the same axis, but the axis lengths are unequal, so ...

  # build the plots
  p1.common.y <- suppressWarnings(ggplot_gtable(ggplot_build(p1.common.y)))
  p2.common.y <- suppressWarnings(ggplot_gtable(ggplot_build(p2.common.y)))
  p3.common.x <- suppressWarnings(ggplot_gtable(ggplot_build(p3.common.x)))

  # copy the plot height from p1 to p2
  p2.common.y$heights <- p1.common.y$heights
  p3.common.x$widths <- p1.common.y$widths
  blankPanel<-grid.rect(gp=gpar(col="white"))

  lay1=rbind(c(1,2),c(3,4))
  suppressWarnings(kable(lay1))
  plots1=list(blankPanel, p3.common.x, p2.common.y, p1.common.y)
  grid.arrange(grobs=plots1, layout_matrix=lay1, widths=c(1,4), heights=c(1,4))
  return(list(Scatterplot.info=AT2,
              Hplus.histogram=BT,
              Hminus.weighted.histogram=CT))
}
