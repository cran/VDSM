#' @title VDSM-heatmap.
#' @description  Plotting the VDSM-heatmap.
#' @param X A m*p matrix which contains m different p-dimensional models. All the elements are either 0 or 1.
#' @param f A vector with m elements which represent each model's frequency in X.
#' @param p The number of variate in the model.
#' @param Anchor.estimate An estimation for the anchor model.
#' @param Anchor.model A vector containing p elements with either 1 or 0 value and must be found in X. Default is the model with the highest frequency.
#' @param xlim A vector with two elements which determine the range of x-axis in the plot.
#' @param ylim A vector with two elements which determine the range of y-axis in the plot.
#' @param fontsize Customize the size of the font in the circles, default is 1.5.
#' @returns A list with components
#' \item{Heatmap.info}{The table includes all the information about each group, i.e., the total possible number of models in the group and the actual existing number of model in the group.}
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
#' Heatmap_example1 = VDSM_heatmap(X,f,p,Anchor.estimate)
#' Heatmap_example2 = VDSM_heatmap(X,f,p,Anchor.estimate,fontsize=3)
#' Heatmap_example3 = VDSM_heatmap(X,f,p,Anchor.estimate,xlim=c(0,5),ylim=c(0,5),fontsize=3)


VDSM_heatmap <- function(X,f,p,Anchor.estimate,xlim=NULL,ylim=NULL,Anchor.model=NULL,fontsize=NULL){
  Rts=VDSM_scatter_heat(X,f,p,Anchor.estimate,Anchor.model)
  AT=Rts[[1]]
  BT=Rts[[2]]
  CT=Rts[[3]]
  true.x.loc=Rts[[4]]
  x.comb=Rts[[5]]
  combination=Rts[[6]]
  hat.beta.true=Rts[[7]]

  if(is.null(fontsize)){
    fontsize = 1.5
  }

  # heat.x=seq(from=0, to=max(x.comb), length.out = 10)
  heat.x=seq(from=0, to=sum(hat.beta.true), length.out = 10)
  heat.yseq=c()
  ynum=p-length(true.x.loc)+1
  heat.yseq[1:ynum]=AT[1:ynum,4]
  for (i in 1:9){
    start=heat.x[i]
    end=heat.x[i+1]
    heat.xloc=which(AT[,1]>start & AT[,1]<=end)
    if (length(heat.xloc)==0){
      heat.yseq[(ynum*i+1):(ynum*(i+1))]=0
    }else if(length(heat.xloc)==ynum){
      heat.yseq[(ynum*i+1):(ynum*(i+1))]=AT[heat.xloc,4]
    }else{
      target1=matrix(AT[heat.xloc,4], nrow = ynum, ncol = length(heat.xloc)/ynum)
      target1[is.na(target1)]=0
      heat.yseq[(ynum*i+1):(ynum*(i+1))]=apply(target1,1,sum)
    }
  }

  heat.yseq[heat.yseq==0]=NA
  DT=data.frame(xloc=rep(1:10, each=ynum),
                yloc=rep(0:(ynum-1), 10),
                frequency=heat.yseq)


  p4 <- ggplot(DT, aes(xloc, yloc)) +
    theme_classic() +
    xlab( expression(paste(H[w]^{"-"})) ) +
    ylab("") +
    theme_bw() +
    theme(axis.ticks = element_blank(),
          axis.line = element_blank()) +
    geom_tile(aes(fill = frequency), colour = "white", na.rm=TRUE) +
    geom_text(aes(x = xloc, y = yloc, label = frequency), size=fontsize, color="white", na.rm=TRUE) +
    scale_fill_continuous(type = "viridis", limits = c(0, 1)) +
    theme(text=element_text(size=8)
          ,axis.text.y=element_blank()
          ,axis.text.x = element_text(size = 6, vjust = 0.5, hjust = 0.5, angle = 60)
          ,plot.margin=unit(c(0.2,0.5,0.5,0), "cm"))

  p5<-
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

  bottom.6y=c()
  for(i in 1:10){
    b6.loc=which(DT[,1]==i)
    freqsum=DT[b6.loc,3]
    freqsum[is.na(freqsum)]=0
    bottom.6y[i]=sum(freqsum)
  }

  ET=data.frame(xloc=1:10, yloc=bottom.6y)

  p6<-ggplot(ET, aes(x=xloc, y=yloc)) +
    geom_bar(stat = "identity") +
    xlab("") +
    ylab("") +
    scale_y_continuous(limits=c(0,1), breaks=seq(0, 1, by=0.5)) +
    theme_bw() +
    theme(text=element_text(size=8)
          ,axis.text.x=element_blank()
          ,legend.position="none"
          ,axis.text.y = element_text(size = 6)
          ,plot.margin=unit(c(0.5,0.5,-0.2,0.5), "cm")
    )

  heatmap.xname=c("0", rep("",9))
  for(i in 1:9){
    heatmap.xname[i+1]=paste(round(heat.x[c(i,i+1)],2),collapse=",")
  }

  ##########################################
  if(is.null(xlim)){
    limitLeft=0
    xmax=max(DT[which(DT[,3]>0),][,1])
    max.index=min(floor(sum(hat.beta.true)/2),xmax)
    if(max.index<10){
      limitRight=heat.x[max.index+1]
    }else{
      limitRight=heat.x[10]
    }
    xlim=c(limitLeft,limitRight)
  }else{
    LL=xlim[1]
    RR=xlim[2]
    limitLeft=heat.x[which(heat.x>LL)[1]-1]
    limitRight=heat.x[which(heat.x>=RR)[1]]
  }


  if(is.null(ylim)){
    limitLow=0
    ymax=max(DT[which(DT[,3]>0),][,2])
    limitHigh=min(10,ymax+1)
    ylim=c(limitLow,limitHigh)
  }else{
    limitLow=ylim[1]
    limitHigh=ylim[2]
  }

  limits2 <- c(limitLow-0.5, limitHigh+0.5)
  breaks2 <- seq(limits2[1]+0.5, limits2[2]-0.5, by=1)

  limits3 <- c(limitLeft+1-0.5, limitRight+0.5)
  breaks3 <- seq(limits3[1]+0.5, limits3[2]-0.5, by=1)


  heatmap.xname=paste("(",heatmap.xname, sep = "")
  heatmap.xname=paste(heatmap.xname, "]", sep = "")
  heatmap.xname[1]="0"

  # assign common axis to both plots
  p4.common.y <- p4 + scale_y_continuous(limits=limits2, breaks=breaks2) +
    scale_x_continuous(limits=limits3, breaks=breaks3, labels = heatmap.xname[1:max(breaks3)])
  p5.common.y <- p5 + scale_x_continuous(limits=limits2, breaks=breaks2)
  p6.common.x <- p6 + scale_x_continuous(limits=limits3, breaks=breaks3)

  # At this point, they have the same axis, but the axis lengths are unequal, so ...

  # build the plots
  p4.common.y <- suppressWarnings(ggplot_gtable(ggplot_build(p4.common.y)))
  p5.common.y <- suppressWarnings(ggplot_gtable(ggplot_build(p5.common.y)))
  p6.common.x <- suppressWarnings(ggplot_gtable(ggplot_build(p6.common.x)))

  # copy the plot height from p1 to p2
  p5.common.y$heights <- p4.common.y$heights
  p6.common.x$widths <- p4.common.y$widths
  blankPanel<-grid.rect(gp=gpar(col="white"))
  lay2=rbind(c(1,2),c(3,4))
  suppressWarnings(kable(lay2))
  plots2=list(blankPanel, p6.common.x, p5.common.y, p4.common.y)
  grid.arrange(grobs=plots2, layout_matrix=lay2, widths=c(1,4), heights=c(1,4))
  return(list(Heatmap.info=DT,
              Hplus.histogram=BT,
              Hminus.weighted.histogram=ET))
}
