#' @title Gplot.
#' @description Plotting Gplot.
#' @param X A m*p matrix which contains m different p-dimensional models. All the elements are either 0 or 1.
#' @param f A vector with m elements which represent each model's frequency in X.
#' @param p The number of variate in the model.
#' @param Anchor.model A vector containing p elements with either 1 or 0 value and must be found in X. Default is the model with the highest frequency.
#' @param xlim A vector with two elements which determine the range of x-axis in the plot.
#' @param ylim A vector with two elements which determine the range of y-axis in the plot.
#' @param circlesize customize the size of the circle in the plot, default is 10.
#' @param linewidth Customize the width of the line in the plot, default is 1.
#' @param fontsize Customize the size of the font in the circles, default is 1.5.
#' @returns A list with components
#' \item{Gplot.info}{The table includes all the information about each group, i.e., the total possible number of models in the group and the actual existing number of model in the group.}
#' \item{MC.histogram}{The frequency of model complexity.}
#' \item{HD.histogram}{The frequency of Hamming distance.}
#' @export
#' @importFrom plyr count
#' @importFrom dplyr arrange full_join left_join
#' @importFrom ggplot2 ggplot geom_line geom_text geom_point theme_bw theme geom_bar coord_flip scale_y_reverse scale_y_continuous coord_cartesian scale_x_continuous ggplot_gtable ggplot_build aes element_text element_blank xlab ylab theme_classic geom_tile scale_fill_continuous unit
#' @importFrom viridis scale_color_viridis
#' @importFrom grid grid.rect gpar
#' @importFrom gridExtra grid.arrange
#' @importFrom knitr kable
#' @importFrom stats aggregate coef frequency
#' @examples
#' data(exampleX)
#' X=exampleX
#' data(examplef)
#' f=examplef
#' p=8
#' G_example1 = Gplot(X,f,p)
#' G_example2 = Gplot(X,f,p,xlim=c(0,7),ylim=c(3,8))
#' G_example3 = Gplot(X,f,p,xlim=c(0,7),ylim=c(3,8),circlesize=15,linewidth=2,fontsize=3)

Gplot <- function(X,f,p,Anchor.model=NULL,xlim=NULL,ylim=NULL,circlesize=NULL,linewidth=NULL,fontsize=NULL){
  Rts=Groupinfo(X,f,p,Anchor.model)
  final=Rts[[1]]
  FreqModel=Rts[[2]]
  MC.anchor=Rts[[3]]

  ########################################################################
  ## Modify the user-specified value of ylim, default is MC.anchor+-5
  if(is.null(ylim)){
    MClength=5
    if(MC.anchor + MClength > p){
      limitH=p
      limitL=p-2*MClength
    }else if(MC.anchor - MClength < 0){
      limitH=2*MClength
      limitL=0
    }else{
      limitH=MC.anchor + MClength
      limitL=MC.anchor - MClength
    }
  }else if(length(ylim)!=2){
    stop("ylim should be a vector with 2 elements or leave it NULL.")
  }else{
    limitL=ylim[1]
    limitH=ylim[2]
  }

  if(limitL<0){
    limitL=0
  }
  if(limitH>p){
    limitH=p
  }
  ylim=c(limitL,limitH)

  ########################################################################
  ## Modify the user-specified value of xlim
  if(is.null(xlim)){
    limitLeft=0
    limitRight=limitH-limitL
  }else if(length(xlim)!=2){
    stop("xlim should be a vector with 2 elements or leave it NULL.")
  }else{
    limitLeft=xlim[1]
    limitRight=xlim[2]
  }

  if(limitLeft<0){
    limitLeft=0
  }
  if(limitRight>p){
    limitRight=p
  }
  xlim=c(limitLeft,limitRight)

  ########################################################################
  if(limitL==0){
    set.L=0
  }else{
    set.L=limitL-1
  }

  if(limitH==p){
    set.H=p
  }else{
    set.H=limitH+1
  }

  if(limitLeft==0){
    set.Left=0
  }else{
    set.Left=limitLeft-1
  }

  if(limitRight==p){
    set.Right=p
  }else{
    set.Right=limitRight+1
  }
  ########################################################################
  need.loc=which(final[,2]>=set.L & final[,2]<=set.H & final[,3]>=set.Left & final[,3]<=set.Right)
  Gplot.info=final[need.loc,]

  Gleft.hist=c()
  for(i in set.L:set.H){
    left.loc=which(FreqModel[,1]==i)
    Gleft.hist[i+1]=sum(FreqModel[left.loc,3])
  }
  G.left=data.frame(yloc=set.L:set.H,
                    xloc=Gleft.hist[(set.L+1):(set.H+1)])


  Gup.hist=c()
  for(i in set.Left:set.Right){
    up.loc=which(FreqModel[,2]==i)
    Gup.hist[i+1]=sum(FreqModel[up.loc,3])
  }

  G.up=data.frame(xloc=set.Left:set.Right,
                  yloc=Gup.hist[(set.Left+1):(set.Right+1)])





  Gline.x=c()
  Gline.y=c()
  G.line.group=c()
  for (i in set.Left:set.Right){
    start.loc=which(Gplot.info[,3]==i)
    end.loc=which(Gplot.info[,3]==(i+1))
    end.poss.y=Gplot.info[end.loc,2]
    for (j in 1:length(start.loc)){
      vec1x=c(Gplot.info[start.loc,3][j])
      vec1y=c(Gplot.info[start.loc,2][j])
      if(any(end.poss.y-vec1y==1)==TRUE){
        vec2x=c(Gplot.info[end.loc[which(end.poss.y-1==vec1y)],3])
        vec2y=c(Gplot.info[end.loc[which(end.poss.y-1==vec1y)],2])
        Gline.x=c(Gline.x,vec1x,vec2x)
        Gline.y=c(Gline.y,vec1y,vec2y)
      }
      if(any(end.poss.y-vec1y==-1)==TRUE){
        vec3x=c(Gplot.info[end.loc[which(end.poss.y+1==vec1y)],3])
        vec3y=c(Gplot.info[end.loc[which(end.poss.y+1==vec1y)],2])
        Gline.x=c(Gline.x,vec1x,vec3x)
        Gline.y=c(Gline.y,vec1y,vec3y)
      }
    }
  }

  G.line.group=rep(1:(length(Gline.x)/2), each=2)
  G.line=data.frame(linex=Gline.x, liney=Gline.y, linegroup=G.line.group)

  if(is.null(circlesize)){
    circlesize = 10
  }
  if(is.null(linewidth)){
    linewidth = 1
  }
  if(is.null(fontsize)){
    fontsize = 1.5
  }

  p1<-
    ggplot(Gplot.info, aes(x=HD, y=MC, colour=freq)) +
    geom_line(data = G.line, aes(linex, liney, group = linegroup),
              size = linewidth, alpha = 0.8, color = '#8BB6D6') +
    geom_point(shape=16, size=circlesize) +
    geom_text(aes(x = HD, y = MC, label = freq), size=fontsize, color="white", na.rm=TRUE) +
    xlab("HD") +
    ylab("") +
    theme_bw()+
    scale_color_viridis("frequency", option = "D", limits = c(0, 1)) +
    theme(text=element_text(size=8)
          ,axis.text.y=element_blank()
          ,plot.margin=unit(c(0.2,0.5,0.5,0), "cm"))


  p2<-ggplot(G.left, aes(x=yloc, y=xloc)) +
    geom_bar(stat = "identity") +
    xlab("MC") +
    ylab("frequency") +
    theme_bw() +
    coord_flip()+scale_y_reverse(limits=c(1,0), breaks=seq(1, 0, by=-0.5)) +
    theme(text=element_text(size=8)
          ,legend.position="none"
          ,axis.text.x = element_text(size = 6)
          ,plot.margin=unit(c(0.5,0.1,0.5,0.5), "cm")
    )



  p3<-ggplot(G.up, aes(x=xloc, y=yloc)) +
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

  limits <- c(limitL-0.5, limitH+0.5)
  breaks <- seq(limits[1]+0.5, limits[2]-0.5, by=1)

  limits1 <- c(limitLeft-0.5, limitRight+0.5)
  breaks1 <- seq(limits1[1]+0.5, limits1[2]-0.5, by=1)

  # assign common axis to both plots
  p1.common.y <- p1 + coord_cartesian(ylim = limits, xlim = limits1) +
    scale_y_continuous(breaks=breaks) +
    scale_x_continuous(breaks=breaks1)
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
  grid.arrange(grobs=plots1, layout_mGrix=lay1, widths=c(1,4), heights=c(1,4))
  return(list(Gplot.info=Gplot.info,
              MC.histogram=G.left,
              HD.histogram=G.up))
}
