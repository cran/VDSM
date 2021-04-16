#' @title Hplot.
#' @description Plotting Hplot.
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
#' \item{Hplot.info}{The table includes all the information about each group, i.e., the total possible number of models in the group and the actual existing number of model in the group.}
#' \item{Hplus.histogram}{The frequency of Hamming distance plus.}
#' \item{Hminus.histogram}{The frequency of Hamming distance minus.}
#' @export
#' @examples
#' data(exampleX)
#' X=exampleX
#' data(examplef)
#' f=examplef
#' p=8
#' H_example1 = Hplot(X,f,p)
#' H_example2 = Hplot(X,f,p,xlim=c(0,4),ylim=c(0,2))
#' H_example3 = Hplot(X,f,p,xlim=c(0,4),ylim=c(0,2),circlesize=15,linewidth=2,fontsize=3)

Hplot <- function(X,f,p,Anchor.model=NULL,xlim=NULL,ylim=NULL,circlesize=NULL,linewidth=NULL,fontsize=NULL){
  Rts=Groupinfo(X,f,p,Anchor.model)
  MC.anchor=Rts[[3]]
  Summary1=Rts[[4]]


  Hplot.info=data.frame(lose=rep(0:MC.anchor, each=(p-MC.anchor+1)),
                        add=rep(0:(p-MC.anchor), (MC.anchor+1)))
  Hplot.info=left_join(Hplot.info,Summary1,by=c("lose","add"))

  Hleft.hist=c()
  for(i in 0:(p-MC.anchor)){
    left.loc=which(Summary1[,2]==i)
    Hleft.hist[i+1]=sum(Summary1[left.loc,4])
  }
  H.left=data.frame(yloc=0:(p-MC.anchor),
                    xloc=Hleft.hist)
  Hup.hist=c()
  for(i in 0:MC.anchor){
    bottom.loc=which(Summary1[,1]==i)
    Hup.hist[i+1]=sum(Summary1[bottom.loc,4])
  }

  H.up=data.frame(xloc=0:MC.anchor,
                  yloc=Hup.hist)




  if(MC.anchor==p){
    Hline.x=c(rep(rep(0:MC.anchor, c(1,rep(2,MC.anchor-1),1)), (p-MC.anchor+1)), rep(0:MC.anchor, each=(p-MC.anchor)*2))
    Hline.y=rep(0, length(Hline.x))
  }else if(MC.anchor==0){
    Hline.y=c(rep(0:(p-MC.anchor), each=MC.anchor*2), rep(rep(0:(p-MC.anchor), c(1,rep(2,p-MC.anchor-1),1)), (MC.anchor+1)))
    Hline.x=rep(0, length(Hline.y))
  }else{
    Hline.x=c(rep(rep(0:MC.anchor, c(1,rep(2,MC.anchor-1),1)), (p-MC.anchor+1)), rep(0:MC.anchor, each=(p-MC.anchor)*2))
    Hline.y=c(rep(0:(p-MC.anchor), each=MC.anchor*2), rep(rep(0:(p-MC.anchor), c(1,rep(2,p-MC.anchor-1),1)), (MC.anchor+1)))
  }
  H.line.group=rep(1:(length(Hline.x)/2), each=2)
  H.line=data.frame(linex=Hline.x, liney=Hline.y, linegroup=H.line.group)


  if(is.null(circlesize)){
    circlesize = 10
  }
  if(is.null(linewidth)){
    linewidth = 1
  }
  if(is.null(fontsize)){
    fontsize = 1.5
  }


  p4<-
    ggplot(Hplot.info, aes(x=lose, y=add, colour=freq)) +
    geom_line(data = H.line, aes(linex, liney, group = linegroup),
              size = linewidth, alpha = 0.8, color = '#8BB6D6') +
    geom_point(shape=16, size=circlesize) +
    geom_text(aes(x = lose, y = add, label = freq), size=fontsize, color="white", na.rm=TRUE) +
    xlab( expression(paste(H^{"-"})) ) +
    ylab("") +
    theme_bw()+
    scale_color_viridis("frequency", option = "D", limits = c(0, 1)) +
    theme(text=element_text(size=8)
          ,axis.text.y=element_blank()
          ,plot.margin=unit(c(0.2,0.5,0.5,0), "cm"))


  p5<-
    ggplot(H.left, aes(x=yloc, y=xloc)) +
    geom_bar(stat = "identity") +
    xlab(expression(paste(H^{"+"}))) +
    ylab("frequency") +
    theme_bw() +
    coord_flip()+scale_y_reverse(limits=c(1,0), breaks=seq(1, 0, by=-0.5)) +
    theme(text=element_text(size=8)
          ,legend.position="none"
          ,axis.text.x = element_text(size = 6)
          ,plot.margin=unit(c(0.5,0.1,0.5,0.5), "cm")
    )


  p6<-
    ggplot(H.up, aes(x=xloc, y=yloc)) +
    geom_bar(stat = "identity") +
    xlab("") +
    ylab("") +
    theme_bw() +
    scale_y_continuous(limits=c(0,1), breaks=seq(0, 1, by=0.5)) +
    theme(text=element_text(size=8)
          ,axis.text.x=element_blank()
          ,axis.text.y = element_text(size = 6)
          ,legend.position="none"
          ,plot.margin=unit(c(0.5,0.5,-0.2,0.5), "cm")
    )

  ## Modify the user-specified value of ylim, default is c(0,10)
  if(is.null(ylim)){
    limitL=0
    limitH=10
    ylim=c(limitL,limitH)
  }else if(length(ylim)!=2){
    stop("ylim should be a vector with 2 elements or leave it NULL.")
  }else{
    limitL=ylim[1]
    limitH=ylim[2]
  }

  if(is.null(xlim)){
    limitLeft=0
    limitRight=10
    xlim=c(limitLeft,limitRight)
  }else if(length(xlim)!=2){
    stop("xlim should be a vector with 2 elements or leave it NULL.")
  }else{
    limitLeft=xlim[1]
    limitRight=xlim[2]
  }

  limits2 <- c(limitL-0.5, limitH+0.5)
  breaks2 <- seq(limits2[1]+0.5, limits2[2]-0.5, by=1)

  limits3 <- c(limitLeft-0.5, limitRight+0.5)
  breaks3 <- seq(limits3[1]+0.5, limits3[2]-0.5, by=1)

  # assign common axis to both plots
  p4.common.y <- p4 + coord_cartesian(ylim = limits2, xlim = limits3) +
    scale_y_continuous(breaks=breaks2) +
    scale_x_continuous(breaks=breaks3)
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
  return(list(Hplot.info=Hplot.info,
              Hplus.histogram=H.left,
              Hminus.histogram=H.up))
}
