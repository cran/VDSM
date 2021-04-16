#' DSM_plot
#' plot the naive visualization of the distribution of selected model
#' @param X A m*p matrix which contains m different p-dimensional models. All the elements are either 0 or 1.
#' @param f A vector with m elements which represent each model's frequency in X.
#' @param p The number of variate in the model
#' @param Anchor.model A vector containing p elements with either 1 or 0 value and must be found in X. Default is the model with the highest frequency.
#' @param circlesize customize the size of the circle in the plot, default is 10.
#' @param linewidth Customize the width of the line in the plot, default is 1.
#' @param fontsize Customize the size of the font in the circles, default is 1.5.
#' @return A summarized information of the grouped models.
#' @importFrom dplyr arrange
#' @importFrom ggplot2 ggplot geom_line geom_text geom_point theme_bw theme geom_bar coord_flip scale_y_reverse scale_y_continuous coord_cartesian scale_x_continuous ggplot_gtable ggplot_build aes element_text element_blank xlab ylab theme_classic geom_tile scale_fill_continuous unit
#' @importFrom viridis scale_color_viridis
#' @export
#' @examples
#' data(exampleX)
#' X=exampleX
#' data(examplef)
#' f=examplef
#' p=8
#' DSM_example1 = DSM_plot(X,f,p)

DSM_plot<-function(X,f,p,Anchor.model=NULL,circlesize=NULL,linewidth=NULL,fontsize=NULL){
  NewX=as.data.frame(CheckInput(X,f,p))
  Sorted.X=arrange(NewX,-frequency)

  m=nrow(Sorted.X)
  Sorted.X$MC=rep(0,m)
  Sorted.X[,(p+2)]=apply(Sorted.X[,1:p],1,sum)
  if(is.null(Anchor.model)){
    Anchor.model=as.numeric(Sorted.X[1,1:p])
  }else if(all(apply(abs(sweep(X,2,Anchor.model)),1,sum)!=0)==TRUE){
    Anchor.model=as.numeric(Sorted.X[1,1:p])
    print(paste("The anchor model user entered doesn't exist in X and has beeb automatically changed to the one with the highest frequency. User can also re-enter a new anchor model and try again."))
  }

  Anchor.model=t(as.matrix(Anchor.model))
  Table5=as.matrix(arrange(Sorted.X, MC))
  possible.width=unique(Table5[,(p+2)])

  Col=matrix(0,nrow(Table5),1)
  colnames(Col)="color"
  Loc=matrix(0,nrow(Table5),1)
  colnames(Loc)="location"
  Table5=cbind(Table5,Col,Loc)

  #  thresh.num=rpercentile*M
  loc=1
  Table6=Table5
  Table6[,]=0

  ############################################3
  ## the full MSD plot
  Table7=Table5
  for(i in 1:length(possible.width)){
    wd=possible.width[i]
    aaa=which(Table7[,(p+2)]==wd)
    Table7[aaa,(p+4)]=1:length(aaa)
  }
  small.num=min(possible.width)
  large.num=max(possible.width)
  loc.true1=apply(Anchor.model,1,function(x) which(duplicated(rbind(x,Table7[,1:p])))-1)
  Table7[loc.true1,p+3]=3

  AT=as.data.frame(Table7)
  ATline.x=c()
  ATline.y=c()
  AT.line.group=c()
  for (i in small.num:(large.num-1)){
    start.loc=which(AT[,p+2]==i)
    end.loc=which(AT[,p+2]==(i+1))
    len1=length(start.loc)
    len2=length(which(AT[,p+2]==(i+1)))
    for (j in 1:len1){
      row1=which(AT[,p+2]==i & AT[,p+4]==j)
      for (k in 1:len2){
        row2=which(AT[,p+2]==i+1 & AT[,p+4]==k)
        if(all(AT[row2,(1:p)]-AT[row1,(1:p)]>=0)){
          vec1x=j
          vec1y=i
          vec2x=k
          vec2y=i+1
          ATline.x=c(ATline.x,vec1x,vec2x)
          ATline.y=c(ATline.y,vec1y,vec2y)
        }
      }
    }
  }


  AT.line.group=rep(1:(length(ATline.x)/2), each=2)
  AT.line=data.frame(linex=ATline.x, liney=ATline.y, linegroup=AT.line.group)

  if(is.null(circlesize)){
    circlesize = 13
  }
  if(is.null(linewidth)){
    linewidth = 1
  }
  if(is.null(fontsize)){
    fontsize = 3
  }

  plot1<-ggplot(AT, aes(x=location, y=MC, colour=frequency)) +
    geom_line(data = AT.line, aes(linex, liney, group = linegroup),
              size = linewidth, alpha = 0.8, color = '#8BB6D6') +
    geom_point(shape=16, size=circlesize) +
    geom_text(aes(x=location, y=MC, label=frequency), size=fontsize, color="white") +
    xlab("") +
    ylab("Model Complexity") +
    theme_bw()+
    scale_color_viridis("frequency", option = "D", limits = c(0, 1))
  print(plot1)
}
