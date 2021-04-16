#' VDSM-Scatter-heatmap-info
#'
#' Report VDSM-Scatter-heatmap-infomation
#' @param X A m*p matrix which contains m different p-dimensional models. All the elements are either 0 or 1.
#' @param f A vector with m elements which represent each model's frequency in X.
#' @param p The number of variate in the model
#' @param Anchor.model A vector containing p elements with either 1 or 0 value and must be found in X. Default is the model with the highest frequency.
#' @param Anchor.estimate An estimation for the anchor model
#' @return A list of information which helps to plot VDSM-Scatter-heatmap.
# @importFrom plyr count
# @importFrom dplyr arrange full_join left_join


VDSM_scatter_heat <- function(X,f,p,Anchor.estimate,Anchor.model=NULL){
  NewX=as.data.frame(CheckInput(X,f,p))
  Sorted.X=arrange(NewX,-frequency)

  m=nrow(Sorted.X)
  Sorted.X$xloc=rep(0,m)
  Sorted.X$yloc=rep(0,m)
  Table2=Sorted.X[,1:(p+1)]

  if(is.null(Anchor.model)){
    Anchor.model=as.numeric(Sorted.X[1,1:p])
  }else if(all(apply(abs(sweep(X,2,Anchor.model)),1,sum)!=0)==TRUE){
    Anchor.model=as.numeric(Sorted.X[1,1:p])
    print(paste("The anchor model user entered doesn't exist in X and has beeb automatically changed to the one with the highest frequency. User can also re-enter a new anchor model and try again."))
  }

  check1=which(Anchor.model!=0)
  check2=which(Anchor.estimate!=0)
  if(identical(check1,check2)==FALSE){
    stop("The entered Anchor.estimate is not the estimate for the anchor model. Please check and try again.")
  }
  MC.anchor=sum(Anchor.model)
  true.x.loc=which(Anchor.model!=0)

  W=Sorted.X
  if(MC.anchor==0){
    W[,(p+2)]=0
    W[,(p+3)]=rowSums(Sorted.X[,1:p])
    summary=W[,(p+1):(p+3)]
    df1=as.data.frame(summary)
    FreqModel=aggregate(frequency~xloc+yloc, data=df1, sum)
    x.comb=0
    combination=0
    hat.beta.true=0
  }else{
    # true.x.select=paste("x",true.x.loc,sep = "")
    # beta.hat.true=matrix(0,true.mode.freq,p)
    #
    # truemodel.mode.loc=which(apply(abs(sweep(M.models[,1:p],2,Anchor.model)),1,sum)==0)
    #
    # for(i in 1:true.mode.freq){
    #   l1=truemodel.mode.loc[i]
    #   X1=fixdata[[l1]][[2]]
    #   d=as.data.frame(X1)
    #   y.need=fixdata[[l1]][[1]]
    #   lm.need=lm(reformulate(true.x.select, "y.need"), data = d)
    #   beta.hat.true[i,true.x.loc]=lm.need$coefficients[-1]
    # }
    # hat.beta.true=apply(beta.hat.true,2,mean)
    hat.beta.true=Anchor.estimate

    tp1.loc=which(Anchor.model==1)
    tp2.loc=which(Anchor.model==0)
    SD.true.M.info=matrix(0,nrow(Table2),(p+1))
    SD.true.M.info[,(p+1)]=Table2[,(p+1)]

    for(i in 1:nrow(Table2)){
      true.part1=tp1.loc[which(Table2[i,tp1.loc]==0)]
      SD.true.M.info[i,true.part1]=hat.beta.true[true.part1]
      true.part2=tp2.loc[which(Table2[i,tp2.loc]==1)]
      SD.true.M.info[i,true.part2]=1
    }

    if(m==1){
      W[1,(p+2)]=0
      W[1,(p+3)]=0
      summary=W[,(p+1):(p+3)]
      FreqModel=as.data.frame(t(summary))
    }else{
      W[,(p+2)]=apply(abs(SD.true.M.info[,tp1.loc,drop=F]),1,sum)
      W[,(p+3)]=apply(SD.true.M.info[,tp2.loc,drop=F],1,sum)
      summary=W[,(p+1):(p+3)]
      df1=as.data.frame(summary)
      FreqModel=aggregate(frequency~xloc+yloc, data=df1, sum)
    }

    selectMatrix <- BitMatrix(length(true.x.loc))
    combination <-apply(selectMatrix, 1, function(x){ true.x.loc[which(x==1)] })
    x.comb=c()
    SD1.beta=abs(hat.beta.true)
    for (i in 1:(2^length(true.x.loc))){
      x.comb[i]=sum(SD1.beta[combination[[i]]])
    }
  }

  AT=data.frame(xloc=rep(x.comb,each=(p-length(true.x.loc)+1)),
                yloc=rep(0:(p-length(true.x.loc)), length(x.comb)),
                comb.loc=rep(1:length(x.comb), each=(p-length(true.x.loc)+1)))
  AT=left_join(AT,FreqModel,by=c("xloc","yloc"))

  left.x=c()
  for(i in 0:(p-length(true.x.loc))){
    left.loc=which(FreqModel[,2]==i)
    left.x[i+1]=sum(FreqModel[left.loc,3])
  }
  BT=data.frame(yloc=0:(p-length(true.x.loc)),
                xloc=left.x)
  bottom.y=c()
  for(i in 1:length(x.comb)){
    bottom.loc=which(FreqModel[,1]==x.comb[i])
    bottom.y[i]=sum(FreqModel[bottom.loc,3])
  }

  CT=data.frame(xloc=x.comb,
                yloc=bottom.y)
  return(list(AT, BT, CT, true.x.loc, x.comb, combination, hat.beta.true))

}


BitMatrix <- function(k){
  set <- 0:(2^k-1)
  rst <- matrix(0,ncol = k,nrow = 2^k)
  for (i in 1:k){
    rst[, i] = ifelse((set-rowSums(rst*rep(c(2^((k-1):0)), each=2^k)))/(2^(k-i))>=1, 1, 0)
  }
  return(rst)
}
