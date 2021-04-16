#' Group the models according to their Hamming distance and Model complexity to the anchor model
#'
#' Group the given models
#' @param X A m*p matrix which contains m different p-dimensional models. All the elements are either 0 or 1.
#' @param f A vector with m elements which represent each model's frequency in X.
#' @param p The number of variate in the model
#' @param Anchor.model A vector containing p elements with either 1 or 0 value and must be found in X. Default is the model with the highest frequency.
#' @return A summarized information of the grouped models.
# @importFrom plyr count
# @importFrom dplyr arrange full_join left_join

Groupinfo <- function(X,f,p,Anchor.model=NULL){
  NewX=as.data.frame(CheckInput(X,f,p))
  Sorted.X=arrange(NewX,-frequency)

  m=nrow(Sorted.X)
  Sorted.X$MC=rep(0,m)
  Sorted.X$HD=rep(0,m)

  if(is.null(Anchor.model)){
    Anchor.model=as.numeric(Sorted.X[1,1:p])
  }else if(all(apply(abs(sweep(X,2,Anchor.model)),1,sum)!=0)==TRUE){
    Anchor.model=as.numeric(Sorted.X[1,1:p])
    print(paste("The anchor model user entered doesn't exist in X and has beeb automatically changed to the one with the highest frequency. User can also re-enter a new anchor model and try again."))
  }

  Sorted.X[,(p+2)]=apply(Sorted.X[,1:p],1,sum)
  Sorted.X[,(p+3)]=apply(abs(sweep(Sorted.X[,1:p],2,Anchor.model,'-')),1,sum)
  summary=Sorted.X[,(p+1):(p+3)]
  colnames(summary)[1]="f"
  df1=as.data.frame(summary)
  NumModel=plyr::count(df1,c("MC","HD"))
  FreqModel=aggregate(f~MC+HD, data=df1, sum)
  Summary=full_join(NumModel,FreqModel,by=c("MC","HD"))
  colnames(Summary)=c("MC","HD","#exist","freq")

  MC.anchor=sum(Anchor.model)
  a=MC.anchor
  b=p-a

  Target=matrix(0,(p+1),3)
  Target[,2]=0:p
  colnames(Target)=c("#total","MC","HD")
  for(i in 1:(p+1)){
    if(Target[i,2] <= MC.anchor){
      Target[i,3]= MC.anchor-Target[i,2]
      c=Target[i,3]
      d=0
      Target[i,1]=choose(a,c)*choose(b,d)
    }else{
      Target[i,3]= Target[i,2] - MC.anchor
      c=0
      d=Target[i,3]
      Target[i,1]=choose(a,c)*choose(b,d)
    }
    ####################################################
    extracircle.number=min((a-c),(b-d))
    if(extracircle.number==0){
      Target=Target
    }else{
      start=Target[i,3]
      end=start+2*extracircle.number
      showlocation=seq(from=start, to=end, by=2)
      needlocation=showlocation[which(showlocation<=(p+1))]
      needlocation=needlocation[-1]
      addnum=length(needlocation)
      add.circle=matrix(0,addnum,3)
      add.circle[,2]=Target[i,2]
      add.circle[,3]=needlocation
      for(j in 1:addnum){
        add.circle[j,1]=choose(a,(c+j))*choose(b,(d+j))
      }
      Target=rbind(Target,add.circle)
    }
  }

  df2=as.data.frame(Target)
  final=left_join(df2,Summary,by=c("MC","HD"))

  select.x1=which(Anchor.model==1)
  select.x0=which(Anchor.model==0)
  col3=matrix(0,m,1)
  colnames(col3)="lose"
  col4=matrix(0,m,1)
  colnames(col4)="add"
  WU=cbind(Sorted.X[,1:(p+1)],col3,col4)
  if(m==1){
    WU[1,(p+2)]=0
    WU[1,(p+3)]=0

    T2=WU[,(p+1):(p+3)]
    colnames(T2)[1]="f"
    df3=as.data.frame(T2)
    NumModel1=plyr::count(df3,c("lose","add"))
    FreqModel1=as.data.frame(t(T2))
    Summary1=full_join(NumModel1,FreqModel1,by=c("lose","add"))
    colnames(Summary1)=c("lose","add","#exist","freq")
  }else{
    if(MC.anchor==1){
      WU[,(p+2)]=rep(1,MC.anchor)-WU[,select.x1]
      WU[,(p+3)]=apply(WU[,select.x0],1,sum)
    }else if(MC.anchor==0){
      WU[,(p+2)]=0
      WU[,(p+3)]=apply(WU[,select.x0],1,sum)
    }else if(MC.anchor==(p-1)){
      WU[,(p+2)]=apply((rep(1,MC.anchor)-WU[,select.x1]),1,sum)
      WU[,(p+3)]=WU[,select.x0]
    }else{
      WU[,(p+2)]=apply((rep(1,MC.anchor)-WU[,select.x1]),1,sum)
      WU[,(p+3)]=apply(WU[,select.x0],1,sum)
    }


    T2=WU[,(p+1):(p+3)]
    colnames(T2)[1]="f"
    df3=as.data.frame(T2)
    NumModel1=plyr::count(df3,c("lose","add"))
    FreqModel1=aggregate(f~lose+add, data=df3, sum)
    Summary1=full_join(NumModel1,FreqModel1,by=c("lose","add"))
    colnames(Summary1)=c("lose","add","#exist","freq")
  }

  return(list(final, FreqModel, MC.anchor, Summary1, Sorted.X))
}
