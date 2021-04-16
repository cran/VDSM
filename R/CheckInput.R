#' Check if the input is valid or not
#'
#' Input a valid matrix
#' @param X A m*p matrix which each row represents one unique model with the elements either 0 or 1.
#' @param f A vector with m elements contain each model's frequency in X.
#' @param p The number of variate in the model
#' @return The standardized matrix
#' @export

CheckInput<-function(X,f,p){
  nrow.X=nrow(X)
  nrow.f=length(f)
  if(nrow.X!=nrow.f){
    stop("The number of unique models in X matrix is different from that in f vector.")
  }else{
    if(any(f>=1)==TRUE){
    f=f/sum(f)
    }
   Input=cbind(X,f)
   colnames(Input)[p+1]="frequency"
   colnames(Input)[1:p]=paste("x",1:p,sep="")
  }
  return(Input)
}

