#' function for rank 2 inverse to update the inverse of X'X matrix for each coordinate exchange iteration
#' @param x the row vector with the new coordinate value
#' @param x_i the current row vector that needs to be exchanged
#' @param D Current (X'X)inv matrix
#' @return A matrix with the updated (X'X)inv if the exhange leads to a non singular X else a character string "Singular"
#'  @examples
#'   rank2_inverse(x,x_i,D)

rank2_inverse= function(x,x_i,D){
  f_x=matrix(x, ncol=1)
  f_x_i=matrix(x_i, ncol=1)
  U=cbind(f_x,f_x_i)
  V=rbind(t(f_x),-t(f_x_i))
  I_2=diag(1,nrow=2)
  r2=I_2+V%*%D%*%U
  #rank 2 update is performed only if the exhchange leads to a non--singular matrix
  #condition number of the matrix is used to check for singularity
  if(rcond(r2)>1e-8){
    #print(rcond(r2))
    D_new=D-D%*%U%*%solve(r2)%*%V%*%D
    if(rcond(D_new)>1e-15){
      #print(rcond(D_new))
      #return(round(D_new,8))
      return(D_new)
    }
    else{
      #print(paste("The matrix is singuar with condition number",rcond(D_new)))
      return("Singular")
    }
  }else{
    #print(paste("The matrix is singuar with condition number",rcond(r2)))
    return("Singular")
  }
}
