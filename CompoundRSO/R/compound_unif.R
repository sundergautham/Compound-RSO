#' computes the W2 uniform discrepancy criteria update for each coordinate exchange
#' @param X settings matrix
#' @param i row being updated
#' @param j column being updated
#' @param xij_new new factor level at the i,j index
#' @param disc_matrix discrepancy matrix for W2 discrepancy measure
#' @param update  boolean value indicating whether to update or create a new discrepancy matrix
#' @return returns a discrepancy matrix
#' @description computes the W2 uniform discrepancy criteria update for each coordinate exchange
#' @examples
#' \dontrun{
#' discW2_update(X,i,j,xij_new,disc_matrix,update=FALSE)
#' }
#' @export

discW2_update=function(X,i,j,xij_new,disc_matrix,update=FALSE){
  if(!update){
    #change domain to [0,1]
    X=(X-(-1))/2
    disc_matrix=matrix(0,nrow=nrow(X),ncol=nrow(X))
    for(i in seq(nrow(X))){
      for(k in i:nrow(X)){
        temp=prod(3/2-abs(X[i,]-X[k,])+abs(X[i,]-X[k,])**2)
        disc_matrix[i,k]=temp
      }
    }
    disc_matrix[lower.tri(disc_matrix)]=t(disc_matrix)[lower.tri(disc_matrix)]
    return(disc_matrix)
  }
  if(update){
    X_new=X
    X_new[i,j]=xij_new
    #change domain to [0,1]
    X=(X-(-1))/2
    X_new=(X_new-(-1))/2
    xij_new=X_new[i,j]
    disc_matrix_new=disc_matrix
    for(col_up in 1:nrow(X)){
      up_value=(3/2-abs(xij_new-X_new[col_up,j])+abs(xij_new-X_new[col_up,j])**2)/(3/2-abs(X[i,j]-X[col_up,j])+abs(X[i,j]-X[col_up,j])**2)
      disc_matrix_new[i,col_up]=disc_matrix[i,col_up]*up_value
      disc_matrix_new[col_up,i]=disc_matrix[col_up,i]*up_value
    }
    return(disc_matrix_new)
  }
}


#' computes the M2 uniform discrepancy criteria update for each exchange
#' @param X settings matrix
#' @param i row being updated
#' @param j column being updated
#' @param xij_new new factor level at the i,j index
#' @param disc_matrix discrepancy matrix for M2 discrepancy measure
#' @param disc_vector discrepancy vector for M2 discrepancy measure
#' @param update  boolean value indicating whether to update or create a new discrepancy matrix
#' @param rescale boolean value indicating whether to rescale the design matrix to a (0,1) experimental region
#' @return returns a discrepancy matrix
#' @description computes the M2 uniform discrepancy criteria update for each exchange
#' @examples
#' \dontrun{
#' discM2_update(X,i,j,xij_new,disc_matrix,disc_vector,update=FALSE,rescale=TRUE)
#' }
#' @export

discM2_update=function(X,i,j,xij_new,disc_matrix,disc_vector,update=FALSE,rescale=TRUE){
  #update it to make use of the symmetric property
  if(!update){
    if(rescale){
      X=(X+1)/2
    }
    disc_matrix=matrix(0,nrow=nrow(X),ncol=nrow(X))
    disc_vector=rep(0,nrow(X))
    for(i in seq(nrow(X))){
      disc_vector[i]=prod(5/3-(1/4)*abs(X[i,]-0.5)-(1/4)*abs(X[i,]-0.5)**2)
      for(k in i:nrow(X)){
        temp=prod(15/8-(1/4)*abs(X[i,]-0.5)-(1/4)*abs(X[k,]-0.5)-0.75*abs(X[i,]-X[k,])+0.5*abs(X[i,]-X[k,])**2)
        disc_matrix[i,k]=temp
      }
    }
    disc_matrix[lower.tri(disc_matrix)]=t(disc_matrix)[lower.tri(disc_matrix)]
    return(list(disc_matrix,disc_vector))
  }
  if(update){
    X_new=X
    X_new[i,j]=xij_new
    if(rescale){
      X=(X+1)/2
      X_new=(X_new+1)/2
    }
    xij_new=X_new[i,j]
    disc_matrix_new=disc_matrix
    disc_vector_new=disc_vector
    #print(paste("All relevant intermediate values",disc_vector_new[i],(5/3-(1/4)*abs(xij_new-0.5)-(1/4)*abs(xij_new-0.5)**2),(5/3-(1/4)*abs(X[i,j]-0.5)-(1/4)*abs(X[i,j]-0.5)**2)))
    disc_vector_new[i]=disc_vector_new[i]*(5/3-(1/4)*abs(xij_new-0.5)-(1/4)*abs(xij_new-0.5)**2)/(5/3-(1/4)*abs(X[i,j]-0.5)-(1/4)*abs(X[i,j]-0.5)**2)
    for(col_up in 1:nrow(X)){
      up_value=(15/8-(1/4)*abs(xij_new-0.5)-(1/4)*abs(X_new[col_up,j]-0.5)-0.75*abs(xij_new-X_new[col_up,j])+0.5*abs(xij_new-X_new[col_up,j])**2)/(15/8-(1/4)*abs(X[i,j]-0.5)-(1/4)*abs(X[col_up,j]-0.5)-0.75*abs(X[i,j]-X[col_up,j])+0.5*abs(X[i,j]-X[col_up,j])**2)
      disc_matrix_new[i,col_up]=disc_matrix[i,col_up]*up_value
      disc_matrix_new[col_up,i]=disc_matrix[col_up,i]*up_value
    }
    return(list(disc_matrix_new,disc_vector_new))
  }
}

#Description:
##function for computing the compound uniform design using coordinate exchange

#Usage:
##coor_exch_compound_unif(x,x_matrix,C,model_order,C_dtype,freeze_rows=0,x_matrix_Dopt,x_unifopt,unif_crit="MD2",w,K=0,me_index_daug, qe_index_daug,two_fi_index_daug,telescoping=FALSE)

#Arguments:
## x: a starting matrix of factor levels in the design
## x_matrix: a starting design matrix of the specified order including the intercept
## C_dtype: data type of the factors in the design either "cont"(continuous) or "cat"(categorical)
## C: list of candidate points for all factors in the design
## model_order: string indicating linear,quadratic (Main Effects+Quadratic Effects), 2FI (Main Effects + TwoFactor Interactions),full, non_standard for a user defined order
## freeze_rows: the rows of the matrix to freeze while augmenting
## x_matrix_Dopt: the D-optimal design matrix
## D_opt_x: the factor levels matrix of the D-optimal design
## x_unifopt: the Uniform-optimal design matrix
## unif_crit: uniform criterion one of "MD2" (default) or "WD2"
## w: weight for computing the compound optimal design
## K: a diagonal matrix of prior variance for the Bayesian D optimal design
## me_index_daug: vector of index of Main Effects to be included in the D-optimal design for model_order "non_standard"
## qe_index_daug: vector index of Quadratic Effects to be included in the D-optimal design for model_order "non_standard"
## two_fi_index_daug: list of index of Two Factor Interactions to be included in the D-optimal design for model_order "non_standard"
## telescoping: a boolean value indicating whether telecoping (fine tuning) the current design or not

#Value: return a list of x_matrix,prev,x_best,deff_best,Unif_best
## x_matrix: the design matrix for a given model order
## prev: the compound optimal criteria
## x_best: matrix of factor levels for the design
## deff_best: the d-efficiency criteria
## Unif_best: the uniform efficiency criteria

# coor_exch_compound_unif=function(x,x_matrix,C,model_order,
#                                  C_dtype, freeze_rows=0,
#                                  x_matrix_Dopt,
#                                  x_unifopt,unif_crit="MD2",w,K=0,
#                                  me_index_daug, qe_index_daug,two_fi_index_daug,
#                                  telescoping=FALSE){
#   print("Executing R function")
#   print(paste("freeze rows:",freeze_rows))
#   x_start=x
#   p=ncol(x_matrix)
#   p_x=ncol(x_unifopt)
#   #p_iopt=ncol(x_matrix_Iopt)
#   #computes the inverse of X'x
#   D= solve(t(x_matrix)%*%x_matrix+K)
#   n=nrow(x_matrix)
#   #saves a copy of D
#   D_local=D
#   dopt=det(t(x_matrix_Dopt)%*%x_matrix_Dopt + K)
#   if(unif_crit=="WD2"){
#     unifopt=discrepancyCriteria(x_unifopt)$DisW2
#   }else if(unif_crit=="MD2"){
#     unifopt=discM2_update(X=x_unifopt,update=FALSE,rescale=FALSE)
#     disc_matrix=unifopt[[1]]
#     disc_vector=unifopt[[2]]
#     unifopt=sqrt((19/12)**(ncol(x_unifopt))-(2/n)*sum(disc_vector)+(1/n**2)*sum(disc_matrix))
#   }else{
#     stop("Error: Discrepancy criteria should be one of WD2 or MD2")
#   }
#   x_matrix_local=x_matrix
#
#   #epsilon for monitoring convergence
#   epsilon=1
#
#   #Current L-optimality score
#   M=det(t(x_matrix_local)%*%x_matrix_local + K)
#   D_eff=(M/dopt)**(1/p)
#   if(unif_crit=="WD2"){
#     disc_matrix=discW2_update(X=x_start,update=FALSE)
#     disc_val=-(4/3)**(ncol(x_start))+(1/nrow(x_start)**2)*sum(disc_matrix)
#     unif_disc=sqrt(disc_val)
#   }else if(unif_crit=="MD2"){
#     unif_disc=discM2_update(X=x_start,update=FALSE)
#     disc_matrix=unif_disc[[1]]
#     disc_vector=unif_disc[[2]]
#     unif_disc=sqrt((19/12)**(ncol(x_start))-(2/n)*sum(disc_vector)+(1/n**2)*sum(disc_matrix))
#   }else{
#     stop("Error: Discrepancy criteria should be one of WD2 or MD2")
#   }
#   #print(paste("Unif val from update",unif_disc))
#   #unif_disc=discrepancyCriteria(x_start)$DisW2
#   #print(paste("Unif desc from package",unif_disc))
#   Unif_eff=(unifopt/unif_disc)
#   current=w*D_eff+(1-w)*Unif_eff
#   print(paste("Ueff", Unif_eff))
#   print(paste("D eff",D_eff))
#   print(paste("starting objective:",current))
#   count_convergence=0
#   while(epsilon>1e-3){
#     count_convergence=count_convergence+1
#     # at the end of one complete exhchange of x_matrix_local, the updated  matrix is stored
#     x_matrix=x_matrix_local
#     prev=current
#     Unif_best=Unif_eff
#     deff_best=D_eff
#     x_best=x
#     for(i in (freeze_rows+1):n){
#       #print(paste("exchange of row number:", i))
#       for(j in 1:length(C)){
#         #print(paste("exchange of col number:", j))
#         #det_val=c()
#         #current value of trace
#         #I_val_current=sum(diag(D_local%*%W))
#         #D_val_current=M
#         if(telescoping){
#           point=x_start[i,j]
#           candidate_points=seq(max(-1,point-0.09),min(1,point+0.09),0.01)
#         }else{
#           candidate_points=C[[j]]
#         }
#         det_val=vector("numeric",length(candidate_points))
#         for(k in 1:length(candidate_points)){
#           #print(paste("Coordinate Exchange for row",i,"column",j,"candidate_point",k))
#           x_unif_local=x
#           x_unif_local[i,j]=candidate_points[k]
#           #disc_matrix_local=disc_matrix
#           x_local=x[i,]
#           #point exchange
#           x_local[j]=candidate_points[k]
#           x_expl=f_x(x_local,model_order="non_standard",C_dtype,
#                      me_index = me_index_daug,qe_index = qe_index_daug,
#                      two_fi_index = two_fi_index_daug)
#           D_new=rank2_inverse(x_expl,x_matrix_local[i,],D_local)
#           #print(D_new)
#           #check if singular
#           if(D_new[1]!="Singular"){
#             #M_local=det(t(x_new)%*%x_new)#M*delta_D(x_expl,x_matrix_local[i,],D_local)
#             print(paste("Value of M",M,delta_D(x_expl,x_matrix_local[i,],D_local)))
#             M_local=M*delta_D(x_expl,x_matrix_local[i,],D_local)
#             if(unif_crit=="WD2"){
#               disc_matrix_local=discW2_update(X=x,i=i,j=j,xij_new=candidate_points[k],disc_matrix=disc_matrix,update=TRUE)
#               disc_val_updated=-(4/3)**(ncol(x))+(1/nrow(x)**2)*sum(disc_matrix_local)
#               unif_disc_local=sqrt(disc_val_updated)
#             }else if(unif_crit=="MD2"){
#               unif_disc_local=discM2_update(X=x,i=i,j=j,xij_new=candidate_points[k],disc_matrix=disc_matrix,disc_vector=disc_vector,update=TRUE)
#               disc_matrix_local=unif_disc_local[[1]]
#               disc_vector_local=unif_disc_local[[2]]
#               unif_disc_local=sqrt((19/12)**(ncol(x))-(2/n)*sum(disc_vector_local)+(1/n**2)*sum(disc_matrix_local))
#             }else{
#               stop("Error: Discrepancy criteria should be one of WD2 or MD2")
#             }
#             #print(paste("Unif val from update",unif_disc_local))
#             #unif_disc_local_=discrepancyCriteria(x_unif_local)$DisW2
#             #print(paste("Unif desc from package",unif_disc_local_))
#             #if(!isTRUE(all.equal(unif_disc_local,unif_disc_local_,tolerance = 1e-6))){
#             #    print(all.equal(unif_disc_local,unif_disc_local_))
#             #    disc_matrix=discW2_update(X=x_unif_local,update=FALSE)
#             #    disc_val=-(4/3)**(ncol(x_unif_local))+(1/nrow(x_unif_local)**2)*sum(disc_matrix)
#             #    print(paste("Discrepancy from method 3:",sqrt(disc_val)))
#             #    stop("The two discrepancies don't match")
#             #}
#             #print(paste("D_eff and Unif_eff",(M_local/dopt)**(1/p),unifopt/unif_disc_local,unif_disc_local,M_local))
#             if(M_local>0){
#               local_val=w*((M_local/dopt)**(1/p))+(1-w)*(unifopt/unif_disc_local)
#               #det_val=c(det_val,local_val)
#               det_val[k]=local_val
#             }else{
#               #det_val=c(det_val,NaN)
#               det_val[k]=NaN
#             }
#           }else{
#             print("Singular hence adding NaN to det val")
#             #det_val=c(det_val,NaN)
#             det_val[k]=NaN
#           }
#           #print(paste("det_val",det_val))
#         }
#         #det_val=round(det_val,8)
#
#         #print("Done with det_val chekc")
#         if(any(is.na(det_val))){
#           print(paste("change in objective function at iteration",count_convergence))
#           print(det_val)
#         }
#         #det_val=ifelse(abs(det_val)==Inf,NaN,det_val)
#         #print(det_val)
#         #if(!any(is.na(det_val)) & !any(abs(det_val)==Inf)){
#         #if all of de_val is not NA
#         if(!all(is.na(det_val))){
#           #print("Starting the det_val chekc")
#           if(any(det_val>1.5,na.rm=TRUE)){
#             stop("Error: Efficiency is greater than 1, use a higher value of K")
#           }
#           ind=which(det_val==max(det_val,na.rm=TRUE))[1]
#           #print(paste("max=",ind))
#           #make the swap only if the delta is greater than 0 (there's an improvement)
#           if(det_val[ind]-current>0){
#             if(unif_crit=="WD2"){
#               disc_matrix=discW2_update(X=x,i=i,j=j,xij_new=candidate_points[ind],disc_matrix=disc_matrix,update=TRUE)
#               disc_val_updated=-(4/3)**(ncol(x))+(1/nrow(x)**2)*sum(disc_matrix)
#               Unif_val=sqrt(disc_val_updated)
#             }else if(unif_crit=="MD2"){
#               Unif_val=discM2_update(X=x,i=i,j=j,xij_new=candidate_points[ind],disc_matrix=disc_matrix,disc_vector=disc_vector,update=TRUE)
#               disc_matrix=Unif_val[[1]]
#               disc_vector=Unif_val[[2]]
#               Unif_val=sqrt((19/12)**(ncol(x))-(2/n)*sum(disc_vector)+(1/n**2)*sum(disc_matrix))
#             }else{
#               stop("Error: Discrepancy criteria should be one of WD2 or MD2")
#             }
#             x[i,j]=candidate_points[ind]
#             x_expl=f_x(x[i,],model_order="non_standard",C_dtype,
#                        me_index = me_index_daug,qe_index = qe_index_daug,
#                        two_fi_index = two_fi_index_daug)
#             M=M*delta_D(x_expl,x_matrix_local[i,],D_local)
#             D_local=rank2_inverse(x_expl,x_matrix_local[i,],D_local)
#             x_matrix_local[i,]=x_expl
#             #print(paste("Unif val from update",Unif_val))
#             #Unif_val=discrepancyCriteria(x)$DisW2
#             #print(paste("Unif desc from package",Unif_val))
#             Unif_eff=(unifopt/Unif_val)
#             D_eff=(M/dopt)**(1/p)
#             #print(paste("D efficiency",D_eff))
#             #print(paste("Uniform Efficiency",Unif_eff))
#             current=w*(D_eff)+(1-w)*Unif_eff
#             #print("############################")
#           }
#           #else{
#           #  print("No improvement in current coordinate update")
#           #}
#           #print("x matrix after interchange")
#           #print(x_matrix_local)
#           #print(paste("det of inf matrix after exchange:",determinant(t(x_matrix_local)%*%x_matrix_local)$modulus[1]))
#           #print("___________________________________")
#         }else{
#           #return(list(x_matrix,100000,x_best))
#           print(x_matrix_local)
#           print(j)
#           print(i)
#         }
#       }
#       #print(paste("det of inf matrix:",determinant(t(x_matrix_local)%*%x_matrix_local)$modulus[1]))
#       #current=sum(diag(D_local%*%W))
#     }
#     print(paste("Ueff", Unif_eff))
#     print(paste("D eff",D_eff))
#     print(paste("updated objective:",current))
#     #calculates epsilon
#     epsilon=current-prev
#     #print(paste("prev=",prev))
#     print(epsilon[1])
#     # if(current>1){
#     #   print(w)
#     #   print((M_local/dopt)**(1/p))
#     #   print((iopt/I_val))
#     # }
#     #D_eff=det(t(x_matrix_local)%*%x_matrix_local)/dopt
#     #print(D_eff)
#   }
#   print(paste("Number of iterations to convergence is:",count_convergence))
#   return(list(x_matrix,prev,x_best,deff_best,Unif_best))
# }


#Usage:
##
#Arguments:
## x: a starting matrix of factor levels in the design
## x_matrix: a starting design matrix of the specified order including the intercept
## C_dtype: data type of the factors in the design either "cont"(continuous) or "cat"(categorical)
## C: list of candidate points for all factors in the design
## model_order: string indicating linear,quadratic (Main Effects+Quadratic Effects), 2FI (Main Effects + TwoFactor Interactions),full, non_standard for a user defined order
## freeze_rows: the rows of the matrix to freeze while augmenting
## x_matrix_Dopt: the D-optimal design matrix
## D_opt_x: the factor levels matrix of the D-optimal design
## x_unifopt: the Uniform-optimal design matrix
## unif_crit: uniform criterion one of "MD2" (default) or "WD2"
## w: weight for computing the compound optimal design
## K: a diagonal matrix of prior variance for the Bayesian D optimal design
## me_index_daug: vector of index of Main Effects to be included in the D-optimal design for model_order "non_standard"
## qe_index_daug: vector index of Quadratic Effects to be included in the D-optimal design for model_order "non_standard"
## two_fi_index_daug: list of index of Two Factor Interactions to be included in the D-optimal design for model_order "non_standard"
## telescoping: a boolean value indicating whether telecoping (fine tuning) the current design or not

#Value: return a list of x_matrix,prev,x_best,deff_best,Unif_best
## x_matrix: the design matrix for a given model order
## prev: the compound optimal criteria
## x_best: matrix of factor levels for the design
## deff_best: the d-efficiency criteria
## Unif_best: the uniform efficiency criteria

#' function for computing the compound uniform design using coordinate exchange with foldover constraint
#' @param x settings matrix
#' @param x_matrix design matrix of the specified order including the intercept
#' @param C list of candidate space for each experimental factor in the design
#' @param model_order a string indicating the order of the regression model. "linear","quadratic" (Main Effects+Quadratic Effects), "2FI" (Main Effects + TwoFactor Interactions),"full", "non_standard"- for a user defined model order order
#' @param C_dtype a string vector indicating the data type of the experimental factors in the design either "cont"(continuous) or "cat"(categorical)
#' @param freeze_rows the rows of the matrix to freeze while augmenting, 0 when type is "non_augment"
#' @param x_matrix_Dopt the D-optimal design matrix
#' @param x_unifopt the Uniform-optimal design matrix
#' @param unif_crit string indicating the uniform optimal criterion. One of "MD2" (default) or "WD2"
#' @param w weight for computing the compound optimal design
#' @param K a diagonal matrix of prior variance for the Bayesian D optimal design
#' @param me_index_daug vector of index of Main Effects to be included in the D-optimal design for model_order "non_standard"
#' @param qe_index_daug vector index of Quadratic Effects to be included in the D-optimal design for model_order "non_standard"
#' @param two_fi_index_daug list of index of Two Factor Interactions to be included in the D-optimal design for model_order "non_standard"
#' @param telescoping a boolean value taking TRUE or FALSE. Telescoping reduces the candidate factors the local region of the optimum design
#' @description function for computing the compound uniform design using coordinate exchange with foldover constraint
#' @return a list of
#'  \itemize{
#'        \item Compound_unif_design- Comound Uniform optimal design matrix of the specified order including the intercept
#'        \item Compound_unif_obj- The value of the Comound Uniform Optimal objective value
#'        \item Compound_unif_x- a matrix with factor levels of the Comound Uniform optimal design
#'        \item Compound_Deff- Bayesian D effciency of the compound uniform design
#'        \item Compound_Ueff- Uniform Effciency of the compound uniform design
#'        }
#'  @examples
#'  \dontrun{
#'   coor_exch_compound_unif_foldover(x,x_matrix,C,model_order,C_dtype,freeze_rows=0,x_matrix_Dopt,x_unifopt,unif_crit="MD2",w,K=0,me_index_daug, qe_index_daug,two_fi_index_daug,telescoping=FALSE)
#'  }
#'  @export

coor_exch_compound_unif_foldover=function(x,x_matrix,C,model_order,
                                          C_dtype, freeze_rows=0,
                                          x_matrix_Dopt,
                                          x_unifopt,unif_crit="MD2",w,K=0,
                                          me_index_daug, qe_index_daug,two_fi_index_daug,
                                          telescoping=FALSE){
  #print(freeze_rows)
  x_start=x
  p=ncol(x_matrix)
  #p_iopt=ncol(x_matrix_Iopt)
  #computes the inverse of X'x
  D= solve(t(x_matrix)%*%x_matrix+K)
  n=nrow(x_matrix)
  #saves a copy of D
  D_local=D
  dopt=det(t(x_matrix_Dopt)%*%x_matrix_Dopt + K)
  if(unif_crit=="WD2"){
    unifopt=DiceDesign::discrepancyCriteria(x_unifopt)$DisW2
  }else if(unif_crit=="MD2"){
    unifopt=discM2_update(X=x_unifopt,update=FALSE,rescale=FALSE)
    disc_matrix=unifopt[[1]]
    disc_vector=unifopt[[2]]
    unifopt=sqrt((19/12)**(ncol(x_unifopt))-(2/n)*sum(disc_vector)+(1/n**2)*sum(disc_matrix))
  }else{
    stop("Error: Discrepancy criteria should be one of WD2 or MD2")
  }
  x_matrix_local=x_matrix

  #epsilon for monitoring convergence
  epsilon=1

  #Current L-optimality score
  M=det(t(x_matrix_local)%*%x_matrix_local + K)
  D_eff=(M/dopt)**(1/p)
  if(unif_crit=="WD2"){
    disc_matrix=discW2_update(X=x_start,update=FALSE)
    disc_val=-(4/3)**(ncol(x_start))+(1/nrow(x_start)**2)*sum(disc_matrix)
    unif_disc=sqrt(disc_val)
  }else if(unif_crit=="MD2"){
    unif_disc=discM2_update(X=x_start,update=FALSE)
    disc_matrix=unif_disc[[1]]
    disc_vector=unif_disc[[2]]
    unif_disc=sqrt((19/12)**(ncol(x_start))-(2/n)*sum(disc_vector)+(1/n**2)*sum(disc_matrix))
  }else{
    stop("Error: Discrepancy criteria should be one of WD2 or MD2")
  }
  #print(paste("Unif val from update",unif_disc))
  #unif_disc=discrepancyCriteria(x_start)$DisW2
  #print(paste("Unif desc from package",unif_disc))
  Unif_eff=(unifopt/unif_disc)
  current=w*D_eff+(1-w)*Unif_eff
  print(paste("Ueff", Unif_eff))
  print(paste("D eff",D_eff))
  print(paste("starting objective:",current))
  count_convergence=0
  while(epsilon>1e-3){
    count_convergence=count_convergence+1
    # at the end of one complete exhchange of x_matrix_local, the updated  matrix is stored
    x_matrix=x_matrix_local
    prev=current
    Unif_best=Unif_eff
    deff_best=D_eff
    x_best=x
    for(i in seq((freeze_rows+1),n,2)){
      #print(paste("exchange of row number:", i))
      for(j in 1:length(C)){
        #print(paste("exchange of col number:", j))
        det_val=c()
        #current value of trace
        #I_val_current=sum(diag(D_local%*%W))
        #D_val_current=M
        if(telescoping){
          point=x_start[i,j]
          candidate_points=seq(max(-1,point-0.09),min(1,point+0.09),0.01)
        }else{
          candidate_points=C[[j]]
        }
        for(k in candidate_points){
          print(paste("Coordinate Exchange for row",i,"column",j,"candidate_point",k))
          x_unif_local=x
          x_unif_local[i,j]=k
          #disc_matrix_local=disc_matrix
          x_local=x[i,]
          #point exchange
          x_local[j]=k
          x_expl=f_x(x_local,model_order="non_standard",C_dtype,
                     me_index = me_index_daug,qe_index = qe_index_daug,
                     two_fi_index = two_fi_index_daug)
          D_new=rank2_inverse(x_expl,x_matrix_local[i,],D_local)
          #check if singular
          if(D_new[1]!="Singular"){
            #M_local=det(t(x_new)%*%x_new)#M*delta_D(x_expl,x_matrix_local[i,],D_local)
            M_local=M*delta_D(x_expl,x_matrix_local[i,],D_local)
            if(unif_crit=="WD2"){
              disc_matrix_local=discW2_update(X=x,i=i,j=j,xij_new=k,disc_matrix=disc_matrix,update=TRUE)
              disc_val_updated=-(4/3)**(ncol(x))+(1/nrow(x)**2)*sum(disc_matrix_local)
              unif_disc_local=sqrt(disc_val_updated)
            }else if(unif_crit=="MD2"){
              unif_disc_local=discM2_update(X=x,i=i,j=j,xij_new=k,disc_matrix=disc_matrix,disc_vector=disc_vector,update=TRUE)
              disc_matrix_local=unif_disc_local[[1]]
              disc_vector_local=unif_disc_local[[2]]
              unif_disc_local=sqrt((19/12)**(ncol(x))-(2/n)*sum(disc_vector_local)+(1/n**2)*sum(disc_matrix_local))
            }else{
              stop("Error: Discrepancy criteria should be one of WD2 or MD2")
            }
            #print(paste("Unif val from update",unif_disc_local))
            #unif_disc_local_=discrepancyCriteria(x_unif_local)$DisW2
            #print(paste("Unif desc from package",unif_disc_local_))
            #if(!isTRUE(all.equal(unif_disc_local,unif_disc_local_,tolerance = 1e-6))){
            #    print(all.equal(unif_disc_local,unif_disc_local_))
            #    disc_matrix=discW2_update(X=x_unif_local,update=FALSE)
            #    disc_val=-(4/3)**(ncol(x_unif_local))+(1/nrow(x_unif_local)**2)*sum(disc_matrix)
            #    print(paste("Discrepancy from method 3:",sqrt(disc_val)))
            #    stop("The two discrepancies don't match")
            #}
            #x_unif_local[(i+1),j]=-k
            #disc_matrix_local=disc_matrix
            x_local=x[(i+1),]
            #point exchange
            x_local[j]=-k
            x_expl=f_x(x_local,model_order="non_standard",C_dtype,
                       me_index = me_index_daug,qe_index = qe_index_daug,
                       two_fi_index = two_fi_index_daug)
            D_new=rank2_inverse(x_expl,x_matrix_local[(i+1),],D_new)
            if(D_new[1]!="Singular"){
              #print("Non Singular")
              #M_local=det(t(x_new)%*%x_new)#M*delta_D(x_expl,x_matrix_local[i,],D_local)
              M_local=M_local*delta_D(x_expl,x_matrix_local[(i+1),],D_new)
              if(unif_crit=="WD2"){
                disc_matrix_local=discW2_update(X=x_unif_local,i=(i+1),j=j,xij_new=-k,disc_matrix=disc_matrix_local,update=TRUE)
                disc_val_updated=-(4/3)**(ncol(x))+(1/nrow(x)**2)*sum(disc_matrix_local)
                unif_disc_local=sqrt(disc_val_updated)
              }else if(unif_crit=="MD2"){
                unif_disc_local=discM2_update(X=x_unif_local,i=(i+1),j=j,xij_new=-k,disc_matrix=disc_matrix_local,disc_vector=disc_vector_local,update=TRUE)
                disc_matrix_local=unif_disc_local[[1]]
                disc_vector_local=unif_disc_local[[2]]
                unif_disc_local=sqrt((19/12)**(ncol(x))-(2/n)*sum(disc_vector_local)+(1/n**2)*sum(disc_matrix_local))
              }else{
                stop("Error: Discrepancy criteria should be one of WD2 or MD2")
              }
              local_val=w*((M_local/dopt)**(1/p))+(1-w)*(unifopt/unif_disc_local)
              #print(M_local)
              #print((M_local/dopt)**(1/p))
              #print(unif_disc_local)
              det_val=c(det_val,local_val)
            }else{
              #print("singular")
              det_val=c(det_val,NaN)
            }
          }
        }
        #det_val=round(det_val,8)
        if(any(is.na(det_val))){
          print(paste("change in objective function at iteration",i,j))
          print(det_val)
        }
        #det_val=ifelse(abs(det_val)==Inf,NaN,det_val)
        #print(det_val)
        #if(!any(is.na(det_val)) & !any(abs(det_val)==Inf)){
        #if all of de_val is not NA
        if(!all(is.na(det_val))){
          ind=which(det_val==max(det_val,na.rm=TRUE))[1]
          #print(paste("max=",ind))
          #make the swap only if the delta is greater than 0
          if(det_val[ind]-current>0){
            if(unif_crit=="WD2"){
              x_unif_local=x
              disc_matrix=discW2_update(X=x,i=i,j=j,xij_new=candidate_points[ind],disc_matrix=disc_matrix,update=TRUE)
              x_unif_local[i,j]=candidate_points[ind]
              disc_matrix=discW2_update(X=x_unif_local,i=(i+1),j=j,xij_new=-candidate_points[ind],disc_matrix=disc_matrix,update=TRUE)
              disc_val_updated=-(4/3)**(ncol(x))+(1/nrow(x)**2)*sum(disc_matrix)
              Unif_val=sqrt(disc_val_updated)
            }else if(unif_crit=="MD2"){
              x_unif_local=x
              Unif_val=discM2_update(X=x,i=i,j=j,xij_new=candidate_points[ind],disc_matrix=disc_matrix,disc_vector=disc_vector,update=TRUE)
              disc_matrix=Unif_val[[1]]
              disc_vector=Unif_val[[2]]
              x_unif_local[i,j]=candidate_points[ind]
              Unif_val=discM2_update(X=x_unif_local,i=(i+1),j=j,xij_new=-candidate_points[ind],disc_matrix=disc_matrix,disc_vector=disc_vector,update=TRUE)
              disc_matrix=Unif_val[[1]]
              disc_vector=Unif_val[[2]]
              Unif_val=sqrt((19/12)**(ncol(x))-(2/n)*sum(disc_vector)+(1/n**2)*sum(disc_matrix))
            }else{
              stop("Error: Discrepancy criteria should be one of WD2 or MD2")
            }
            x[i,j]=candidate_points[ind]
            M=M*delta_D(f_x(x[i,],model_order="non_standard",C_dtype,
                            me_index = me_index_daug,qe_index = qe_index_daug,
                            two_fi_index = two_fi_index_daug),x_matrix_local[i,],D_local)
            D_local=rank2_inverse(f_x(x[i,],model_order="non_standard",C_dtype,
                                      me_index = me_index_daug,qe_index = qe_index_daug,
                                      two_fi_index = two_fi_index_daug),x_matrix_local[i,],D_local)
            x_matrix_local[i,]=f_x(x[i,],model_order="non_standard",C_dtype,
                                   me_index = me_index_daug,qe_index = qe_index_daug,
                                   two_fi_index = two_fi_index_daug)
            x[(i+1),j]=-candidate_points[ind]
            M=M*delta_D(f_x(x[(i+1),],model_order="non_standard",C_dtype,
                            me_index = me_index_daug,qe_index = qe_index_daug,
                            two_fi_index = two_fi_index_daug),x_matrix_local[(i+1),],D_local)
            D_local=rank2_inverse(f_x(x[(i+1),],model_order="non_standard",C_dtype,
                                      me_index = me_index_daug,qe_index = qe_index_daug,
                                      two_fi_index = two_fi_index_daug),x_matrix_local[(i+1),],D_local)
            x_matrix_local[(i+1),]=f_x(x[(i+1),],model_order="non_standard",C_dtype,
                                       me_index = me_index_daug,qe_index = qe_index_daug,
                                       two_fi_index = two_fi_index_daug)
            #print(paste("Unif val from update",Unif_val))
            #Unif_val=discrepancyCriteria(x)$DisW2
            #print(paste("Unif desc from package",Unif_val))
            Unif_eff=(unifopt/Unif_val)
            D_eff=(M/dopt)**(1/p)
            #print(paste("D efficiency",D_eff))
            #print(paste("Uniform Efficiency",Unif_eff))
            current=w*(D_eff)+(1-w)*Unif_eff
            #print("############################")
          }
          #else{
          #  print("No improvement in current coordinate update")
          #}
          #print("x matrix after interchange")
          #print(x_matrix_local)
          #print(paste("det of inf matrix after exchange:",determinant(t(x_matrix_local)%*%x_matrix_local)$modulus[1]))
          #print("___________________________________")
        }else{
          #return(list(x_matrix,100000,x_best))
          print(x_matrix_local)
          print(j)
          print(i)
        }
      }
      #print(paste("det of inf matrix:",determinant(t(x_matrix_local)%*%x_matrix_local)$modulus[1]))
      #current=sum(diag(D_local%*%W))
    }
    print(paste("Ueff", Unif_eff))
    print(paste("D eff",D_eff))
    print(paste("updated objective:",current))
    #calculates epsilon
    epsilon=current-prev
    #print(paste("prev=",prev))
    print(epsilon[1])
    # if(current>1){
    #   print(w)
    #   print((M_local/dopt)**(1/p))
    #   print((iopt/I_val))
    # }
    #D_eff=det(t(x_matrix_local)%*%x_matrix_local)/dopt
    #print(D_eff)
  }
  print(paste("Number of iterations to convergence is:",count_convergence))
  return(list("Compound_unif_design"=x_matrix,"Compound_unif_obj"=prev,"Compound_unif_x"=x_best,"Compound_Deff"=deff_best,"Compound_Ueff"=Unif_best))
}

#' function for initiating individual starts of the compound uniform design using coordinate exchange.This is the parallel version and works only on linux machines.
#' @param i index of the random start
#' @param x_matrix_start the design matrix of the specified order including the intercept of the random starting
#' @param x_start the settings matrix of the random starting design
#' @param n total number of runs in the design (including the augmentations)
#' @param C list of candidate space for each experimental factor in the design
#' @param model_order a string indicating the order of the regression model. "linear","quadratic" (Main Effects+Quadratic Effects), "2FI" (Main Effects + TwoFactor Interactions),"full", "non_standard"- for a user defined model order order
#' @param C_dtype a string vector indicating the data type of the experimental factors in the design either "cont"(continuous) or "cat"(categorical)
#' @param type a string indicating the type of design augmentation. "non_augment" or "augment" (augment an existing design)
#' @param freeze_rows numeric value indicating the rows to not update
#' @param aug_type a string indicating whether the design has constraint or not. "foldover" or "not_foldover"
#' @param D_opt_design the D-optimal design matrix
#' @param D_opt_x the settings matrix of the D-optimal design
#' @param x_unifopt the settings matrix of the Uniform-optimal design
#' @param unif_crit a string indicating the discrepancy measure for uniform criterion. "MD2" or "WD2"
#' @param D_opt_type a string indicating whether the D-optimal design type is "Bayesian" or "Non_Bayesian"
#' @param w a numeric value of weight used to compute the compute the compound optimal design
#' @param me_index_daug vector of index of Main Effects to be included in the D-optimal design for model_order "non_standard"
#' @param qe_index_daug vector index of Quadratic Effects to be included in the D-optimal design for model_order "non_standard"
#' @param two_fi_index_daug list of index of Two Factor Interactions to be included in the D-optimal design for model_order "non_standard"
#' @param K a diagonal matrix specifying the prior variance in Bayesian D-optimal design set to 0.001 for potential terms and 0 for primary terms
#' @param use_cpp a boolean value indicating whether to use the CPP version of coordinate exchange
#' @param update_all_rows a boolean value indicating whether all the rows should be updated or only the augmented rows should be updated
#' @description function for initiating individual starts of the compound uniform design using coordinate exchange.This is the parallel version and works only on linux machines.
#' @return a list of
#'  \itemize{
#'        \item Compound_unif_design- Comound Uniform optimal design matrix of the specified order including the intercept
#'        \item Compound_unif_obj- The value of the Comound Uniform Optimal objective value
#'        \item Compound_unif_x- a matrix with factor levels of the Comound Uniform optimal design
#'        \item Compound_Deff- Bayesian D effciency of the compound uniform design
#'        \item Compound_Ueff- Uniform Effciency of the compound uniform design
#'        }
#'  @examples
#'  \dontrun{
#'   individual_starts_compound_unif(i,x_matrix_start,x_start,n,C_dtype,C,model_order, type,freeze_rows,D_opt_design,D_opt_x,x_unifopt,unif_crit,w,use_cpp=FALSE,K,D_opt_type="Bayesian",aug_type="not_foldover",me_index_daug, qe_index_daug,two_fi_index_daug,update_all_rows=FALSE)
#'    }
#'  @export
individual_starts_compound_unif=function(i,x_matrix_start,x_start,n,C_dtype,C,
                                         model_order, type,freeze_rows,
                                         D_opt_design,D_opt_x,x_unifopt,unif_crit,
                                         w,use_cpp=FALSE,K,D_opt_type="Bayesian",aug_type="not_foldover",
                                         me_index_daug, qe_index_daug,two_fi_index_daug,update_all_rows=FALSE){
  #if(i<2){
  #  print("Compound Optimal with Bayesian D-Optimal Design as starting design")
  #  x_matrix=D_opt_design
  #  x=D_opt_x
  #}else{
  #print("aug_type")
  #print(aug_type)
  singular=TRUE
  while(singular){
    if(aug_type=="not_foldover"){
      if(D_opt_type=="Bayesian"){
        x=gen_x(n=n-freeze_rows,C,model_order="non_standard",C_dtype,
                type="augment",
                me_index = me_index_daug,qe_index = qe_index_daug,
                two_fi_index = two_fi_index_daug,design="Bayesian")
      }else{
        x=gen_x(n=n-freeze_rows,C,model_order=model_order,C_dtype,
                type=type,
                me_index = me_index_daug,qe_index = qe_index_daug,
                two_fi_index = two_fi_index_daug)
      }
      x_matrix=x[[2]]
      #print(dim(x_matrix))
      #print(dim(x_matrix_start))
      x_matrix=rbind(x_matrix_start,x_matrix)
      x=x[[1]]
      if(class(x)!="matrix"){x=matrix(x)}
      #print(x)
      #print(x_matrix)
      x=rbind(x_start,x)
    }else if(aug_type=="foldover"){
      #print("Entered foldover condition")
      if(n%%2!=0){
        stop("Error: For a foldover design total number of runs should be a even")
      }
      if(D_opt_type=="Bayesian"){
        x=gen_x(n=(n-freeze_rows)/2,C,model_order="non_standard",C_dtype,
                type="augment",
                me_index = me_index_daug,qe_index = qe_index_daug,
                two_fi_index = two_fi_index_daug,design="Bayesian")
      }else{
        x=gen_x(n=n-freeze_rows,C,model_order=model_order,C_dtype,
                type=type,
                me_index = me_index_daug,qe_index = qe_index_daug,
                two_fi_index = two_fi_index_daug)
      }
      x_matrix=x[[2]]
      x_matrix_=rbind(x_matrix[1,],-x_matrix[1,])
      x=x[[1]]
      x_=rbind(x[1,],-x[1,])
      for(i in 1:nrow(x_matrix)){
        if(i>1){
          x_matrix_=rbind(x_matrix_,x_matrix[i,],-x_matrix[i,])
          x_=rbind(x_,x[i,],-x[i,])
        }
      }
      #print(dim(x_matrix))
      #print(dim(x_matrix_start))
      x_matrix=x_matrix_
      x=x_
      x_matrix=rbind(x_matrix_start,x_matrix)

      if(class(x)!="matrix"){x=matrix(x)}
      #print(x)
      #print(x_matrix)
      x=rbind(x_start,x)
    }else{
      stop("Aug tye should be one of foldover or not_foldover")
    }
    if(rcond(t(x_matrix)%*%x_matrix+K)>1e-15){
      singular=FALSE
    }else{
      print("Design Obtained is singular")
      print(rcond(t(x_matrix)%*%x_matrix+K))
    }
  }
  #print(x)
  #x_matrix=cnn_doe_mat
  #x=x_cnn_doe
  if(class(x)!="matrix"){
    x=matrix(x)
  }
  #}
  print(paste("Coordinate exchange for try",i))
  print(paste("Dim of x is",dim(x)))
  #####################################################
  if(update_all_rows){
    freeze_rows=0
  }
  if(use_cpp){
    me_index_daug_=me_index_daug-1
    qe_index_daug_=qe_index_daug-1
    two_fi_index_daug_=lapply(two_fi_index_daug,FUN=function(z){z-1})
  }else{
    me_index_daug_=me_index_daug
    qe_index_daug_=qe_index_daug
    two_fi_index_daug_=two_fi_index_daug
  }
  #print(me_index_daug_)
  #print(qe_index_daug_)
  #print(two_fi_index_daug_)
  if(aug_type=="not_foldover"){
    des=coor_exch_compound_unif(x=x,x_matrix=x_matrix,C=C,model_order =model_order,
                                C_dtype=C_dtype,x_matrix_Dopt = D_opt_design,
                                w=w, x_unifopt= x_unifopt,unif_crit=unif_crit,
                                freeze_rows = freeze_rows,
                                K=K,
                                me_index_daug=me_index_daug_,
                                qe_index_daug=qe_index_daug_,
                                two_fi_index_daug=two_fi_index_daug_,telescoping=FALSE)
    print("Starting telescoping")
    des=coor_exch_compound_unif(x=des[[3]],x_matrix=des[[1]],C=C,model_order =model_order,
                                C_dtype=C_dtype,x_matrix_Dopt = D_opt_design,
                                w=w, x_unifopt= x_unifopt,unif_crit=unif_crit,
                                freeze_rows = freeze_rows,
                                K=K,
                                me_index_daug=me_index_daug_,
                                qe_index_daug=qe_index_daug_,
                                two_fi_index_daug=two_fi_index_daug_,
                                telescoping=TRUE)
  }else if(aug_type=="foldover"){
    #print(freeze_rows)
    des=coor_exch_compound_unif_foldover(x=x,x_matrix=x_matrix,C=C,model_order =model_order,
                                         C_dtype=C_dtype,x_matrix_Dopt = D_opt_design,
                                         w=w, x_unifopt= x_unifopt,unif_crit=unif_crit,
                                         freeze_rows = freeze_rows,
                                         K=K,
                                         me_index_daug=me_index_daug_,
                                         qe_index_daug=qe_index_daug_,
                                         two_fi_index_daug=two_fi_index_daug_,telescoping=FALSE)
    print("Starting telescoping")
    des=coor_exch_compound_unif_foldover(x=des[[3]],x_matrix=des[[1]],C=C,model_order =model_order,
                                         C_dtype=C_dtype,x_matrix_Dopt = D_opt_design,
                                         w=w, x_unifopt= x_unifopt,unif_crit=unif_crit,
                                         freeze_rows = freeze_rows,
                                         K=K,
                                         me_index_daug=me_index_daug_,
                                         qe_index_daug=qe_index_daug_,
                                         two_fi_index_daug=two_fi_index_daug_,
                                         telescoping=TRUE)
  }
  return(des)
}

#' function for initiating multiple starts of an compound uniform design using coordinate exchange.This is the parallel version and works only on linux machines.
#' @param x_matrix_start the design matrix of the specified order including the intercept of the random starting
#' @param x_start the settings matrix of the random starting design
#' @param n total number of runs in the design
#' @param C list of candidate space for each experimental factor in the design
#' @param model_order a string indicating the order of the regression model. "linear","quadratic" (Main Effects+Quadratic Effects), "2FI" (Main Effects + TwoFactor Interactions),"full", "non_standard"- for a user defined model order order
#' @param C_dtype a string vector indicating the data type of the experimental factors in the design either "cont"(continuous) or "cat"(categorical)
#' @param type a string indicating the type of design augmentation. "non_augment" or "augment" (augment an existing deisgn)
#' @param aug_type a string indicating whether the design has constraint or not. "foldover" or "not_foldover"
#' @param D_opt_design the D-optimal design matrix
#' @param D_opt_x the settings matrix of the D-optimal design
#' @param x_unifopt the settings matrix of the Uniform-optimal design
#' @param unif_crit a string indicating the discrepancy measure for uniform criterion. "MD2" or "WD2"
#' @param D_opt_type a string indicating whether the D-optimal design type is "Bayesian" or "Non_Bayesian"
#' @param w_candidate a vector of w values to compute the compound uniform design for
#' @param use_w a boolean value indicating whether w_candidate should be used or the bisection search method for finding the compound optimal designs
#' @param me_index_daug vector of index of Main Effects to be included in the D-optimal design for model_order "non_standard"
#' @param qe_index_daug vector index of Quadratic Effects to be included in the D-optimal design for model_order "non_standard"
#' @param two_fi_index_daug list of index of Two Factor Interactions to be included in the D-optimal design for model_order "non_standard"
#' @param K a diagonal matrix specifying the prior variance in Bayesian D-optimal design set to 0.001 for potential terms and 0 for primary terms
#' @param use_cpp a boolean value indicating whether to use the CPP version of coordinate exchange
#' @param update_all_rows a boolean value indicating whether all the rows should be updated or only the augmented rows should be updated
#' @param n_starts number of random starts, default at 100
#' @param D_opt_tresh a numeric value providing the lower bound of the treshold in bisection search for the D-optimal criteria. Default value is 0.8
#' @param U_opt_tresh a numeric value providing the lower bound of the treshold in bisection search for the U-optimal criteria. Default value is 0.8
#' @param w_tresh the threshold for difference in weights w for the bisection search to converge. The default value is 0.02
#' @param rng seed for random number generation for the random design construction starts
#' @param num_cores number of cores for parallel processing. Default value is NA, in which case Sys.getenv("PBS_RESC_TOTAL_PROCS") is used to get the number of cores in a HPC. Otherwise, the user needs to provide a value. Use 1 for serial execution.
#' @description function for initiating multiple starts of an compound uniform design using coordinate exchange.This is the parallel version and works only on linux machines.
#' @return a list of compound optimal designs for individual weights w
#'  @examples
#'   construct_des_compound_unif(x_matrix_start,x_start,n,C_dtype,C,model_order, type= "augment",aug_type="not_foldover",D_opt_design,D_opt_x,x_unifopt,unif_crit,n_starts=100,D_opt_type="Bayesian",w_candidate=seq(0.05,0.95,0.05),use_w=FALSE,use_cpp=FALSE,K=0,me_index_daug, qe_index_daug,two_fi_index_daug,update_all_rows=FALSE)
#' @import doRNG
#' @export

construct_des_compound_unif=function(x_matrix_start,x_start,n,C_dtype,C,
                                     model_order, type= "augment",aug_type="not_foldover",
                                     D_opt_design,D_opt_x,x_unifopt,unif_crit,
                                     n_starts=500,D_opt_type="Bayesian",
                                     w_candidate=seq(0.05,0.95,0.05),use_w=FALSE,use_cpp=FALSE,
                                     K=0,me_index_daug, qe_index_daug,two_fi_index_daug,update_all_rows=FALSE,
                                     D_opt_tresh=0.8, U_opt_tresh=0.8, w_tresh=0.02, rng=c(),num_cores=NA){
  w_eff=list()
  #print(aug_type)
  #cont_index=which(C_dtype=="cont")
  #for(i in cont_index){
  #  C[[i]]=seq(-1,1,0.01)
  #}
  if(use_w){
    freeze_rows=nrow(x_matrix_start)
    print(paste("Freeze rows value is",freeze_rows))
    for(w in w_candidate){
      best_compound=0
      #cl=makeCluster(detectCores()-1, type="FORK")
      if(is.na(num_cores)){
        cl=parallel::makeCluster(as.numeric(Sys.getenv("PBS_RESC_TOTAL_PROCS")) - 1, type="FORK",outfile="")
      }else{
        cl=parallel::makeCluster(num_cores, type="FORK",outfile="")
      }
      doParallel::registerDoParallel(cl)
      #registerDoRNG(seed = 1234)
      #registerDoParallel(cores=detectCores()-1)
      #registerDoParallel(cores=as.numeric(Sys.getenv("PBS_RESC_TOTAL_PROCS")) - 1)
      mcoptions = list(preschedule=FALSE)
      chunks = foreach::getDoParWorkers()
      print(paste("Number of workers is:",chunks))
      multistart_results=foreach::foreach(i=1:n_starts, .options.multicore=mcoptions,.options.RNG=rng) %dorng% individual_starts_compound_unif(i=i,x_matrix_start=x_matrix_start,x_start=x_start,
                                                                                                                                       n=n,C_dtype=C_dtype,C=C,
                                                                                                                                       model_order=model_order, type=type, freeze_rows=freeze_rows,
                                                                                                                                       D_opt_design=D_opt_design,D_opt_x=D_opt_x,
                                                                                                                                       x_unifopt=x_unifopt,unif_crit=unif_crit,
                                                                                                                                       w=w,use_cpp=use_cpp,K=K,D_opt_type=D_opt_type,aug_type=aug_type,
                                                                                                                                       me_index_daug=me_index_daug, qe_index_daug=qe_index_daug,two_fi_index_daug=two_fi_index_daug,update_all_rows=update_all_rows)
      #stopImplicitCluster()
      parallel::stopCluster(cl)
      print(paste("Length of output is:",length(multistart_results)))
      for(i in 1:length(multistart_results)){
        current=multistart_results[[i]][[2]]
        x_current=multistart_results[[i]][[3]]
        des=multistart_results[[i]][[1]]
        if(current>best_compound){
          best_compound=current
          best_des_compound=des
          x_best=x_current
          deff_best=multistart_results[[i]][[4]]
          ueff_best=multistart_results[[i]][[5]]
        }
      }
      number_best_converegence=0
      criteria_list=c()
      for(i in 1:length(multistart_results)){
        criteria_list=c(criteria_list,multistart_results[[i]][[2]])
        if(isTRUE(all.equal(best_compound,multistart_results[[i]][[2]]))){
          number_best_converegence=number_best_converegence+1
        }
      }
      w_eff[[as.character(w)]]=list(x_best,best_des_compound,c("Weighted criteria"=best_compound,"D_eff"=deff_best,"U_eff"=ueff_best))
      print(paste("Number of iterations converged to the best value is",number_best_converegence))
      graphics::plot(x=1:length(multistart_results),y=criteria_list, main=paste("Plot of compound criteria",w), xlab="Iteration",ylab="Criteria_val")
      graphics::abline(h=best_compound, col="red")
      #print("Top 50 results")
      #print(criteria_list[order(criteria_list,decreasing=TRUE)[1:50]])
      #print("Compound Optimal designs")
      #print(w_eff)
    }
  }else{
    freeze_rows=nrow(x_matrix_start)
    continue=TRUE
    w_zero=FALSE
    w_start=0
    w_end=1
    w=(w_start+w_end)/2
    iter_w=0
    w_candidate=c(0.1,0.2,0.3,0.4)
    while(continue){
      iter_w=iter_w+1
      best_compound=0
      #cl=makeCluster(detectCores()-1, type="FORK")
      if(is.na(num_cores)){
        cl=parallel::makeCluster(as.numeric(Sys.getenv("PBS_RESC_TOTAL_PROCS")) - 1, type="FORK",outfile="")
      }else{
        cl=parallel::makeCluster(num_cores, type="FORK",outfile="")
      }
      doParallel::registerDoParallel(cl)
      #registerDoRNG(seed = 1234)
      #registerDoParallel(cores=detectCores()-1)
      #registerDoParallel(cores=as.numeric(Sys.getenv("PBS_RESC_TOTAL_PROCS")) - 1)
      mcoptions = list(preschedule=FALSE)
      chunks = foreach::getDoParWorkers()
      print(paste("Number of workers is:",chunks))
      multistart_results=foreach::foreach(i=1:n_starts, .options.multicore=mcoptions,.options.RNG=rng) %dorng% individual_starts_compound_unif(i=i,x_matrix_start=x_matrix_start,x_start=x_start,
                                                                                                                                       n=n,C_dtype=C_dtype,C=C,
                                                                                                                                       model_order=model_order, type=type, freeze_rows=freeze_rows,
                                                                                                                                       D_opt_design=D_opt_design,D_opt_x=D_opt_x,
                                                                                                                                       x_unifopt=x_unifopt,unif_crit=unif_crit,
                                                                                                                                       w=w,use_cpp=use_cpp,K=K,D_opt_type=D_opt_type,aug_type=aug_type,
                                                                                                                                       me_index_daug=me_index_daug, qe_index_daug=qe_index_daug,two_fi_index_daug=two_fi_index_daug)
      #stopImplicitCluster()
      parallel::stopCluster(cl)
      print(paste("Length of output is:",length(multistart_results)))
      for(i in 1:length(multistart_results)){
        current=multistart_results[[i]][[2]]
        x_current=multistart_results[[i]][[3]]
        des=multistart_results[[i]][[1]]
        if(current>best_compound){
          best_compound=current
          best_des_compound=des
          x_best=x_current
          deff_best=multistart_results[[i]][[4]]
          ueff_best=multistart_results[[i]][[5]]
        }
      }
      number_best_converegence=0
      criteria_list=c()
      for(i in 1:length(multistart_results)){
        criteria_list=c(criteria_list,multistart_results[[i]][[2]])
        if(isTRUE(all.equal(best_compound,multistart_results[[i]][[2]]))){
          #if(all(x_best==multistart_results[[i]][[3]])){
          number_best_converegence=number_best_converegence+1
          #}
        }
      }
      w_eff[[as.character(w)]]=list(x_best,best_des_compound,c("Weighted criteria"=best_compound,"D_eff"=deff_best,"U_eff"=ueff_best))
      print(paste("Number of iterations converged to the best value is",number_best_converegence))
      graphics::plot(x=1:length(multistart_results),y=criteria_list, main=paste("Plot of compound criteria",w), xlab="Iteration",ylab="Criteria_val")
      graphics::abline(h=best_compound, col="red")
      #print("Top 50 results")
      #print(criteria_list[order(criteria_list,decreasing=TRUE)[1:50]])
      #update w_start and w_end
      print(paste("W start",w_start))
      print(paste("W end", w_end))
      print(paste(w))
      print(paste("U efficiency is,",ueff_best))
      print(paste("D efficiency is,",deff_best))
      if(deff_best>=0.8 & ueff_best>=0.8){
        w_end=w
        w_start=w_start
        w=(w_start+w_end)/2
      }else if(deff_best>=D_opt_tresh & ueff_best<U_opt_tresh){
        w_end=w
        w_start=w_start
        w=(w_start+w_end)/2
      }else if(deff_best<D_opt_tresh & ueff_best>=U_opt_tresh){
        w_start=w
        w_end=w_end
        w=(w_start+w_end)/2
      }else if(deff_best<D_opt_tresh & ueff_best<U_opt_tresh){
        w_start=w
        w_end=w_end
        w=(w_start+w_end)/2
      }
      if(round(w_end-w_start,2)<=w_tresh | length(w_candidate)==0 | iter_w>10 | w_zero){
        if(!w_zero){
          w=0
          w_zero=TRUE
        }else{
          continue=FALSE
        }
      }
    }
    if(unif_crit=="WD2"){
      uopt=DiceDesign::discrepancyCriteria(x_unifopt)$DisW2/DiceDesign::discrepancyCriteria(D_opt_x)$DisW2
    }else if(unif_crit=="MD2"){
      unifopt=discM2_update(X=x_unifopt,update=FALSE,rescale=FALSE)
      disc_matrix=unifopt[[1]]
      disc_vector=unifopt[[2]]
      unifopt=sqrt((19/12)**(ncol(x_unifopt))-(2/nrow(x_unifopt))*sum(disc_vector)+(1/nrow(x_unifopt)**2)*sum(disc_matrix))
      dopt_unif=discM2_update(X=D_opt_x,update=FALSE,rescale=FALSE)
      disc_matrix=dopt_unif[[1]]
      disc_vector=dopt_unif[[2]]
      dopt_unif=sqrt((19/12)**(ncol(D_opt_x))-(2/nrow(D_opt_x))*sum(disc_vector)+(1/nrow(D_opt_x)**2)*sum(disc_matrix))
      uopt=unifopt/dopt_unif
    }
    w_eff[["1"]]=list(D_opt_design,c("Weighted criteria"=1,"D_eff"=1,"U_eff"=uopt))
    #print(w_eff)
  }
  return(w_eff)
}

