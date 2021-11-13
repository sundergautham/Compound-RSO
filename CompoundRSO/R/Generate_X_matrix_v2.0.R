#' function to sample from the candidate space
#' @description function to sample from the candidate space
#' @param C List of candidate space for each experimental factor
#' @return A row of the factor levels matrix sampled from the candidate space
#' @export
pick_rand= function(C){
  x=c()
  for(i in 1:length(C)){
    x=c(x,sample(C[[i]],1))
  }
  return(x)
}

#' function to create a starting design of run size n and specified model order
#' @description function to create a starting design of run size n and specified model order
#' @param n run size of the design
#' @param C List of candidate space for each experimental factor
#' @param model_order a string indicating the order of the regression model. "linear","quadratic" (Main Effects+Quadratic Effects), "2FI" (Main Effects + TwoFactor Interactions),"full", "non_standard"- for a user defined model order order
#' @param C_dtype a string vector indicating the data type of the experimental factors in the design either "cont"(continuous) or "cat"(categorical)
#' @param type a string indicating the type of design "non_augment" or "augment" (augment an existing deisgn)
#' @param me_index vector of index of Main Effects to be included in the design for model_order "non_standard"
#' @param qe_index vector index of Quadratic Effects to be included in the design for model_order "non_standard"
#' @param two_fi_index list of index of Two Factor Interactions to be included in the design for model_order "non_standard"
#' @param cubic_index vector of index of Cubic Effects to be included in the design for model_order "non_standard" (although a cubic design can be sampled, current design algorithms don't support cubic models)
#' @param quatric_index vector of index of Quatric Effects to be included in the design for model_order "non_standard" (although a quatric design can be sampled, current design algorithms don't support cubic models)
#' @param design a string indicating the type of design constructed. "Non-Bayesian" or "Bayesian" optimal design. The major difference is in how signularity of the matrix is evaluated. In a Bayesian optimal design X'X+K(prior variance matrix) is checked for singularity
#' @return a list of factor levels matrix and design matrix of the specified model order
#'  \itemize{
#'        \item x- a matrix of factor levels in the design
#'        \item x_matrix- a design matrix of the specified order including the intercept
#'        }
#'  @examples
#'   gen_x(n,C,model_order="full",C_dtype,type="non_augment",me_index=c(),qe_index=c(),two_fi_index=c(),cubic_index=c(),quatric_index=c(),design="Non-Bayesian")
#' @export
gen_x=function(n,C,model_order,C_dtype,type="non_augment",
               me_index=c(),qe_index=c(),two_fi_index=c(),
               cubic_index=c(),quatric_index=c(),
               design="Non-Bayesian"){
  x=matrix(ncol=length(C))
  if(model_order=="linear"){
    p=1+length(C)
  }else if(model_order=="quad"){
    p=1+length(C)
    q=length(which(C_dtype=="cont"))
    p=p+q
  }else if(model_order=="2FI"){
    p=length(C)
    p=1+p+p*(p-1)/2
  }else if(model_order=="full"){
    #model is full
    p=length(C)
    p=1+p+p*(p-1)/2
    q=length(which(C_dtype=="cont"))
    p=p+q
  }else if(model_order=="non_standard"){
    p=1+length(me_index)+length(qe_index)+length(two_fi_index)+
      length(cubic_index)+length(quatric_index)
  }
  x_matrix=matrix(ncol=p)
  for(i in 1:n){
    if(model_order=="non_standard"){
      x=rbind(pick_rand(C)[1:length(me_index)],x)
    }else{
      x=rbind(pick_rand(C),x)
    }
    row=f_x(x=x[1,],model_order=model_order,C_dtype=C_dtype,
            me_index=me_index,qe_index = qe_index,
            two_fi_index = two_fi_index,
            cubic_index=cubic_index,
            quatric_index=quatric_index)
    x_matrix=rbind(row,x_matrix)
  }
  x=x[-nrow(x),]
  x_matrix=x_matrix[-nrow(x_matrix),]
  r_cond=rcond(t(x_matrix)%*%x_matrix)
  #condition number of the X'X matrix is used to generate full column rank starting designs
  if(r_cond>1e-15 | type=="augment" | design=="Bayesian"){
    return(list("x"=x,"x_matrix"=x_matrix))
  }else{
    print("Signular Design")
    gen_x(n=n,C=C,model_order=model_order,C_dtype=C_dtype,type=type,
          me_index=me_index,qe_index=qe_index,two_fi_index=two_fi_index,
          cubic_index=cubic_index,quatric_index=quatric_index,
          design=design)
  }
}

#' function to create a row of the design matrix of a specified model order
#' @description function to create a row of the design matrix of a specified model order
#' @param x a vector or row matrix of factor levels to generate the design matrix for
#' @param model_order a string indicating the order of the regression model. "linear","quadratic" (Main Effects+Quadratic Effects), "2FI" (Main Effects + TwoFactor Interactions),"full", "non_standard"- for a user defined model order order
#' @param C_dtype a string vector indicating the data type of the experimental factors in the design either "cont"(continuous) or "cat"(categorical)
#' @param me_index vector of index of Main Effects to be included in the design for model_order "non_standard"
#' @param qe_index vector index of Quadratic Effects to be included in the design for model_order "non_standard"
#' @param two_fi_index list of index of Two Factor Interactions to be included in the design for model_order "non_standard"
#' @param cubic_index vector of index of Cubic Effects to be included in the design for model_order "non_standard" (although a cubic design can be sampled, current design algorithms don't support cubic models)
#' @param quatric_index vector of index of Quatric Effects to be included in the design for model_order "non_standard" (although a quatric design can be sampled, current design algorithms don't support cubic models)
#' @return a row of the design matrix for the specified factor levels and  model_order
#'  @examples
#'   f_x(x,model_order="full",C_dtype, me_index=c(),qe_index=c(),two_fi_index=c(),cubic_index=c(),quatric_index=c())
#' @export
f_x= function(x,model_order,C_dtype, me_index=c(),qe_index=c(),
              two_fi_index=c(),cubic_index=c(),quatric_index=c()){
  x_copy=x
  if(model_order=="linear"){
    fx_prime=c(1)
    fx_prime=c(fx_prime,x)
  }
  else if(model_order=="2FI"){
    fx_prime=x
    l=length(fx_prime)
    for(i in 1:l-1){
      for(j in (i+1):l){
        fx_prime=c(fx_prime,x[i]*x[j])
      }
    }
    fx_prime=c(1,fx_prime)
  }
  else if(model_order=='full'){
    fx_prime=x
    l=length(fx_prime)
    for(i in 1:l-1){
      for(j in (i+1):l){
        fx_prime=c(fx_prime,x[i]*x[j])
      }
    }
    ind=which(C_dtype=="cont")
    quad=sapply(ind,FUN=function(j){x_copy[j]**2})
    fx_prime=c(1,fx_prime,quad)
  }
  else if(model_order=='quad'){
    fx_prime=x
    ind=which(C_dtype=="cont")
    quad=sapply(ind,FUN=function(j){x_copy[j]**2})
    fx_prime=c(1,fx_prime,quad)
  }
  else if(model_order=="non_standard"){
    fx_prime=c(1)
    for(i in me_index){
      fx_prime=c(fx_prime,x[i])
    }
    for(i in two_fi_index){
      fx_prime=c(fx_prime,x[i[1]]*x[i[2]])
    }
    for(i in qe_index){
      fx_prime=c(fx_prime,x[i]**2)
    }
    for(i in cubic_index){
      fx_prime=c(fx_prime,x[i]**3)
    }
    for(i in quatric_index){
      fx_prime=c(fx_prime,x[i]**4)
    }
  }
  fx_prime=matrix(fx_prime,nrow=1)
  return(fx_prime)
}

#' an internal function that computes the intermediate steps in D-optimal designs
#' @description an internal function that computes the intermediate steps in computing the change in determinant value for each coordinate exchange in the construction of D-optimal designs
#' @param x_i a row of the design matrix
#' @param x_j a row of the design matrix
#' @param D current X'X inverse matrix
#' @return a matrix
#' @export
v_x=function(x_i,x_j,D){
  v=matrix(x_i,nrow=1)%*%D%*%matrix(x_j,ncol=1)
  #print(paste("Value of v",v))
  #return(round(v,8))
  return(v)
}

#' an internal function that computes the intermediate steps in I-optimal design
#' @description an internal function that computes the intermediate steps in computing the change in determinant value for each coordinate exchange in the construction of I-optimal designs
#' @param x_i a row of the design matrix
#' @param x_j a row of the design matrix
#' @param D current X'X inverse matrix
#' @param W Moment matrix
#' @return a matrix
#' @export
phi_x=function(x_i,x_j,D,W){
  phi=matrix(x_i,nrow=1)%*%D%*%W%*%D%*%matrix(x_j,ncol=1)
  return(round(phi,8))
}


#' an internal function to compute the change in determinant value for a coordinate exhange in a D-optimal design
#' @description an internal function to compute the change in determinant value for a coordinate exhange in a D-optimal design
#' @param x the row vector with the new coordinate value
#' @param x_i the current row vector that needs to be exchanged
#' @param D current X'X inverse values
#' @export
delta_D=function(x,x_i,D){
  vx=v_x(x,x,D)
  v_x_xi=v_x(x,x_i,D)
  vi=v_x(x_i,x_i,D)
  del=1+(vx-vi)+(v_x_xi^2-vx*vi)
  return(del)
}
