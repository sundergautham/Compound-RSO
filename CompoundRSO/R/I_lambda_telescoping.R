#' function to compute the optimal x* for a given regression function
#' @param i is the bootstrap sample index
#' @param b vector that is a sample from the multivariate normal distribution of the regression parameters
#' @param fx function to compute regression prediction t(x)%*%b
#' @param model_order a string indicating the order of the regression model. "linear","quadratic" (Main Effects+Quadratic Effects), "2FI" (Main Effects + TwoFactor Interactions),"full", "non_standard"- for a user defined model order order
#' @param C_dtype a string vector indicating the data type of the experimental factors in the design either "cont"(continuous) or "cat"(categorical)
#' @param me_index vector of index of Main Effects to be included in the design for model_order "non_standard"
#' @param qe_index vector index of Quadratic Effects to be included in the design for model_order "non_standard"
#' @param two_fi_index list of index of Two Factor Interactions to be included in the design for model_order "non_standard"
#' @param lower a vector of lower limits of the experimental region for all factors (lower level constraint to the optimization)
#' @param upper a vector of upper limits of the experimental region for all factors (upper level constraint to the optimization)
#' @param nfactor number of experimental factors
#' @param C list of candidate space for each experimental factor in the design
#' @param minimize a boolean value indicating whether to minimize or maximize. Default value is minimize
#' @return returns a vector of optimum x values for the constrained optimization within the experimental region
#' @description function to compute the optimal x* for a given regression function
#' @export
individual_bootstrap=function(i,b,fx,model_order="non_standard",C_dtype,me_index,qe_index,two_fi_index,lower,upper,nfactor,C,minimize=TRUE){
  print(paste("Random start",i))
  parmat=sapply(C,FUN=sample,size=20,replace=TRUE)
  if(minimize){
    opt_multiplier=1
  }else{
    opt_multiplier=-1
  }
  opt_vals=optimr::multistart(parmat = parmat, fn=fx,method = "L-BFGS-B",
                      b=opt_multiplier*matrix(b),model_order=model_order,C_dtype=C_dtype,
                      me_index=me_index,
                      qe_index=qe_index,
                      two_fi_index=two_fi_index,
                      lower = lower,upper=upper)
  opt_vals=opt_vals[opt_vals[,"convergence"]==0,]
  x_star=as.matrix(opt_vals[which(opt_vals$value==min(opt_vals$value)),1:nfactor])
  print(paste("x_star",x_star))
  return(x_star)
}


#' function to compute the bootstrap confidence region for the location of the optimal response
#' @param beta current estimate of the regression parameters.
#' @param Sigma the variance of the regression parameters estimated.
#' @param B Number of bootstrap samples
#' @param model_order a string indicating the order of the regression model. "linear","quadratic" (Main Effects+Quadratic Effects), "2FI" (Main Effects + TwoFactor Interactions),"full", "non_standard"- for a user defined model order order
#' @param C_dtype a string vector indicating the data type of the experimental factors in the design either "cont"(continuous) or "cat"(categorical)
#' @param me_index vector of index of Main Effects to be included in the design for model_order "non_standard"
#' @param qe_index vector index of Quadratic Effects to be included in the design for model_order "non_standard"
#' @param two_fi_index list of index of Two Factor Interactions to be included in the design for model_order "non_standard"
#' @param lower a vector of lower limits of the experimental region for all factors (lower level constraint to the optimization)
#' @param upper a vector of upper limits of the experimental region for all factors (upper level constraint to the optimization)
#' @param nfactor number of experimental factors
#' @param C list of candidate space for each experimental factor in the design
#' @param fx function to compute regression prediction t(x)%*%b
#' @param minimize a boolean value indicating whether to minimize or maximize. Default value is minimize
#' @param rng seed for random number generation
#' @param num_cores number of cores for parallel processing. Default value is NA, in which case Sys.getenv("PBS_RESC_TOTAL_PROCS") is used to get the number of cores in a HPC. Otherwise, the user needs to provide a value. Use 1 for serial execution.
#' @description function to compute the bootstrap confidence region for the optimal x values. The function uses bootstrap samples to sample from the multivariate normal distribution of the regression parameters.The confidence region is the subspace of the experimental region that contains the optimal x values for each sample in the bootstrap. The confidence region can be the entire experimental region when there is high uncertainty in the regression estimate.The code uses parallel computing to compute the confidence region and currently works only on linux machines.
#' @return returns a list of
#'    \itemize{
#'        \item confidence_set_x_star- the x_star values for each bootstrap sample of regression parameters
#'        \item confidence_set_beta- the bootstrap sample of regression parameters
#'        }
#' @examples
#' \dontrun{
#'   bootstrap_confidence_set(beta,Sigma,B=1000,me_index,qe_index,two_fi_index,C_dtype,nfactor,C,model_order="non_standard", lower, upper, fx, rng)
#' }
#' @export
#' @import doRNG
bootstrap_confidence_set=function(beta,Sigma,B=1000,me_index,qe_index,two_fi_index,C_dtype,nfactor,C,model_order="non_standard", lower, upper, fx,minimize=TRUE, rng=c(), num_cores=NA){
  i=NULL
  B_samples=mvtnorm::rmvnorm(n = B, mean=beta, sigma=Sigma, method="svd")
  print(dim(B_samples))
  chi_stat=stats::qchisq(0.05,length(beta),lower.tail = FALSE)
  Sigma_inv=solve(Sigma)
  #removes the extreme values in the bootstrap sample
  conf_95=apply(B_samples,MARGIN=1,FUN=function(x){matrix((x-beta),nrow=1)%*%Sigma_inv%*%matrix((x-beta)) <=chi_stat})
  confidence_set_beta=B_samples[conf_95,]
  print("Dimension of Bootstrap samples")
  print(dim(confidence_set_beta))
  print("Dimension of unique Bootstrap samples")
  print(dim(unique(confidence_set_beta)))
  parmat=sapply(C,FUN=sample,size=20,replace=TRUE)
  #print(dim(parmat))
  #print(me_index)
  #print(qe_index)
  #print(two_fi_index)
  #print(C_dtype)
  #print(nfactor)
  #print(beta)
  #print(paste("length of model terms",length(me_index)+length(qe_index)+length(two_fi_index)+1))
  if(minimize){
    opt_multiplier=1
  }else{
    opt_multiplier=-1
  }
  opt_vals=optimr::multistart(parmat = parmat, fn=fx,method = "L-BFGS-B",
                      b=opt_multiplier*matrix(beta),model_order="non_standard",C_dtype=C_dtype,
                      me_index=me_index,
                      qe_index=qe_index,
                      two_fi_index=two_fi_index,
                      lower = rep(-1,nfactor),upper=rep(1,nfactor))
  opt_vals=opt_vals[opt_vals[,"convergence"]==0,]
  x_star=as.matrix(opt_vals[which.min(opt_vals$value),1:nfactor])
  confidence_set_x_star=x_star
  #for(i in 1:nrow(confidence_set_beta)){
  #  parmat=sapply(C,FUN=sample,size=20,replace=TRUE)
  #  opt_vals=multistart(parmat = parmat, fn=fx,method = "L-BFGS-B",
  #                      b=matrix(confidence_set_beta[i,]),model_order="non_standard",C_dtype=C_dtype,
  #                      me_index=me_index,
  #                      qe_index=qe_index,
  #                      two_fi_index=two_fi_index,
  #                      lower = rep(-1,nfactor),upper=rep(1,nfactor))
  #  opt_vals=opt_vals[opt_vals[,"convergence"]==0,]
  #  x_star=as.matrix(opt_vals[which(opt_vals$value==min(opt_vals$value)),1:nfactor])
  #  confidence_set_x_star=rbind(confidence_set_x_star,x_star)
  #}
  if(is.na(num_cores)){
    cl=parallel::makeCluster(as.numeric(Sys.getenv("PBS_RESC_TOTAL_PROCS")) - 1, type="FORK",outfile="")
  }else{
    cl=parallel::makeCluster(num_cores, type="FORK",outfile="")
  }
  doParallel::registerDoParallel(cl)
  #registerDoParallel(cores=detectCores()-1)
  mcoptions = list(preschedule=TRUE)
  chunks = foreach::getDoParWorkers()
  print(paste("Number of workers is:",chunks))
  #print(nfactor)
  confidence_set_x_star=foreach::foreach(i=1:nrow(confidence_set_beta), .options.multicore=mcoptions,.combine=rbind,.options.RNG=rng) %dorng%  individual_bootstrap(i=i,b=matrix(confidence_set_beta[i,]),fx=fx,model_order=model_order,
                                                                                                                                                            C_dtype=C_dtype,me_index=me_index,qe_index=qe_index,
                                                                                                                                                            two_fi_index=two_fi_index,lower=lower,
                                                                                                                                                            upper=upper,nfactor=nfactor,C=C,minimize=minimize)
  #stopImplicitCluster()
  parallel::stopCluster(cl)
  #print(confidence_set_x_star)
  confidence_set_x_star=rbind(x_star,confidence_set_x_star)
  print(paste("Dim of output is:",dim(confidence_set_x_star)))
  return(list("confidence_set_x_star"=confidence_set_x_star,"confidence_set_beta"=confidence_set_beta))
}


#' function to compute the x'b value for any given x and b
#' @param x vector of factor levels in a row of the design
#' @param b vector of regression parameters
#' @param model_order a string indicating the order of the regression model. "linear","quadratic" (Main Effects+Quadratic Effects), "2FI" (Main Effects + TwoFactor Interactions),"full", "non_standard"- for a user defined model order order
#' @param C_dtype a string vector indicating the data type of the experimental factors in the design either "cont"(continuous) or "cat"(categorical)
#' @param me_index vector of index of Main Effects to be included in the design for model_order "non_standard"
#' @param qe_index vector index of Quadratic Effects to be included in the design for model_order "non_standard"
#' @param two_fi_index list of index of Two Factor Interactions to be included in the design for model_order "non_standard"
#' @param cubic_index vector of index of Cubic Effects to be included in the design for model_order "non_standard"
#' @param quatric_index vector of index of Quatric Effects to be included in the design for model_order "non_standard"
#' @description function to compute the x'b value for any given x and b
#' @return returns a numeric value numeric value of t(x)%*%b
#' @import doRNG
#' @export
fx=function(x,b, model_order,C_dtype,me_index=c(),qe_index,two_fi_index=c(),
            cubic_index=c(),quatric_index=c()){
  x=f_x(x=x,model_order = model_order,C_dtype = C_dtype,
        me_index = me_index,qe_index = qe_index,two_fi_index = two_fi_index,
        cubic_index = cubic_index,quatric_index = quatric_index)
  return(x%*%b)
}

#' internal function to compute the change in the objective function of I-optimal design for each coordinate exchange
#' @param x the row vector with the new coordinate value
#' @param x_i the current row vector that needs to be exchanged
#' @param D Current (X'X)inv matrix
#' @param W Moment matrix computed for a given lambda
#' @description internal function to compute the change in the objective function of I-optimal design for each coordinate exchange
#' @return returns the the change in I optimal objective function for the coordinate exchange
#' @export
# @examples
#   delta_I(x,x_i,D,W)

delta_I=function(x,x_i,D,W){
  phix=phi_x(x,x,D,W)
  phixi=phi_x(x_i,x_i,D,W)
  phi_x_xi=phi_x(x,x_i,D,W)
  vx=v_x(x,x,D)
  v_x_xi=v_x(x,x_i,D)
  vi=v_x(x_i,x_i,D)
  del_D=1+(vx-vi)+(v_x_xi^2-vx*vi)
  del_I=(1-vi)*phix+2*(v_x_xi*phi_x_xi)-((1-vx)*phixi)
  del_I=del_I/round(del_D,8)
  return(del_I)
}


#' function for constructing an I-optimal design using coordinate exchange
#' @param x settings matrix of the design
#' @param x_matrix design matrix of the specified order including the intercept
#' @param C list of candidate space for each experimental factor in the design
#' @param model_order a string indicating the order of the regression model. "linear","quadratic" (Main Effects+Quadratic Effects), "2FI" (Main Effects + TwoFactor Interactions),"full", "non_standard"- for a user defined model order order
#' @param C_dtype a string vector indicating the data type of the experimental factors in the design either "cont"(continuous) or "cat"(categorical)
#' @param type a string indicating the type of design augmentation. "non_augment" or "augment" (augment an existing design)
#' @param W Moment matrix computed for a given lambda
#' @param freeze_rows the rows of the matrix to freeze while augmenting, 0 when type is "non_augment"
#' @param me_index vector of index of Main Effects to be included in the design for model_order "non_standard"
#' @param qe_index vector index of Quadratic Effects to be included in the design for model_order "non_standard"
#' @param two_fi_index list of index of Two Factor Interactions to be included in the design for model_order "non_standard"
#' @param telescoping a boolean value taking TRUE or FALSE. Telescoping reduces the candidate factors the local region of the optimum design
#' @description function for constructing an I-optimal design using coordinate exchange
#' @return a list of
#'  \itemize{
#'        \item Iopt_design- I-optimal design matrix of the specified order including the intercept from the random starts
#'        \item Iopt_obj- The value of the I-optimal objective value
#'        \item Iopt_x- a matrix with factor levels of the I-optimal design
#'        }
#'  @examples
#'   coor_exch_I_augment(x,x_matrix,C,model_order,C_dtype, W, freeze_rows=0,type="non_augment",me_index=c(),qe_index=c(),two_fi_index=c(),telescoping=FALSE)
#' @export
coor_exch_I_augment_R= function(x,x_matrix,C,model_order,C_dtype, W, freeze_rows=0,
                              type="non_augment",
                              me_index=c(),qe_index=c(),two_fi_index=c(),telescoping=FALSE){
  x_start=x
  #computes the inverse of X'X
  D= solve(t(x_matrix)%*%x_matrix)
  n=dim(x_matrix)[1]
  #saves a copy of D
  D_local=D
  #computes phi(x) for each row
  start=freeze_rows+1
  phi_i= sapply(seq_along(start:n),
                FUN=function(j){phi_x(x_i=x_matrix[j,],
                                      x_j=x_matrix[j,],D=D_local,W=W)})
  #computes v(x) for each row
  v_i= sapply(seq_along(start:n),
              FUN=function(j){v_x(x_i=x_matrix[j,],x_j=x_matrix[j,],D=D_local)})

  #computes the ratio phi/1-v
  phi_i=phi_i/(1-v_i)

  #orders the rows of the x matrix in ascending order for exchange, the rows with smallest values are deleted first.
  ord=order(x=phi_i, decreasing = FALSE)

  ord=freeze_rows+ord

  #print(phi_i)
  #print("order of excecution")
  #print(ord)

  #saves a local copy of x_matrix
  x_matrix_local=x_matrix

  #epsilon for monitoring convergence
  epsilon=1

  #Current L-optimality score
  current=sum(diag(D_local%*%W))
  print(paste("current objective:",current))

  while(epsilon>1e-15){
    # at the end of one complete exhchange of x_matrix_local, the updated  matrix is stored
    x_matrix=x_matrix_local
    prev=current
    x_best=x
    #for(i in ord){
    for(i in (freeze_rows+1):dim(x_matrix_local)[1]){
      #print(paste("exchange of row number:", i))
      for(j in 1:length(C)){
        #print(paste("exchange of col number:", j))
        det_val=c()
        #current value of trace
        tr_current=sum(diag(D_local%*%W))
        if(telescoping){
          point=x_start[i,j]
          candidate_points=seq(max(-1,point-0.09),min(1,point+0.09),0.01)
        }else{
          candidate_points=C[[j]]
        }
        for(k in candidate_points){
          x_local=x[i,]
          #point exchange
          x_local[j]=k
          x_expl=f_x(x_local,model_order,C_dtype,
                     me_index=me_index,qe_index=qe_index,
                     two_fi_index=two_fi_index)
          #x_new=x_matrix_local
          #x_new[i,]=x_expl
          #if(determinant(t(x_new)%*%x_new)$modulus[1]!=-Inf){
          #d_val=delta_I(x=x_expl,x_i = x_matrix_local[i,],
          #        D=D_local,W=W)
          #del_D=d_val[[2]]
          #print(paste("del_D",del_D))
          #d_val=d_val[[1]]
          #if(del_D!=0){
          #det_val=c(det_val,d_val)
          #}else{det_val=c(det_val,NaN)}
          #}else{
          #  det_val=c(det_val,NaN)
          #}
          #x_new=x_matrix_local
          #x_new[i,]=x_expl
          #compute the D_inverse after exchange
          D_new=rank2_inverse(x_expl,x_matrix_local[i,],D_local)
          #check if singular
          if(D_new[1]!="Singular"){
            det_val=c(det_val,tr_current-sum(diag(D_new%*%W)))
          }else{
            det_val=c(det_val,NaN)
          }
          #if(is.finite(d_val)){
          #  det_val=c(det_val,d_val)
          #}else{
          #  det_val=c(det_val,NaN)
          #}
        }
        det_val=round(det_val,8)
        #print("change in objective function")
        #print(det_val)
        #det_val=ifelse(abs(det_val)==Inf,NaN,det_val)
        #print(det_val)
        #if(!any(is.na(det_val)) & !any(abs(det_val)==Inf)){
        #if all of de_val is not NA
        if(!all(is.na(det_val))){
          ind=which(det_val==max(det_val,na.rm=TRUE))[1]
          #print(paste("max=",ind))
          #make the swap only if the delta is greater than 0
          if(det_val[ind]>0){
            x[i,j]=candidate_points[ind]
            D_local=rank2_inverse(f_x(x[i,],model_order,C_dtype,
                                      me_index=me_index,qe_index=qe_index,
                                      two_fi_index=two_fi_index),
                                  x_matrix_local[i,],D_local)
            x_matrix_local[i,]=f_x(x[i,],model_order,C_dtype,
                                   me_index=me_index,qe_index=qe_index,
                                   two_fi_index=two_fi_index)
          }
          #print("x matrix after interchange")
          #print(x_matrix_local)
          #print(paste("det of inf matrix after exchange:",determinant(t(x_matrix_local)%*%x_matrix_local)$modulus[1]))
          #print("___________________________________")
        }else{
          #return(list(x_matrix,100000,x_best))
          print(x_matrix_local)
          print(j)
          print(i)
          stop("Singular Design Obtained")
        }
        # if variance after update is
        #trych=tryCatch({
        #  D_local=solve(t(x_matrix_local)%*%x_matrix_local)
        #},error=function(err){
        #  print(x_matrix_local)
        #  print(which(det_val==min(det_val)))
        #  print("det value is:")
        #  print(det_val)
        #  print(ind)
        #  print(err)
        #  #return(list(x_matrix,100000,x_best))
        #break
        #  NULL
        #})
        #print(paste("trych=",trych))
        #if(is.null(trych)){
        #  return(list(x_matrix_local,10000,x))
        #}
      }
      #print(paste("det of inf matrix:",determinant(t(x_matrix_local)%*%x_matrix_local)$modulus[1]))
      #if(determinant(t(x_matrix_local)%*%x_matrix_local)$modulus[1]==-Inf){
      #  #return(list(x_matrix,100000,x_best))
      #  print(x_matrix_local)
      #}
      current=sum(diag(D_local%*%W))
    }
    print(paste("current objective:",current))
    #calculates epsilon
    epsilon=prev-current
    print(paste("prev=",prev))
    print(epsilon)
    #epsilon=abs((current - prev) / prev)

    phi_i= sapply(seq_along(start:n),
                  FUN=function(j){phi_x(x_i=x_matrix[j,],x_j=x_matrix[j,],D=D_local,W=W)})

    v_i= sapply(seq_along(start:n),
                FUN=function(j){v_x(x_i=x_matrix[j,],x_j=x_matrix[j,],D=D_local)})

    phi_i=phi_i/(1-v_i)

    ord=order(x=phi_i, decreasing = FALSE)
    ord=freeze_rows+ord
    #print("next iteration")
    #print(phi_i)
    #print("order of excecution")
    #print(ord)
  }
  return(list("Iopt_design"=x_matrix,"Iopt_obj"=prev,"Iopt_x"=x_best))
}

#' function for initiating the coordinate exchange for I optimal design
#' @param x_matrix_start the design matrix of the specified order including the intercept of the random starting
#' @param x_start the settings matrix of the random starting design
#' @param W_mat Moment matrix computed for a given lambda
#' @param n total number of runs in the design
#' @param C list of candidate space for each experimental factor in the design
#' @param model_order a string indicating the order of the regression model. "linear","quadratic" (Main Effects+Quadratic Effects), "2FI" (Main Effects + TwoFactor Interactions),"full", "non_standard"- for a user defined model order order
#' @param C_dtype a string vector indicating the data type of the experimental factors in the design either "cont"(continuous) or "cat"(categorical)
#' @param type a string indicating the type of design augmentation. "non_augment" or "augment" (augment an existing design)
#' @param me_index vector of index of Main Effects to be included in the design for model_order "non_standard"
#' @param qe_index vector index of Quadratic Effects to be included in the design for model_order "non_standard"
#' @param two_fi_index list of index of Two Factor Interactions to be included in the design for model_order "non_standard"
#' @param n_starts number of random starts, default at 100
#' @description function for initiating the coordinate exchange for I optimal design for a given number of random starts and returns the best design. This is a serial version of the code. By default the parallel version of the code will be run.
#' @return a list of
#'  \itemize{
#'        \item I_opt_design- the best I-optimal design matrix of the specified order including the intercept from the random starts
#'        \item Iopt_x- a matrix with factor levels of the best I-optimal design from the random starts
#'        \item Iopt_obj- The value of best objective value of the I-optimal designs from the random starts
#'        }
#'  @examples
#'   construct_des_I_lambda_serial(x_matrix_start,x_start,W_mat,n,C_dtype, C, model_order, type= "augment", me_index = c(),qe_index = c(),two_fi_index = c(), n_starts=100)
#' @export
construct_des_I_lambda_serial=function(x_matrix_start,x_start,W_mat,n,C_dtype, C,
                                model_order, type= "augment",
                                me_index = c(),qe_index = c(),
                                two_fi_index = c(), n_starts=100){
  #best/sum(diag(x_matrix_start%*%W_mat))
  #x_start=x_dsd
  if(type != "augment"){
    best_augment=100000
    for(i in 1:n_starts){
      print(paste("random start",i))
      x=gen_x(n=n,C=C,model_order=model_order,C_dtype=C_dtype,
              type="non_augment",
              me_index = me_index,qe_index = qe_index,
              two_fi_index = two_fi_index)
      #print(x)
      x_matrix=x[[2]]
      x=x[[1]]
      if(class(x)!="matrix"){x=matrix(x)}

      print(paste("Coordinate exchange for try",i))
      des=coor_exch_I_augment(x=x,x_matrix=x_matrix,C=C,model_order =model_order,
                              C_dtype=C_dtype,W=W_mat,freeze_rows = 0,
                              me_index=me_index,qe_index=qe_index,
                              two_fi_index=two_fi_index)
      current=des[[2]]
      x_current=des[[3]]
      des=des[[1]]
      if(current<best_augment){
        best_augment=current
        best_des_augment=des
        x_best_augment=x_current
        #break
      }
      print("------------------------------------")
    }
  }else{
    freeze_rows=nrow(x_matrix_start)
    best_augment=100000
    for(i in 1:n_starts){
      n=n
      print(paste("random start",i))
      singular=TRUE
      while(singular){
        x=gen_x(n=n-freeze_rows,C,model_order=model_order,C_dtype,
                type="augment",
                me_index = me_index,qe_index = qe_index,
                two_fi_index = two_fi_index)
        #print(x)
        x_matrix=x[[2]]
        x_matrix=rbind(x_matrix_start,x_matrix)
        x=x[[1]]
        if(class(x)!="matrix"){x=matrix(x)}
        #print(x)
        #print(x_matrix)
        x=rbind(x_start,x)
        if(rcond(t(x_matrix)%*%x_matrix)>1e-5){
          singular=FALSE
        }
      }
      #print(x)
      #x_matrix=cnn_doe_mat
      #x=x_cnn_doe
      if(class(x)!="matrix"){
        x=matrix(x)
      }
      print(paste("Coordinate exchange for try",i))
      des=coor_exch_I_augment(x=x,x_matrix=x_matrix,C=C,model_order =model_order,
                              C_dtype=C_dtype,W=W_mat,freeze_rows = freeze_rows,
                              me_index=me_index,qe_index=qe_index,
                              two_fi_index=two_fi_index)
      current=des[[2]]
      x_current=des[[3]]
      des=des[[1]]
      if(current<best_augment ){
        best_augment=current
        best_des_augment=des
        x_best_augment=x_current
        #break
      }
      print("------------------------------------")
    }
  }
  return(list("I_opt_design"=best_des_augment,"Iopt_x"=x_best_augment,"Iopt_obj"=best_augment))
}


#' function for computing the moment matrix for a given lambda using Monte Carlo sampling
#' @param low a vector of lower bound on the experimental region to compute the moment matrix
#' @param high a vector of upper bound on the experimental region to compute the moment matrix
#' @param shape1 shape of the beta distribution (set at a default of 1)
#' @param shape2 shape of the beta distribution (set at a default of 1)
#' @param num the number of samples in the Monte Carlo estimate
#' @param C_dtype a string vector indicating the data type of the experimental factors in the design either "cont"(continuous) or "cat"(categorical)
#' @param model_order a string indicating the order of the regression model. "linear","quadratic" (Main Effects+Quadratic Effects), "2FI" (Main Effects + TwoFactor Interactions),"full", "non_standard"- for a user defined model order order
#' @param me_index vector of index of Main Effects to be included in the design for model_order "non_standard"
#' @param qe_index vector index of Quadratic Effects to be included in the design for model_order "non_standard"
#' @param two_fi_index list of index of Two Factor Interactions to be included in the design for model_order "non_standard"
#' @description function for computing the moment matrix for a given lambda using Monte Carlo sampling. The specific Lambda used by default is uniform (tophat prior) weights in a subregion of the experimental region and 0 elsewhere.
#' @return returns a moment matrix matrix W_mat
#'  @examples
#'   E_W(low, high, shape1,shape2,num=100000,C_dtype,model_order,me_index=c(),qe_index=c(),two_fi_index=c())
#' @export
E_W= function(low, high, shape1=1,shape2=1,
              num=100000,C_dtype,model_order,
              me_index=c(),qe_index=c(),two_fi_index=c()){
  sum=0
  inc=0.001
  for(i in 1:num){
    x=c()#rbeta(2,shape1=shape1[1],shape2=shape2[1])*(high-low)+low
    #print(x)
    for(j in 1:length(shape1)){
      x_j=stats::rbeta(1,shape1=shape1[j],shape2=shape2[j])*(high[j]-low[j])+low[j]
      x=c(x,x_j)
    }
    #print(x)
    x_mat=f_x(x,model_order = model_order,
              me_index = me_index,qe_index=qe_index,
              two_fi_index = two_fi_index,
              C_dtype = C_dtype)
    #print(dim(x_mat))
    m=t(x_mat)%*%x_mat
    sum=sum+m
  }
  return(sum/num)
}

#Description:
##function for constructing an compound D-I design using coordinate exchange. The D-I optimal design is the exploration-exploitation tradeoff design

#Usage:
##coor_exch_compound(x,x_matrix,C,model_order,me_index = c(),qe_index = c(),two_fi_index = c(),C_dtype, W, freeze_rows=0,x_matrix_Dopt,x_matrix_Iopt,w,K=0,D_opt_type="Bayesian",me_index_daug, qe_index_daug,two_fi_index_daug,telescoping=FALSE)

#Arguments:
## x: a matrix of factor levels in the design
## x_matrix: design matrix of the specified order including the intercept
## C: list of candidate points for all factors in the design
## model_order: string indicating linear,quadratic (Main Effects+Quadratic Effects), 2FI (Main Effects + TwoFactor Interactions),full, non_standard for a user defined order
## me_index: vector of index of Main Effects to be included in the I-optimal design for model_order "non_standard"
## qe_index: vector index of Quadratic Effects to be included in the I-optimal design for model_order "non_standard"
## two_fi_index: list of index of Two Factor Interactions to be included in the I-optimal design for model_order "non_standard"
## C_dtype: data type of the factors in the design either "cont"(continuous) or "cat"(categorical)
## W: Moment matrix computed for a given lambda
## freeze_rows: the rows of the matrix to freeze while augmenting, 0 when type is "non_augment"
## x_matrix_Dopt: the D-optimal design matrix
## x_matrix_Iopt: the I-optimal design matrix
## w: a numeric value of weight that is used to compute the compound optimal design
## K: a diagonal matrix specifying the prior variance in Bayesian D-optimal design set to 0.001 for potential terms and 0 for primary terms
## D_opt_type: Type of D-optimal design Bayesian or regular D-optimal design
## me_index_daug: vector of index of Main Effects to be included in the D-optimal design for model_order "non_standard"
## qe_index_daug: vector index of Quadratic Effects to be included in the D-optimal design for model_order "non_standard"
## two_fi_index_daug: list of index of Two Factor Interactions to be included in the D-optimal design for model_order "non_standard"
## type:"non_augment" or "augment" (augment an existing design)
## telescoping: a boolean value taking TRUE or FALSE. Telescoping reduces the candidate factors the local region of the optimum design

#Value: return a list of x_matrix,prev,x_best
## x_matrix: a I-optimal design matrix of the specified order including the intercept
## optimal_value: The value of log(det(X'X+K)) for the D-optimal design returned
## x: a matrix with factor levels of the Bayesian D-optimal design

#Compound D-I Optimal designs
# coor_exch_compound=function(x,x_matrix,C,model_order,
#                             me_index = c(),qe_index = c(),
#                             two_fi_index = c(),
#                             C_dtype, W, freeze_rows=0,
#                             x_matrix_Dopt,
#                             x_matrix_Iopt,w,K=0,
#                             D_opt_type="Bayesian",me_index_daug, qe_index_daug,two_fi_index_daug,
#                             telescoping=FALSE){
#   x_start=x
#   p=ncol(x_matrix)
#   #p_iopt=ncol(x_matrix_Iopt)
#   #computes the inverse of X'x
#   D= solve(t(x_matrix)%*%x_matrix+K)
#   n=nrow(x_matrix)
#   #saves a copy of D
#   D_local=D
#   dopt=det(t(x_matrix_Dopt)%*%x_matrix_Dopt + K)
#   iopt= sum(diag(solve(t(x_matrix_Iopt)%*%x_matrix_Iopt)%*%W))
#   #print(iopt)
#   #saves a local copy of x_matrix
#   x_matrix_local=x_matrix
#
#   #epsilon for monitoring convergence
#   epsilon=1
#
#   #Current L-optimality score
#   M=det(t(x_matrix_local)%*%x_matrix_local + K)
#   D_eff=(M/dopt)**(1/p)
#   if(D_opt_type=="Bayesian"){
#     x_matrix_local_I=t(apply(x,MARGIN=1,FUN=f_x,C_dtype=C_dtype,
#                              model_order=model_order,
#                              me_index = me_index,qe_index = qe_index,
#                              two_fi_index = two_fi_index))
#     print(x_matrix_local_I)
#     D_local_I=solve(t(x_matrix_local_I)%*%x_matrix_local_I)
#     I_val=sum(diag(D_local_I%*%W))
#   }else{
#     I_val=sum(diag(D_local%*%W))
#   }
#   #I_eff=(iopt/I_val)**(1/p_iopt)
#   I_eff=(iopt/I_val)
#   current=w*D_eff+(1-w)*I_eff
#   print(paste("starting objective:",current))
#   count_convergence=0
#   while(epsilon>1e-3){
#     count_convergence=count_convergence+1
#     # at the end of one complete exhchange of x_matrix_local, the updated  matrix is stored
#     x_matrix=x_matrix_local
#     prev=current
#     ieff_best=I_eff
#     deff_best=D_eff
#     x_best=x
#     for(i in (freeze_rows+1):n){
#       #print(paste("exchange of row number:", i))
#       for(j in 1:length(C)){
#         #print(paste("exchange of col number:", j))
#         det_val=c()
#         #current value of trace
#         #I_val_current=sum(diag(D_local%*%W))
#         #D_val_current=M
#         if(telescoping){
#           point=x_start[i,j]
#           candidate_points=seq(max(-1,point-0.09),min(1,point+0.09),0.01)
#         }else{
#           candidate_points=C[[j]]
#         }
#         for(k in candidate_points){
#           x_local=x[i,]
#           #point exchange
#           x_local[j]=k
#           if(D_opt_type=="Bayesian"){
#             x_expl_I=f_x(x_local,model_order,C_dtype,
#                          me_index = me_index,qe_index = qe_index,
#                          two_fi_index = two_fi_index)
#             #x_new_I=x_matrix_local_I
#             #x_new_I[i,]=x_expl_I
#             D_new_I=rank2_inverse(x_expl_I,x_matrix_local_I[i,],D_local_I)
#             x_expl=f_x(x_local,model_order="non_standard",C_dtype,
#                        me_index = me_index_daug,qe_index = qe_index_daug,
#                        two_fi_index = two_fi_index_daug)
#             #x_new=x_matrix_local
#             #x_new[i,]=x_expl
#             D_new=rank2_inverse(x_expl,x_matrix_local[i,],D_local)
#           }else{
#             x_expl=f_x(x_local,model_order,C_dtype,
#                        me_index = me_index,qe_index = qe_index,
#                        two_fi_index = two_fi_index)
#             #x_new=x_matrix_local
#             #x_new[i,]=x_expl
#             #compute the D_inverse after exchange
#             D_new=rank2_inverse(x_expl,x_matrix_local[i,],D_local)
#             D_new_I=0
#           }
#           #check if singular
#           if(D_new[1]!="Singular" && D_new_I[1]!="Singular"){
#             #M_local=det(t(x_new)%*%x_new)#M*delta_D(x_expl,x_matrix_local[i,],D_local)
#             M_local=M*delta_D(x_expl,x_matrix_local[i,],D_local)
#             if(D_opt_type=="Bayesian"){
#               I_val_local=sum(diag(D_new_I%*%W))
#             }else{
#               I_val_local=sum(diag(D_new%*%W))
#             }
#             #local_val=w*((M_local/dopt)**(1/p))+(1-w)*(iopt/I_val_local)**(1/p_iopt)
#             if(M_local>0){
#               local_val=w*((M_local/dopt)**(1/p))+(1-w)*(iopt/I_val_local)
#               det_val=c(det_val,local_val)
#             }else{
#               det_val=c(det_val,NaN)
#             }
#           }else{
#             det_val=c(det_val,NaN)
#           }
#         }
#         #det_val=round(det_val,8)
#         if(any(is.na(det_val))){
#           print(paste("change in objective function at iteration",count_convergence))
#           print(det_val)
#         }
#         #det_val=ifelse(abs(det_val)==Inf,NaN,det_val)
#         #print(det_val)
#         #if(!any(is.na(det_val)) & !any(abs(det_val)==Inf)){
#         #if all of de_val is not NA
#         if(!all(is.na(det_val))){
#           ind=which(det_val==max(det_val,na.rm=TRUE))[1]
#           #print(paste("max=",ind))
#           #make the swap only if the delta is greater than 0
#           if(det_val[ind]-current>0){
#             if(D_opt_type=="Bayesian"){
#               x[i,j]=candidate_points[ind]
#               M=M*delta_D(f_x(x[i,],model_order="non_standard",C_dtype,
#                               me_index = me_index_daug,qe_index = qe_index_daug,
#                               two_fi_index = two_fi_index_daug),x_matrix_local[i,],D_local)
#               D_local=rank2_inverse(f_x(x[i,],model_order="non_standard",C_dtype,
#                                         me_index = me_index_daug,qe_index = qe_index_daug,
#                                         two_fi_index = two_fi_index_daug),x_matrix_local[i,],D_local)
#               x_matrix_local[i,]=f_x(x[i,],model_order="non_standard",C_dtype,
#                                      me_index = me_index_daug,qe_index = qe_index_daug,
#                                      two_fi_index = two_fi_index_daug)
#               D_local_I=rank2_inverse(f_x(x[i,],model_order,C_dtype,
#                                           me_index = me_index,qe_index = qe_index,
#                                           two_fi_index = two_fi_index),x_matrix_local_I[i,],D_local_I)
#               x_matrix_local_I[i,]=f_x(x[i,],model_order,C_dtype,
#                                        me_index = me_index,qe_index = qe_index,
#                                        two_fi_index = two_fi_index)
#               I_val=sum(diag(D_local_I%*%W))
#               #print(I_val)
#               #current=w*((M_local/dopt)**(1/p))+(1-w)*(iopt/I_val)**(1/p_iopt)
#               I_eff=(iopt/I_val)
#               D_eff=(M/dopt)**(1/p)
#               current=w*(D_eff)+(1-w)*I_eff
#             }else{
#               x[i,j]=candidate_points[ind]
#               M=M*delta_D(f_x(x[i,],model_order,C_dtype,
#                               me_index = me_index,qe_index = qe_index,
#                               two_fi_index = two_fi_index),x_matrix_local[i,],D_local)
#               D_local=rank2_inverse(f_x(x[i,],model_order,C_dtype,
#                                         me_index = me_index,qe_index = qe_index,
#                                         two_fi_index = two_fi_index),x_matrix_local[i,],D_local)
#               x_matrix_local[i,]=f_x(x[i,],model_order,C_dtype,
#                                      me_index = me_index,qe_index = qe_index,
#                                      two_fi_index = two_fi_index)
#               I_val=sum(diag(D_local%*%W))
#               #print(I_val)
#               #current=w*((M_local/dopt)**(1/p))+(1-w)*(iopt/I_val)**(1/p_iopt)
#               I_eff=(iopt/I_val)
#               D_eff=(M/dopt)**(1/p)
#               current=w*(D_eff)+(1-w)*I_eff
#             }
#           }#else{
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
#   return(list(x_matrix,prev,x_best,deff_best,ieff_best))
# }

#' function for initiating multiple starts of the Adaptive-RSO designs using coordinate exchange
#' @param x_matrix_start the design matrix of the specified order including the intercept of the random starting
#' @param x_start the settings matrix of the random starting design
#' @param W_mat Moment matrix computed for a given lambda
#' @param n total number of runs in the design
#' @param C list of candidate space for each experimental factor in the design
#' @param model_order a string indicating the order of the regression model. "linear","quadratic" (Main Effects+Quadratic Effects), "2FI" (Main Effects + TwoFactor Interactions),"full", "non_standard"- for a user defined model order order
#' @param C_dtype a string vector indicating the data type of the experimental factors in the design either "cont"(continuous) or "cat"(categorical)
#' @param type a string indicating the type of design augmentation. "non_augment" or "augment" (augment an existing design)
#' @param D_opt_design the D-optimal design matrix
#' @param I_opt_design the I-optimal design matrix
#' @param me_index vector of index of Main Effects to be included in the design for model_order "non_standard"
#' @param qe_index vector index of Quadratic Effects to be included in the design for model_order "non_standard"
#' @param two_fi_index list of index of Two Factor Interactions to be included in the design for model_order "non_standard"
#' @param K a diagonal matrix specifying the prior variance for Bayesian D-optimal design. By default set to 0.001 for potential terms and 0 for primary terms.
#' @param n_starts number of random starts, default at 100
#' @param w_candidate vector of weights to be used in the computation of compound optimal design
#' @param D_opt_type string indicating the type of D-optimal design Bayesian or regular D-optimal design. "Bayesian" or "non_Bayesian"
#' @description function for initiating multiple starts of a compound D-I design using coordinate exchange.This is the serial version of the code, by default a parallel version is used. The D-I optimal design is the exploration-exploitation trade-off design
#' @return return a list of best design matrix for each w in the w_candidate vector
#'  @examples
#'   construct_des_compound_serial=function(x_matrix_start,x_start,W_mat,n,C_dtype,C,model_order, type= "augment",D_opt_design,I_opt_design,me_index = c(),qe_index = c(),two_fi_index = c(), n_starts=100,w_candidate=seq(0.05,0.95,0.05),D_opt_type="non_Bayesian",K=0)
#' @export
construct_des_compound_serial=function(x_matrix_start,x_start,W_mat,n,C_dtype,C,
                                model_order, type= "augment",
                                D_opt_design,I_opt_design,
                                me_index = c(),qe_index = c(),
                                two_fi_index = c(), n_starts=100,
                                w_candidate=seq(0.05,0.95,0.05),
                                D_opt_type="non_Bayesian",
                                K=0){
  w_eff=list()
  junk=1
  w_candidate=w_candidate
  freeze_rows=nrow(x_matrix_start)
  for(w in w_candidate){
    best_compound=0
    for(i in 1:n_starts){
      print(paste("random start",i))
      singular=TRUE
      while(singular){
        if(D_opt_type=="Bayesian"){
          x=gen_x(n=n-freeze_rows,C,model_order="full",C_dtype,
                  type="augment",
                  me_index = c(),qe_index = c(),
                  two_fi_index = c(),design="Bayesian")
        }else{
          x=gen_x(n=n-freeze_rows,C,model_order=model_order,C_dtype,
                  type="augment",
                  me_index = me_index,qe_index = qe_index,
                  two_fi_index = two_fi_index)
        }
        #print(x)
        x_matrix=x[[2]]
        x_matrix=rbind(x_matrix_start,x_matrix)
        x=x[[1]]
        if(class(x)!="matrix"){x=matrix(x)}
        #print(x)
        #print(x_matrix)
        x=rbind(x_start,x)
        if(rcond(t(x_matrix)%*%x_matrix+K)>1e-5){
          singular=FALSE
        }
      }
      #print(x)
      #x_matrix=cnn_doe_mat
      #x=x_cnn_doe
      if(class(x)!="matrix"){
        x=matrix(x)
      }
      print(paste("Coordinate exchange for try",i))
      des=coor_exch_compound(x=x,x_matrix=x_matrix,C=C,model_order =model_order,
                             C_dtype=C_dtype,W=W_mat,x_matrix_Dopt = D_opt_design,
                             w=w,x_matrix_Iopt = I_opt_design,
                             freeze_rows = freeze_rows,
                             me_index = me_index,qe_index = qe_index,
                             two_fi_index = two_fi_index,
                             D_opt_type="Bayesian",
                             K=K)
      current=des[[2]]
      x_current=des[[3]]
      des=des[[1]]
      if(current>best_compound){
        best_compound=current
        best_des_compound=des
        x_best=x_current
        #break
      }
      print("------------------------------------")
    }
    w_eff[[junk]]=best_des_compound
    junk=junk+1
  }
  #best_des_compound
  return(w_eff)
}

#' function for initiating individual coordinate exchange for I optimal design for the parallel code.
#' @param i random start index
#' @param x_matrix_start the design matrix of the specified order including the intercept of the random starting
#' @param x_start the settings matrix of the random starting design
#' @param W_mat Moment matrix computed for a given lambda
#' @param n total number of runs in the design
#' @param C list of candidate space for each experimental factor in the design
#' @param model_order a string indicating the order of the regression model. "linear","quadratic" (Main Effects+Quadratic Effects), "2FI" (Main Effects + TwoFactor Interactions),"full", "non_standard"- for a user defined model order order
#' @param C_dtype a string vector indicating the data type of the experimental factors in the design either "cont"(continuous) or "cat"(categorical)
#' @param type a string indicating the type of design augmentation. "non_augment" or "augment" (augment an existing design)
#' @param me_index vector of index of Main Effects to be included in the design for model_order "non_standard"
#' @param qe_index vector index of Quadratic Effects to be included in the design for model_order "non_standard"
#' @param two_fi_index list of index of Two Factor Interactions to be included in the design for model_order "non_standard"
#' @param use_cpp a boolean value indicating whether to use the CPP version of coordinate exchange
#' @description function for initiating individual coordinate exchange for I optimal design for the parallel code.
#' @return a list of
#'  \itemize{
#'        \item x_matrix- a I-optimal design matrix of the specified order
#'        \item optimal_value- Optimal value of the objective value of the I-optimal design
#'        \item x- a matrix with factor levels of the I-optimal design
#'        }
#'  @examples
#'   individual_starts_i_opt=function(i,x_matrix_start,x_start,W_mat,n,C_dtype, C,model_order,type,me_index,qe_index,two_fi_index)
#' @export
individual_starts_i_opt=function(i,x_matrix_start,x_start,W_mat,n,C_dtype, C,
                                 model_order,type,
                                 me_index,qe_index,
                                 two_fi_index,use_cpp=TRUE){
  if(type != "augment"){
    x=gen_x(n=n,C=C,model_order=model_order,C_dtype=C_dtype,
            type="non_augment",
            me_index = me_index,qe_index = qe_index,
            two_fi_index = two_fi_index)
    x_matrix=x[[2]]
    x=x[[1]]
    if(class(x)!="matrix"){x=matrix(x)}
  }else{
    freeze_rows=nrow(x_matrix_start)
    singular=TRUE
    while(singular){
      x=gen_x(n=n-freeze_rows,C,model_order=model_order,C_dtype,
              type="augment",
              me_index = me_index,qe_index = qe_index,
              two_fi_index = two_fi_index)
      #print(x)
      x_matrix=x[[2]]
      print(dim(x_matrix))
      print(dim(x_matrix_start))
      x_matrix=rbind(x_matrix_start,x_matrix)
      x=x[[1]]
      if(class(x)!="matrix"){x=matrix(x)}
      #print(x)
      #print(x_matrix)
      x=rbind(x_start,x)
      if(rcond(t(x_matrix)%*%x_matrix)>1e-5){
        singular=FALSE
      }
    }
    if(class(x)!="matrix"){
      x=matrix(x)
    }
  }
  print(use_cpp)
  if(use_cpp){
    me_index_=me_index-1
    qe_index_=qe_index-1
    two_fi_index_=lapply(two_fi_index,FUN=function(z){z-1})
  }else{
    me_index_=me_index
    qe_index_=qe_index
    two_fi_index_=two_fi_index
  }
  print(paste("Coordinate exchange for try",i))
  des=coor_exch_I_augment(x=x,x_matrix=x_matrix,C=C,model_order =model_order,
                          C_dtype=C_dtype,W=W_mat,freeze_rows = freeze_rows,
                          me_index=me_index_,qe_index=qe_index_,
                          two_fi_index=two_fi_index_)
  print("Starting Telescoping")
  des=coor_exch_I_augment(x=des[[3]],x_matrix=des[[1]],C=C,model_order =model_order,
                          C_dtype=C_dtype,W=W_mat,freeze_rows = freeze_rows,
                          me_index=me_index_,qe_index=qe_index_,
                          two_fi_index=two_fi_index_,telescoping=TRUE)
  return(des)
}

#' function for initiating the coordinate exchange for I optimal design for a given number of random starts
#' @param x_matrix_start the design matrix of the specified order including the intercept of the random starting
#' @param x_start the settings matrix of the random starting design
#' @param W_mat Moment matrix computed for a given lambda
#' @param n total number of runs in the design
#' @param C list of candidate space for each experimental factor in the design
#' @param model_order a string indicating the order of the regression model. "linear","quadratic" (Main Effects+Quadratic Effects), "2FI" (Main Effects + TwoFactor Interactions),"full", "non_standard"- for a user defined model order order
#' @param C_dtype a string vector indicating the data type of the experimental factors in the design either "cont"(continuous) or "cat"(categorical)
#' @param type a string indicating the type of design augmentation. "non_augment" or "augment" (augment an existing design)
#' @param me_index vector of index of Main Effects to be included in the design for model_order "non_standard"
#' @param qe_index vector index of Quadratic Effects to be included in the design for model_order "non_standard"
#' @param two_fi_index list of index of Two Factor Interactions to be included in the design for model_order "non_standard"
#' @param n_starts number of random starts, default at 100
#' @param rng seed for random number generation
#' @param use_cpp a boolean value indicating whether to use the CPP version of coordinate exchange
#' @param num_cores number of cores for parallel processing. Default value is NA, in which case Sys.getenv("PBS_RESC_TOTAL_PROCS") is used to get the number of cores in a HPC. Otherwise, the user needs to provide a value. Use 1 for serial execution.
#' @description function for initiating the coordinate exchange for I optimal design a given number of random starts and returns the best design. This is the parallel version of the code and works only on linux machines.
#' @return a list of
#'  \itemize{
#'        \item best_I_opt_design- the best I-optimal design matrix of the specified order from the random starts
#'        \item best_Iopt_x- a matrix of factor levels for the best I-optimal design amongst the random starts
#'        \item best_I_obj- objective value of the best I-optimal design from the random starts
#'        }
#'  @examples
#'   construct_des_I_lambda(x_matrix_start,x_start,W_mat,n,C_dtype, C, model_order, type= "augment", me_index = c(),qe_index = c(),two_fi_index = c(), n_starts=100)
#' @import doRNG
#' @export
construct_des_I_lambda=function(x_matrix_start,x_start,W_mat,n,C_dtype, C,
                                model_order, type= "augment",
                                me_index = c(),qe_index = c(),
                                two_fi_index = c(), n_starts=100, rng=c(),use_cpp=TRUE,num_cores=NA){
  #best/sum(diag(x_matrix_start%*%W_mat))
  #x_start=x_dsd
  best_augment=1e8
  #cl=makeCluster(detectCores()-1, type="FORK")
  if(is.na(num_cores)){
    cl=parallel::makeCluster(as.numeric(Sys.getenv("PBS_RESC_TOTAL_PROCS")) - 1, type="FORK",outfile="")
  }else{
    cl=parallel::makeCluster(num_cores, type="FORK",outfile="")
  }
  doParallel::registerDoParallel(cl)
  #registerDoRNG(seed = 1234)
  #registerDoParallel(cores=detectCores())
  #registerDoParallel(cores=as.numeric(Sys.getenv("PBS_RESC_TOTAL_PROCS")) - 1)
  mcoptions = list(preschedule=FALSE)
  chunks = foreach::getDoParWorkers()
  print(paste("Number of workers is:",chunks))
  multistart_results=foreach::foreach(i=1:n_starts,.options.multicore=mcoptions,.options.RNG=rng) %dorng% individual_starts_i_opt(i,x_matrix_start,x_start,W_mat,n,C_dtype, C,
                                                                                                                          model_order,type,
                                                                                                                          me_index,qe_index,
                                                                                                                          two_fi_index,use_cpp = TRUE)
  #stopImplicitCluster()
  parallel::stopCluster(cl)
  print(paste("Length of output is:",length(multistart_results)))

  for(i in 1:length(multistart_results)){
    current=multistart_results[[i]][[2]]
    x_current=multistart_results[[i]][[3]]
    des=multistart_results[[i]][[1]]
    if(current<best_augment){
      best_augment=current
      best_des_augment=des
      x_best_augment=x_current
    }
  }
  number_best_converegence=0
  criteria_list=c()
  for(i in 1:length(multistart_results)){
    criteria_list=c(criteria_list,multistart_results[[i]][[2]])
    if(isTRUE(all.equal(best_augment,multistart_results[[i]][[2]]))){
      if(all(x_best_augment==multistart_results[[i]][[3]])){
        number_best_converegence=number_best_converegence+1
      }
    }
  }
  print(paste("Number of iterations converged to the best value is",number_best_converegence))
  graphics::plot(x=1:length(multistart_results),y=criteria_list, main="Plot of compound criteria", xlab="Iteration",ylab="Criteria_val")
  graphics::abline(h=best_augment, col="red")
  #print("Top 50 results")
  #print(criteria_list[order(criteria_list,decreasing=FALSE)[1:50]])
  return(list("best_I_opt_design"=best_des_augment,"best_Iopt_x"=x_best_augment,"best_I_obj"=best_augment))
}

#' function for initiating multiple starts of Adaptive-RSO designs for each random start
#' @param i index of the random start
#' @param x_matrix_start the design matrix of the specified order including the intercept of the random starting
#' @param x_start the settings matrix of the random starting design
#' @param W_mat Moment matrix computed for a given lambda
#' @param n total number of runs in the design
#' @param C list of candidate space for each experimental factor in the design
#' @param model_order a string indicating the order of the regression model. "linear","quadratic" (Main Effects+Quadratic Effects), "2FI" (Main Effects + TwoFactor Interactions),"full", "non_standard"- for a user defined model order order
#' @param C_dtype a string vector indicating the data type of the experimental factors in the design either "cont"(continuous) or "cat"(categorical)
#' @param type a string indicating the type of design augmentation. "non_augment" or "augment" (augment an existing design)
#' @param freeze_rows the rows of the matrix to freeze while augmenting, 0 when type is "non_augment"
#' @param D_opt_design the D-optimal design matrix
#' @param D_opt_x the settings matrix of the D-optimal design
#' @param I_opt_design the I-optimal design matrix
#' @param me_index vector of index of Main Effects to be included in the design for model_order "non_standard"
#' @param qe_index vector index of Quadratic Effects to be included in the design for model_order "non_standard"
#' @param two_fi_index list of index of Two Factor Interactions to be included in the design for model_order "non_standard"
#' @param me_index_daug vector of index of Main Effects to be included in the Bayes-D-optimal design for model_order "non_standard"
#' @param qe_index_daug vector index of Quadratic Effects to be included in the Bayes-D-optimal design for model_order "non_standard"
#' @param two_fi_index_daug list of index of Two Factor Interactions to be included in the Bayes-D-optimal design for model_order "non_standard"
#' @param K a diagonal matrix specifying the prior variance in Bayesian D-optimal design set to 0.001 for potential terms and 0 for primary terms
#' @param D_opt_type Type of D-optimal design Bayesian or regular D-optimal design. "Bayesian" or "non_Bayesian"
#' @param use_cpp a boolean value indicating whether to use the CPP version of coordinate exchange
#' @description function for initiating multiple starts of the Adaptive-RSO designs using coordinate exchange algorithm for each random start
#' @return returns a list of compound optimal design and its optimal value for a given w and random start
#'   @examples
#'   individual_starts_compound(i,x_matrix_start,x_start,W_mat,n,C_dtype,C,model_order, type,freeze_rows,D_opt_design,D_opt_x,I_opt_design,me_index,qe_index,two_fi_index,w,D_opt_type,K,me_index_daug, qe_index_daug,two_fi_index_daug,use_cpp=TRUE)
#' @export
individual_starts_compound=function(i,x_matrix_start,x_start,W_mat,n,C_dtype,C,
                                    model_order, type,freeze_rows,
                                    D_opt_design,D_opt_x,I_opt_design,
                                    me_index,qe_index,
                                    two_fi_index,
                                    w,D_opt_type,K,
                                    me_index_daug, qe_index_daug,two_fi_index_daug,use_cpp=TRUE){
  if(i<2){
    print("Compound Optimal with Bayesian D-Optimal Design as starting design")
    x_matrix=D_opt_design
    x=D_opt_x
  }else{
    singular=TRUE
    while(singular){
      if(D_opt_type=="Bayesian"){
        x=gen_x(n=n-freeze_rows,C,model_order="non_standard",C_dtype,
                type="augment",
                me_index = me_index_daug,qe_index = qe_index_daug,
                two_fi_index = two_fi_index_daug,design="Bayesian")
      }else{
        x=gen_x(n=n-freeze_rows,C,model_order=model_order,C_dtype,
                type=type,
                me_index = me_index,qe_index = qe_index,
                two_fi_index = two_fi_index)
      }
      #print(x)
      x_matrix=x[[2]]
      #print(dim(x_matrix))
      #print(dim(x_matrix_start))
      x_matrix=rbind(x_matrix_start,x_matrix)
      x=x[[1]]
      if(class(x)!="matrix"){x=matrix(x)}
      #print(x)
      #print(x_matrix)
      x=rbind(x_start,x)
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
  }
  print(use_cpp)
  if(use_cpp){
    me_index_daug_=me_index_daug-1
    qe_index_daug_=qe_index_daug-1
    two_fi_index_daug_=lapply(two_fi_index_daug,FUN=function(z){z-1})
    me_index_=me_index-1
    qe_index_=qe_index-1
    two_fi_index_=lapply(two_fi_index,FUN=function(z){z-1})
  }else{
    me_index_daug_=me_index_daug
    qe_index_daug_=qe_index_daug
    two_fi_index_daug_=two_fi_index_daug
    me_index_=me_index
    qe_index_=qe_index
    two_fi_index_=two_fi_index
  }
  print(paste("Coordinate exchange for try",i))
  des=coor_exch_compound(x=x,x_matrix=x_matrix,C=C,model_order =model_order,
                         C_dtype=C_dtype,W=W_mat,x_matrix_Dopt = D_opt_design,
                         w=w,x_matrix_Iopt = I_opt_design,
                         freeze_rows = freeze_rows,
                         me_index = me_index_,qe_index = qe_index_,
                         two_fi_index = two_fi_index_,
                         D_opt_type=D_opt_type,
                         K=K,
                         me_index_daug=me_index_daug_,
                         qe_index_daug=qe_index_daug_,
                         two_fi_index_daug=two_fi_index_daug_)
  print("Starting telescoping")
  des=coor_exch_compound(x=des[[3]],x_matrix=des[[1]],C=C,model_order =model_order,
                         C_dtype=C_dtype,W=W_mat,x_matrix_Dopt = D_opt_design,
                         w=w,x_matrix_Iopt = I_opt_design,
                         freeze_rows = freeze_rows,
                         me_index = me_index_,qe_index = qe_index_,
                         two_fi_index = two_fi_index_,
                         D_opt_type=D_opt_type,
                         K=K,
                         me_index_daug=me_index_daug_,
                         qe_index_daug=qe_index_daug_,
                         two_fi_index_daug=two_fi_index_daug_,
                         telescoping=TRUE)
  return(des)
}

#' function for initiating multiple starts of the Adaptive-RSO
#' @param x_matrix_start the design matrix of the specified order including the intercept of the random starting
#' @param x_start the settings matrix of the random starting design
#' @param W_mat Moment matrix computed for a given lambda
#' @param n total number of runs in the design
#' @param C list of candidate space for each experimental factor in the design
#' @param model_order a string indicating the order of the regression model. "linear","quadratic" (Main Effects+Quadratic Effects), "2FI" (Main Effects + TwoFactor Interactions),"full", "non_standard"- for a user defined model order order
#' @param C_dtype a string vector indicating the data type of the experimental factors in the design either "cont"(continuous) or "cat"(categorical)
#' @param type a string indicating the type of design augmentation. "non_augment" or "augment" (augment an existing design)
#' @param D_opt_design the D-optimal design matrix
#' @param D_opt_x the settings matrix of the D-optimal design
#' @param I_opt_design the I-optimal design matrix
#' @param me_index vector of index of Main Effects to be included in the I-optimal design for model_order "non_standard"
#' @param qe_index vector index of Quadratic Effects to be included in the I-optimal design for model_order "non_standard"
#' @param two_fi_index list of index of Two Factor Interactions to be included in the I-optimal design for model_order "non_standard"
#' @param me_index_daug vector of index of Main Effects to be included in the Bayes-Doptimal design for model_order "non_standard"
#' @param qe_index_daug vector index of Quadratic Effects to be included in the Bayes-Doptimal design for model_order "non_standard"
#' @param two_fi_index_daug list of index of Two Factor Interactions to be included in the Bayes-Doptimal design for model_order "non_standard"
#' @param n_starts number of random starts, default at 100
#' @param w_candidate vector of weights to be used in the computation of compound optimal design
#' @param K a diagonal matrix specifying the prior variance in Bayesian D-optimal design set to 0.001 for potential terms and 0 for primary terms
#' @param D_opt_type Type of D-optimal design Bayesian or regular D-optimal design. "Bayesian" or "non_Bayesian"
#' @param use_cpp a boolean value indicating whether to use the CPP version of coordinate exchange
#' @param rng random seed
#' @param num_cores number of cores for parallel processing. Default value is NA, in which case Sys.getenv("PBS_RESC_TOTAL_PROCS") is used to get the number of cores in a HPC. Otherwise, the user needs to provide a value. Use 1 for serial execution.
#' @description function for initiating multiple starts of Adaptove-RSO designs using coordinate exchange algorithm.This is the parallel version and works only on linux machines. The Adaptive-RSO designs are compound optimal designs between Bayesian D-optimal designs and I-optimal designs.
#' @return return a list of best matrix for random starts for each w in the candidate list
#'   @examples
#'   construct_des_compound_wcand(x_matrix_start,x_start,W_mat,n,C_dtype,C, model_order, type= "augment",D_opt_design,D_opt_x,I_opt_design,me_index = c(),qe_index = c(),two_fi_index = c(), n_starts=100,w_candidate=seq(0.05,0.95,0.05),D_opt_type="non_Bayesian",K=0,me_index_daug, qe_index_daug,two_fi_index_daug,use_cpp=TRUE)
#' @import doRNG
#' @export
construct_des_compound_wcand=function(x_matrix_start,x_start,W_mat,n,C_dtype,C,
                                model_order, type= "augment",
                                D_opt_design,D_opt_x,I_opt_design,
                                me_index = c(),qe_index = c(),
                                two_fi_index = c(), n_starts=500,
                                w_candidate=seq(0.05,0.95,0.05),
                                D_opt_type="non_Bayesian",rng=c(),
                                K=0,me_index_daug, qe_index_daug,two_fi_index_daug,use_cpp=TRUE,num_cores=NA){
  w_eff=vector(mode="list",length=length(w_candidate))
  #cont_index=which(C_dtype=="cont")
  #for(i in cont_index){
  #  C[[i]]=seq(-1,1,0.01)
  #}
  freeze_rows=nrow(x_matrix_start)
  print(paste("use_cpp",use_cpp))
  for(w in 1:length(w_candidate)){
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
    #chunks = foreach::getDoParWorkers()
    #print(paste("Number of workers is:",chunks))
    multistart_results=foreach::foreach(i=1:n_starts, .options.multicore=mcoptions,.options.RNG=rng) %dorng% individual_starts_compound(i,x_matrix_start,x_start,W_mat,n,C_dtype,C,
                                                                                                                                model_order, type, freeze_rows=freeze_rows,
                                                                                                                                D_opt_design,D_opt_x,I_opt_design,
                                                                                                                                me_index,qe_index,
                                                                                                                                two_fi_index,
                                                                                                                                w=w_candidate[w],D_opt_type,K,
                                                                                                                                me_index_daug, qe_index_daug,two_fi_index_daug,use_cpp=use_cpp)
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
        ieff_best=multistart_results[[i]][[5]]
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
    w_eff[[as.character(w_candidate[w])]]=list(best_des_compound,c("Weighted criteria"=best_compound,"D_eff"=deff_best,"I_eff"=ieff_best))
    print(paste("Number of iterations converged to the best value is",number_best_converegence))
    graphics::plot(x=1:length(multistart_results),y=criteria_list, main=paste("Plot of compound criteria",w_candidate[w]), xlab="Iteration",ylab="Criteria_val")
    graphics::abline(h=best_compound, col="red")
    print("Top 50 results")
    print(criteria_list[order(criteria_list,decreasing=TRUE)[1:50]])
  }
  x_matrix_I=t(apply(D_opt_x,MARGIN=1,FUN=f_x,C_dtype=C_dtype,
                     model_order="non_standard",
                     me_index = me_index,qe_index = qe_index,
                     two_fi_index = two_fi_index))
  D_I=solve(t(x_matrix_I)%*%x_matrix_I)
  I_val=sum(diag(D_I%*%W_mat))
  iopt=sum(diag(solve(t(I_opt_design)%*%I_opt_design)%*%W_mat))
  w_eff[["1"]]=list(best_des_compound,c("Weighted criteria"=1,"D_eff"=1,"I_eff"=iopt/I_val))
  #print(w_eff)
  return(w_eff)
}


#' function for initiating multiple starts of the Adaptive-RSO designs
#' @param x_matrix_start the design matrix of the specified order including the intercept of the random starting
#' @param x_start the settings matrix of the random starting design
#' @param W_mat Moment matrix computed for a given lambda
#' @param n total number of runs in the design
#' @param C list of candidate space for each experimental factor in the design
#' @param model_order a string indicating the order of the regression model. "linear","quadratic" (Main Effects+Quadratic Effects), "2FI" (Main Effects + TwoFactor Interactions),"full", "non_standard"- for a user defined model order order
#' @param C_dtype a string vector indicating the data type of the experimental factors in the design either "cont"(continuous) or "cat"(categorical)
#' @param type a string indicating the type of design augmentation. "non_augment" or "augment" (augment an existing design)
#' @param D_opt_design the D-optimal design matrix
#' @param D_opt_x the settings matrix of the D-optimal design
#' @param I_opt_design the I-optimal design matrix
#' @param me_index vector of index of Main Effects to be included in the I-optimal design for model_order "non_standard"
#' @param qe_index vector index of Quadratic Effects to be included in the I-optimal design for model_order "non_standard"
#' @param two_fi_index list of index of Two Factor Interactions to be included in the I-optimal design for model_order "non_standard"
#' @param me_index_daug vector of index of Main Effects to be included in the Bayes-D-optimal design for model_order "non_standard"
#' @param qe_index_daug vector index of Quadratic Effects to be included in the Bayes-D-optimal design for model_order "non_standard"
#' @param two_fi_index_daug list of index of Two Factor Interactions to be included in the Bayes-D-optimal design for model_order "non_standard"
#' @param n_starts number of random starts, default at 100
#' @param K a diagonal matrix specifying the prior variance in Bayesian D-optimal design set to 0.001 for potential terms and 0 for primary terms
#' @param D_opt_type Type of D-optimal design Bayesian or regular D-optimal design. "Bayesian" or "non_Bayesian"
#' @param use_cpp a boolean value indicating whether to use the CPP version of coordinate exchange
#' @param D_opt_tresh a numeric value providing the lower bound of the threshold in bisection search for the D-optimal criteria. Default value is 0.8
#' @param I_opt_tresh a numeric value providing the lower bound of the threshold in bisection search for the I-optimal criteria. Default value is 0.8
#' @param w_tresh the threshold for difference in weights w for the bisection search to converge. The default value is 0.02
#' @param rng seed for random number generation for the random design construction starts
#' @param num_cores number of cores for parallel processing. Default value is NA, in which case Sys.getenv("PBS_RESC_TOTAL_PROCS") is used to get the number of cores in a HPC. Otherwise, the user needs to provide a value. Use 1 for serial execution.
#' @description function for initiating multiple starts of an Adaptive-RSO designs using coordinate exchange.This is the parallel version and works only on linux machines. The Adaptive-RSO designs are compound optimal designs between Bayesian D-optimal designs and I-optimal designs and we use bisection search algorithm to determine the optimal weight w.
#' @return return a list of best matrix for random starts for each w in the candidate list
#'   @examples
#'   construct_des_compound(x_matrix_start,x_start,W_mat,n,C_dtype,C, model_order, type= "augment",D_opt_design,D_opt_x,I_opt_design,me_index = c(),qe_index = c(),two_fi_index = c(), n_starts=100,D_opt_type="non_Bayesian",K=0,me_index_daug, qe_index_daug,two_fi_index_daug,use_cpp=TRUE, D_opt_tresh=0.8, I_opt_tresh=0.8)
#' @import doRNG
#' @export

construct_des_compound=function(x_matrix_start,x_start,W_mat,n,C_dtype,C,
                                model_order, type= "augment",
                                D_opt_design,D_opt_x,I_opt_design,
                                me_index = c(),qe_index = c(),
                                two_fi_index = c(), n_starts=100,
                                #w_candidate=seq(0.05,0.95,0.05),
                                D_opt_type="non_Bayesian",
                                K=0,me_index_daug, qe_index_daug,two_fi_index_daug,use_cpp=TRUE,
                                D_opt_tresh=0.8, I_opt_tresh=0.8, w_tresh=0.02, rng=c(),num_cores=NA){
  w_eff=list()
  #cont_index=which(C_dtype=="cont")
  #for(i in cont_index){
  #  C[[i]]=seq(-1,1,0.01)
  #}
  freeze_rows=nrow(x_matrix_start)
  continue=TRUE
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
    multistart_results=foreach::foreach(i=1:n_starts, .options.multicore=mcoptions,.options.RNG=rng) %dorng% individual_starts_compound(i,x_matrix_start,x_start,W_mat,n,C_dtype,C,
                                                                                                                                model_order, type, freeze_rows=freeze_rows,
                                                                                                                                D_opt_design,D_opt_x,I_opt_design,
                                                                                                                                me_index,qe_index,
                                                                                                                                two_fi_index,
                                                                                                                                w=w,D_opt_type,K,
                                                                                                                                me_index_daug, qe_index_daug,two_fi_index_daug,use_cpp=use_cpp)
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
        ieff_best=multistart_results[[i]][[5]]
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
    w_eff[[as.character(w)]]=list(x_best,best_des_compound,c("Weighted criteria"=best_compound,"D_eff"=deff_best,"I_eff"=ieff_best))
    print(paste("Number of iterations converged to the best value is",number_best_converegence))
    graphics::plot(x=1:length(multistart_results),y=criteria_list, main=paste("Plot of compound criteria",w), xlab="Iteration",ylab="Criteria_val")
    graphics::abline(h=best_compound, col="red")
    print("Top 50 results")
    print(criteria_list[order(criteria_list,decreasing=TRUE)[1:50]])
    #update w_start and w_end
    print(paste("W start",w_start))
    print(paste("W end", w_end))
    print(paste(w))
    print(paste("I efficiency is,",ieff_best))
    print(paste("D efficiency is,",deff_best))
    if(deff_best>=D_opt_tresh & ieff_best>=I_opt_tresh){
      w_end=w
      w_start=w_start
      w=(w_start+w_end)/2
    }else if(deff_best>=D_opt_tresh & ieff_best<I_opt_tresh){
      w_end=w
      w_start=w_start
      w=(w_start+w_end)/2
    }else if(deff_best<D_opt_tresh & ieff_best>=I_opt_tresh){
      w_start=w
      w_end=w_end
      w=(w_start+w_end)/2
    }else{
      w=w_candidate[1]
      w_candidate=w_candidate[-1]
    }
    if(round(w_end-w_start,2)<=w_tresh | length(w_candidate)==0 | iter_w>10){
      continue=FALSE
    }
  }
  x_matrix_I=t(apply(D_opt_x,MARGIN=1,FUN=f_x,C_dtype=C_dtype,
                     model_order="non_standard",
                     me_index = me_index,qe_index = qe_index,
                     two_fi_index = two_fi_index))
  D_I=solve(t(x_matrix_I)%*%x_matrix_I)
  I_val=sum(diag(D_I%*%W_mat))
  iopt=sum(diag(solve(t(I_opt_design)%*%I_opt_design)%*%W_mat))
  w_eff[["1"]]=list(D_opt_design,c("Weighted criteria"=1,"D_eff"=1,"I_eff"=iopt/I_val))
  #print(w_eff)
  return(w_eff)
}
