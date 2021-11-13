#' function for constructing a Bayesian D-optimal design using coordinate exchange
#' @param x a settings matrix of the design
#' @param x_matrix design matrix of the specified order
#' @param C list of candidate space for each experimental factor in the design
#' @param model_order a string indicating the order of the regression model. "linear","quadratic" (Main Effects+Quadratic Effects), "2FI" (Main Effects + TwoFactor Interactions),"full", "non_standard"- for a user defined model order order
#' @param C_dtype a string vector indicating the data type of the experimental factors in the design either "cont"(continuous) or "cat"(categorical)
#' @param freeze_rows the rows of the matrix to freeze while augmenting, 0 when type is "non_augment"
#' @param me_index vector of index of Main Effects to be included in the design for model_order "non_standard"
#' @param qe_index vector of index of Quadratic Effects to be included in the design for model_order "non_standard"
#' @param two_fi_index list of index of Two Factor Interactions to be included in the design for model_order "non_standard"
#' @param cubic_index vector of index of Cubic Effects to be included in the design for model_order "non_standard"
#' @param quatric_index vector of index of Quatric Effects to be included in the design for model_order "non_standard"
#' @param type a string indicating the type of design augmentation. "non_augment" or "augment" (augment an existing deisgn)
#' @param design a string indicating the type of design constructed. "Non-Bayesian" or "Bayesian" optimal design. The major difference is in how signularity of the matrix is evaluated. In a Bayesian optimal design X'X+K(prior variance matrix) is checked for singularity
#' @param telescoping a boolean value taking TRUE or FALSE. Telescoping reduces the candidate factors the local region of the optimum design
#' @param K a diagonal matrix specifying the prior variance for Bayesian D-optimal design. By default set to 0.001 for potential terms and 0 for primary terms.
#' @description function for constructing a Bayesian D-optimal design using coordinate exchange
#' @return a list of
#'  \itemize{
#'        \item x_matrix- a Bayesian D-optimal design matrix of the specified order
#'        \item optimal_value- Optimal value of the objective value of the Bayesian D-optimal design
#'        \item x- a matrix with factor levels of the Bayesian D-optimal design
#'        }
#'  @export
coor_exch_D_Bayes_R= function(x,x_matrix,C,model_order,C_dtype,
                            freeze_rows = freeze_rows,
                            me_index=c(),qe_index=c(),
                            two_fi_index=c(),
                            cubic_index = c(),
                            quatric_index = c(),
                            type="non_augment",telescoping=FALSE,
                            K){
  x_start=x
  #comptes the starting values of D-optimality
  D= solve(t(x_matrix)%*%x_matrix + K)
  M=det(t(x_matrix)%*%x_matrix + K)
  n=nrow(x_matrix)
  #if type=="non_augment" change freeze_rows to 0
  if(type=="non_augment"){
    freeze_rows=0
  }
  v_ij=matrix(ncol=(n-freeze_rows))

  for(i in (freeze_rows+1):n){
    v_ij=rbind(v_ij,sapply(seq_along((freeze_rows+1):n),
                           FUN=function(j){v_x(x_i=x_matrix[i,],x_j=x_matrix[j,],D=D)}))
  }
  v_ij=v_ij[-1,]
  v_i=diag(v_ij)
  ord=order(x=v_i, decreasing = FALSE)
  if(type=="augment"){
    ord=freeze_rows+ord
  }
  x_matrix_local=x_matrix
  D_local=D

  epsilon=1
  #current=det(D_local)**(1/ncol(D))
  #current=M**(1/ncol(D))
  current=log(M)
  print(current)
  #print(ord)
  trial=1
  while(epsilon>1e-6){
    x_matrix=x_matrix_local
    prev=current
    x_final=x
    for(i in ord){
      for(j in 1:length(C)){
        det_val=c()
        #if telescoping is TRUE a local region around the current design is the candidate factor for each experimental factor
        if(telescoping){
          point=x_start[i,j]
          candidate_points=seq(max(-1,point-0.09),min(1,point+0.09),0.01)
        }else{
          candidate_points=C[[j]]
        }
        for(k in candidate_points){
          #coordinate exchange step
          x_local=x[i,]
          x_local[j]=k
          x_expl=f_x(x=x_local,model_order=model_order,C_dtype=C_dtype,
                     me_index=me_index,qe_index=qe_index,
                     two_fi_index=two_fi_index,
                     cubic_index = cubic_index,
                     quatric_index = quatric_index)
          det_val=c(det_val,delta_D(x_expl,x_matrix_local[i,],D_local))
        }
        ind=which(det_val==max(det_val))[1]
        x_local[j]=candidate_points[ind]
        x_expl=f_x(x=x_local,model_order=model_order,
                   C_dtype=C_dtype,
                   me_index=me_index,
                   qe_index=qe_index,
                   two_fi_index=two_fi_index,
                   cubic_index = cubic_index,
                   quatric_index = quatric_index)
        #compute the rank2inverse for the new matrix, returns singular if the exchange leads to a singular matrix
        D_local_temp=rank2_inverse(x_expl,
                                   x_matrix_local[i,],
                                   D=D_local)
        #the exchange is carried out only if it leads to a non singular matrix
        if(class(D_local_temp)!="character"){
          D_local=D_local_temp
          x[i,]=x_local
          x_matrix_local[i,]=x_expl
          #update the determinant
          M=M*det_val[ind]
          #print(M)
          #print(det(t(x_matrix_local)%*%x_matrix_local+K))
          #print(round(M,5)==round(det(t(x_matrix_local)%*%x_matrix_local+K),5))
          #D_local=solve(t(x_matrix_local)%*%x_matrix_local)
        }
      }
    }
    #current=det(D_local)**(1/ncol(D))
    #current=M**(1/ncol(D))
    current=log(M)
    if(is.nan(current)|!is.finite(current)){
      print(paste("The current value of M is:",M))
    }
    #current=det(t(x_matrix_local)%*%x_matrix_local)
    #epsilon=prev-current
    #computes the change in objective value
    epsilon=current-prev
    print(current)
    print(paste("Change in objective value is:",epsilon))
    if(epsilon<1e-6){
      trial=trial+1
      if(trial>4){
        print("Early stopping the random trial early since the change in objective is too small for more than 4 iterations")
        epsilon=0
      }
    }else{
      trial=1
    }

    v_ij=matrix(ncol=(n-freeze_rows))

    for(i in (freeze_rows+1):n){
      v_ij=rbind(v_ij,sapply(seq_along((freeze_rows+1):n),
                             FUN=function(j){v_x(x_i=x_matrix[i,],x_j=x_matrix[j,],D=D)}))
    }
    v_ij=v_ij[-1,]
    v_i=diag(v_ij)
    ord=order(x=v_i, decreasing = FALSE)
    if(type=="augment"){
      ord=freeze_rows+ord
    }
  }
  return(list("x_matrix"=x_matrix,"optimal_value"=prev,"x"=x_final))
}

#' function for initiating the coordinate exchange for a D-optimal design for a given number of random starts
#' @param x_matrix_start the design matrix of the specified order including the intercept of the random starting
#' @param x_start the settings matrix of the random starting design
#' @param n total number of runs in the design
#' @param C list of candidate space for each experimental factor in the design
#' @param model_order a string indicating the order of the regression model. "linear","quadratic" (Main Effects+Quadratic Effects), "2FI" (Main Effects + TwoFactor Interactions),"full", "non_standard"- for a user defined model order order
#' @param C_dtype a string vector indicating the data type of the experimental factors in the design either "cont"(continuous) or "cat"(categorical)
#' @param type a string indicating the type of design augmentation. "non_augment" or "augment" (augment an existing deisgn)
#' @param freeze_rows the rows of the matrix to freeze while augmenting, 0 when type is "non_augment"
#' @param me_index vector of index of Main Effects to be included in the design for model_order "non_standard"
#' @param qe_index vector index of Quadratic Effects to be included in the design for model_order "non_standard"
#' @param two_fi_index list of index of Two Factor Interactions to be included in the design for model_order "non_standard"
#' @param cubic_index vector of index of Cubic Effects to be included in the design for model_order "non_standard"
#' @param quatric_index vector of index of Quatric Effects to be included in the design for model_order "non_standard"
#' @param design a string indicating the type of design constructed. "Non-Bayesian" or "Bayesian" optimal design. The major difference is in how singularity of the matrix is evaluated. In a Bayesian optimal design X'X+K(prior variance matrix) is checked for singularity
#' @param K a diagonal matrix specifying the prior variance for Bayesian D-optimal design. By default set to 0.001 for potential terms and 0 for primary terms.
#' @param n_starts number of random starts, default at 100
#' @description function for initiating the coordinate exchange for a D-optimal design for a given number of random starts and returns the best design. This is a serial version of the code. Serial version of the code and will work for Windows as well.
#' @return a list of
#'  \itemize{
#'        \item Bayes_dopt_design- the best Bayesian D-optimal design matrix of the specified order including the intercept from the random starts
#'        \item Bayes_dopt_x- a matrix with factor levels of the best Bayesian D-optimal design from the random starts
#'        \item Bayes_dopt- The value of best objective value of the Bayesian D-optimal design from the random starts
#'        }
#'  @examples
#'  \dontrun{
#'   construct_des_D_Bayes_serial(x_matrix_start,x_start,n,C_dtype, C,model_order, type= "augment", freeze_rows=0,me_index = c(),qe_index = c(),two_fi_index = c(),cubic_index = c(),quatric_index = c(),n_starts=500,K,design = "Bayesian")
#'   }
#' @export
construct_des_D_Bayes_serial=function(x_matrix_start,x_start,n,C_dtype, C,
                               model_order, type= "augment", freeze_rows=0,
                               me_index = c(),qe_index = c(),
                               two_fi_index = c(),
                               cubic_index = c(),
                               quatric_index = c(),
                               n_starts=100,
                               K,design = "Bayesian"){
  best_D=-1e50
  for(i in 1:n_starts){
    #sample the starting design for each random start
    print(paste("random start",i))
    if(type=="augment"){
      freeze_rows=nrow(x_matrix_start)
      singular=TRUE
      while(singular){
        print("generating random x matrix")
        x=gen_x(n=n-freeze_rows,C,model_order=model_order,C_dtype,
                type="augment",
                me_index = me_index,qe_index = qe_index,
                two_fi_index = two_fi_index,
                cubic_index = cubic_index,
                quatric_index = quatric_index,
                design = "Bayesian")
        x_matrix=x[[2]]
        #print(x_matrix)
        #print(dim(x_matrix_start))
        x_matrix=rbind(x_matrix_start,x_matrix)
        x=x[[1]]
        if(class(x)!="matrix"){x=matrix(x,nrow=1)}
        #print(x)
        #print(x_start)
        x=rbind(x_start,x)
        print(rcond(t(x_matrix)%*%x_matrix + K))
        if(rcond(t(x_matrix)%*%x_matrix + K)>1e-15){
          singular=FALSE
        }
      }
    }else if(type=="non_augment"){
      x=gen_x(n=n,C,model_order=model_order,C_dtype,
              type="non_augment",
              me_index = me_index,qe_index = qe_index,
              two_fi_index = two_fi_index,
              cubic_index = cubic_index,
              quatric_index = quatric_index,
              design = "Bayesian")
      x_matrix=x[[2]]
      x=x[[1]]
    }
    if(class(x)!="matrix"){
      x=matrix(x)
    }
    print(paste("Coordinate exchange for try",i))
    #for each random start call the coordinate exchange function
    des=coor_exch_D_Bayes(x=x,x_matrix=x_matrix,C=C,model_order=model_order,C_dtype=C_dtype,
                          freeze_rows = freeze_rows,
                          me_index=me_index,qe_index=qe_index,
                          two_fi_index=two_fi_index,
                          cubic_index = cubic_index,
                          quatric_index = quatric_index,
                          type=type,K=K)
    current=des[[2]]
    x_current=des[[3]]
    des=des[[1]]
    print(paste("Current valus is",current))
    #stores the best design from the random starts
    if(current>best_D){
      print("best_des is getting stored")
      best_D=current
      best_des_D=des
      x_best_D=x_current
      print(paste("Best val is", best_D))
      #break
    }
    print("------------------------------------")
  }
  #if(!exists(best_des_D)){diag(K)}
  return(list("Bayes_dopt_design"=best_des_D,"Bayes_dopt_x"=x_best_D,"Bayes_dopt"=best_D))
}

#' function for initiating the individual coordinate exchange random starts
#' @param i index of random starts
#' @param x_matrix_start the design matrix of the specified order including the intercept of the random starting
#' @param x_start the settings matrix of the random starting design
#' @param n total number of runs in the design
#' @param C list of candidate space for each experimental factor in the design
#' @param model_order a string indicating the order of the regression model. "linear","quadratic" (Main Effects+Quadratic Effects), "2FI" (Main Effects + TwoFactor Interactions),"full", "non_standard"- for a user defined model order order
#' @param C_dtype a string vector indicating the data type of the experimental factors in the design either "cont"(continuous) or "cat"(categorical)
#' @param type a string indicating the type of design augmentation. "non_augment" or "augment" (augment an existing design)
#' @param freeze_rows the rows of the matrix to freeze while augmenting, 0 when type is "non_augment"
#' @param me_index vector of index of Main Effects to be included in the design for model_order "non_standard"
#' @param qe_index vector index of Quadratic Effects to be included in the design for model_order "non_standard"
#' @param two_fi_index list of index of Two Factor Interactions to be included in the design for model_order "non_standard"
#' @param cubic_index vector of index of Cubic Effects to be included in the design for model_order "non_standard"
#' @param quatric_index vector of index of Quatric Effects to be included in the design for model_order "non_standard"
#' @param design a string indicating the type of design constructed. "Non-Bayesian" or "Bayesian" optimal design. The major difference is in how signularity of the matrix is evaluated. In a Bayesian optimal design X'X+K(prior variance matrix) is checked for singularity
#' @param K a diagonal matrix specifying the prior variance for Bayesian D-optimal design. By default set to 0.001 for potential terms and 0 for primary terms.
#' @param use_cpp a boolean value indicating whether to use the CPP version of coordinate exchange
#' @description function for initiating the individual coordinate exchange random starts.
#' @return a list of
#'  \itemize{
#'        \item x_matrix- a Bayesian D-optimal design matrix of the specified order
#'        \item optimal_value- Optimal value of the objective value of the Bayesian D-optimal design
#'        \item x- a matrix with factor levels of the Bayesian D-optimal design
#'        }
#'  @examples
#'  \dontrun{
#'    individual_starts(i,x_matrix_start,x_start,n,C_dtype, C,model_order, type= "augment", freeze_rows=0,me_index = c(),qe_index = c(),two_fi_index = c(),cubic_index = c(),quatric_index = c(),n_starts=500,K,design = "Bayesian")
#'    }
#' @export
individual_starts=function(i,x_matrix_start,x_start,n,C_dtype, C,
                           model_order, type, freeze_rows,
                           me_index ,qe_index ,
                           two_fi_index ,
                           cubic_index ,
                           quatric_index,
                           K,design,use_cpp=TRUE){
  print(paste("random start",i))
  if(type=="augment"){
    freeze_rows=nrow(x_matrix_start)
    singular=TRUE
    while(singular){
      print("generating random x matrix")
      x=gen_x(n=n-freeze_rows,C,model_order=model_order,C_dtype,
              type="augment",
              me_index = me_index,qe_index = qe_index,
              two_fi_index = two_fi_index,
              cubic_index = cubic_index,
              quatric_index = quatric_index,
              design = "Bayesian")
      x_matrix=x[[2]]
      #print(dim(x_matrix))
      #print(dim(x_matrix_start))
      x_matrix=rbind(x_matrix_start,x_matrix)
      x=x[[1]]
      if(class(x)!="matrix"){x=matrix(x,nrow=1)}
      #print(x)
      #print(x_start)
      x=rbind(x_start,x)
      print(rcond(t(x_matrix)%*%x_matrix + K))
      if(rcond(t(x_matrix)%*%x_matrix + K)>1e-15){
        singular=FALSE
      }
    }
  }else if(type=="non_augment"){
    x=gen_x(n=n,C,model_order=model_order,C_dtype,
            type="non_augment",
            me_index = me_index,qe_index = qe_index,
            two_fi_index = two_fi_index,
            cubic_index = cubic_index,
            quatric_index = quatric_index,
            design = "Bayesian")
    x_matrix=x[[2]]
    x=x[[1]]
  }
  if(class(x)!="matrix"){
    x=matrix(x)
  }
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
  #set.seed(i)
  des=coor_exch_D_Bayes(x=x,x_matrix=x_matrix,C=C,model_order=model_order,C_dtype=C_dtype,
                        freeze_rows = freeze_rows,
                        me_index=me_index_,qe_index=qe_index_,
                        two_fi_index=two_fi_index_,
                        #cubic_index = cubic_index,
                        #quatric_index = quatric_index,
                        type=type,K=K,telescoping = FALSE)
  print("Starting Telescoping Stage")
  des=coor_exch_D_Bayes(x=des[[3]],x_matrix=des[[1]],C=C,model_order=model_order,C_dtype=C_dtype,
                        freeze_rows = freeze_rows,
                        me_index=me_index_,qe_index=qe_index_,
                        two_fi_index=two_fi_index_,
                        #cubic_index = cubic_index,
                        #quatric_index = quatric_index,
                        type=type,K=K,telescoping = TRUE)#,
                        #pri_index=pri_index,alpha=alpha,maxmin = maxmin)
  return(des)
}


#' function for initiating the coordinate exchange for a given number of random starts
#' @param x_matrix_start the design matrix of the specified order including the intercept of the random starting
#' @param x_start the settings matrix of the random starting design
#' @param n total number of runs in the design
#' @param C list of candidate space for each experimental factor in the design
#' @param model_order a string indicating the order of the regression model. "linear","quadratic" (Main Effects+Quadratic Effects), "2FI" (Main Effects + TwoFactor Interactions),"full", "non_standard"- for a user defined model order order
#' @param C_dtype a string vector indicating the data type of the experimental factors in the design either "cont"(continuous) or "cat"(categorical)
#' @param type a string indicating the type of design augmentation. "non_augment" or "augment" (augment an existing design)
#' @param freeze_rows the rows of the matrix to freeze while augmenting, 0 when type is "non_augment"
#' @param me_index vector of index of Main Effects to be included in the design for model_order "non_standard"
#' @param qe_index vector index of Quadratic Effects to be included in the design for model_order "non_standard"
#' @param two_fi_index list of index of Two Factor Interactions to be included in the design for model_order "non_standard"
#' @param cubic_index vector of index of Cubic Effects to be included in the design for model_order "non_standard"
#' @param quatric_index vector of index of Quatric Effects to be included in the design for model_order "non_standard"
#' @param design a string indicating the type of design constructed. "Non-Bayesian" or "Bayesian" optimal design. The major difference is in how singularity of the matrix is evaluated. In a Bayesian optimal design X'X+K(prior variance matrix) is checked for singularity
#' @param K a diagonal matrix specifying the prior variance for Bayesian D-optimal design. By default set to 0.001 for potential terms and 0 for primary terms.
#' @param n_starts number of random starts, default at 100
#' @param rng seed for random number generation for the random design construction starts
#' @param use_cpp a boolean value indicating whether to use the CPP version of coordinate exchange
#' @param num_cores number of cores for parallel processing. Default value is NA, in which case Sys.getenv("PBS_RESC_TOTAL_PROCS") is used to get the number of cores in a HPC. Otherwise, the user needs to provide a value. Use 1 for serial execution.
#' @description function for initiating the coordinate exchange for a given number of random starts and returns the best design. This is a parallel version of the code and works only on linux machines. For windows the makeCluster needs to be modified. By default the parallel version of the code will be run.
#' @return a list of
#'  \itemize{
#'        \item Bayes_dopt_design- the best Bayesian D-optimal design matrix of the specified order including the intercept from the random starts
#'        \item Bayes_dopt_x- a matrix with factor levels of the best Bayesian D-optimal design from the random starts
#'        \item Bayes_dopt- The value of best objective value of the Bayesian D-optimal design from the random starts
#'        }
#'  @examples
#'   construct_des_D_Bayes(x_matrix_start,x_start,n,C_dtype, C,model_order, type= "augment", freeze_rows=0,me_index = c(),qe_index = c(),two_fi_index = c(),cubic_index = c(),quatric_index = c(),n_starts=500,K,design = "Bayesian")
#' @import doRNG
#' @importFrom foreach %dopar%
#' @export
construct_des_D_Bayes=function(x_matrix_start,x_start,n,C_dtype, C,
                               model_order, type= "augment", freeze_rows=0,
                               me_index = c(),qe_index = c(),
                               two_fi_index = c(),
                               cubic_index = c(),
                               quatric_index = c(),
                               n_starts=100,
                               K,design = "Bayesian",rng=c(),use_cpp=TRUE, num_cores=NA){#,
                               #pri_index=c(),alpha=c(),maxmin=c()){
  best_D=-1e50
  #cl=makeCluster(detectCores()-1, type="FORK")
  #spawn the cluster for parallel computing, works only on linux
  if(is.na(num_cores)){
    cl=parallel::makeCluster(as.numeric(Sys.getenv("PBS_RESC_TOTAL_PROCS")) - 1, type="FORK",outfile="")
    print(paste("Number of workers is:",as.numeric(Sys.getenv("PBS_RESC_TOTAL_PROCS")) - 1))
  }else{
    cl=parallel::makeCluster(num_cores, type="FORK",outfile="")
    print(paste("Number of workers is:",num_cores))
  }
  doParallel::registerDoParallel(cl)
  #registerDoRNG(seed = 1234)
  #registerDoRNG(1234)
  #registerDoParallel(cores=as.numeric(Sys.getenv("PBS_RESC_TOTAL_PROCS")) - 1)
  mcoptions = list(preschedule=FALSE)
  #chunks = getDoParWorkers()
  #parallel processing for the multistarts
  multistart_results=foreach::foreach(i=1:n_starts, .options.multicore=mcoptions,.options.RNG=rng) %dorng%   individual_starts(i,x_matrix_start,x_start,n,C_dtype, C,
                                                                                                                     model_order, type= "augment", freeze_rows,
                                                                                                                     me_index ,qe_index ,
                                                                                                                     two_fi_index ,
                                                                                                                     cubic_index ,
                                                                                                                     quatric_index,
                                                                                                                     K=K,design = "Bayesian",use_cpp=TRUE)#,
                                                                                                                     #pri_index,alpha,maxmin)
  #stopImplicitCluster()
  parallel::stopCluster(cl)
  print("Done with multistarts")
  print(paste("Length of output is:",length(multistart_results)))

  #stores the best design
  for(i in 1:length(multistart_results)){
    current=multistart_results[[i]][[2]]
    x_current=multistart_results[[i]][[3]]
    des=multistart_results[[i]][[1]]
    if(current>best_D){
      best_D=current
      best_des_D=des
      x_best_D=x_current
    }
  }
  number_best_converegence=0
  criteria_list=c()
  for(i in 1:length(multistart_results)){
    criteria_list=c(criteria_list,multistart_results[[i]][[2]])
    if(isTRUE(all.equal(best_D,multistart_results[[i]][[2]]))){
      if(all(x_best_D==multistart_results[[i]][[3]])){
        number_best_converegence=number_best_converegence+1
      }
    }
  }
  #print(paste("Number of iterations converged to the best value is",number_best_converegence))
  graphics::plot(x=1:length(multistart_results),y=criteria_list, main="Plot of D-Optimal criteria", xlab="Iteration",ylab="Criteria_val")
  graphics::abline(h=best_D, col="red")
  #print("Top 50 results")
  #print(criteria_list[order(criteria_list,decreasing=TRUE)[1:50]])

  #   current=des[[2]]
  #   x_current=des[[3]]
  #   des=des[[1]]
  #   print(paste("Current valus is",current))
  #   if(current>best_D){
  #     print("best_des is getting stored")
  #     best_D=current
  #     best_des_D=des
  #     x_best_D=x_current
  #     print(paste("Best val is", best_D))
  #     #break
  #   }
  #   print("------------------------------------")
  # }
  # #if(!exists(best_des_D)){diag(K)}

  return(list("Bayes_dopt_design"=best_des_D,"Bayes_dopt_x"=x_best_D,"Bayes_dopt"=best_D))
}
