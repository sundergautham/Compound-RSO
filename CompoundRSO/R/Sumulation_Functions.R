#' selects the active first-order and second-order terms for the stage 1 design using JN model selection method
#' @description selects the active first-order and second-order terms for the stage 1 design using JN model selection method
#' @param design_matrix DSD settings matrix with fake factors
#' @param fake_factor_index index of fake factor columns
#' @param C_dtype a string vector indicating the data type of the experimental factors in the design either "cont"(continuous) or "cat"(categorical)
#' @param heredity string indicating the heredity of the model "strong", "weak" or "no"
#' @param model_selection_type string indicating the model selection type "JN_org", "LASSO", "APR" (All Possible Regression) or "StepAICc". We recommend using either JN_org or LASSO as the criterion for model selection.
#' @param rm_valruns index of validation runs to remove while fitting the model, useful when using cross validation approaches
#' @param y the response values for the experiment
#' @return a list of
#'  \itemize{
#'        \item ME- active main effects from DSD model selection
#'        \item 2FI- active 2FI from the DSD model selection
#'        \item QE- active QE from the DSD model selection
#'        }
#' @import doRNG
#' @export
dsd_model_selection=function(design_matrix,fake_factor_index,
                             C_dtype,heredity="strong",model_selection_type="JN_org",rm_valruns=c(),y){
  if(length(rm_valruns)>0){
    print("Removing valruns for screening")
    design_matrix_=design_matrix[-rm_valruns,]
    y=y[-rm_valruns]
    print(dim(design_matrix_))
    print(length(y))
  }else{
    design_matrix_=design_matrix
  }
  JN_result=JN_model_selection(dsd_matrix = design_matrix_, y=y, alpha=0.2,
                               heredity = heredity,fake_factor_index = fake_factor_index,
                               C_dtype = C_dtype,model_selection_type=model_selection_type)

  active_me=colnames(design_matrix)[JN_result[[1]]]
  active_se=colnames(JN_result[[5]])
  qe_index=sapply(active_se,FUN = function(x){stringr::str_split(x,"[*]",simplify = TRUE)[1]==stringr::str_split(x,"[*]",simplify = TRUE)[2]})
  active_qe=active_se[qe_index]
  active_2FI=active_se[!qe_index]
  #print(active_qe)
  if(length(active_qe)==0){
    active_qe=c()
  }else if(is.na(active_qe)|is.null(active_qe)){
    active_qe=c()
  }
  #print(active_2FI)
  if(length(active_2FI)==0){
    active_2FI=c()
  }else if(is.na(active_2FI)|is.null(active_2FI)){
    active_2FI=c()
  }
  return_list=list()
  return_list[["ME"]]=active_me
  return_list[["2FI"]]=active_2FI
  return_list[["QE"]]=active_qe
  return(return_list)
}

#' Augments the starting design DSD with compound uniform design for a specified number of runs
#' @description Augments the starting design DSD with compound uniform design for a specified number of runs
#' @param active_me vector of active ME from the stage 1 regression
#' @param active_2FI vector of active 2FI from the stage 1 regression
#' @param active_qe vector of active QE from the stage 1 regression
#' @param design_matrix DSD matrix with fake factors
#' @param naug number of runs to augment the starting design by
#' @param fake_factor_index index of fake factor columns
#' @param C_dtype data type of the factors in the design either "cont"(continuous) or "cat"(categorical)
#' @param C list of candidate space for each experimental factor in the design
#' @param heredity string indicating the heredity of the model "strong", "weak" or "no". Determines whether inactive factors should be dropped.
#' @param se_names column names of the design matrix
#' @param unif_crit Uniform discrepancy criterion either "MD2" (default) and "WD2"
#' @param design_aug_type if "compound_unif" computes the compound uniform criterion
#' @param n_starts number of random starts when constructing D-optimal designs
#' @param n_starts_compound number of random starts when constructing Compound Uniform (robust) designs
#' @param U_opt_tresh a numeric value providing the lower bound of the threshold in bisection search for the Uniform optimal criteria. Default value is 0.8
#' @param D_opt_tresh a numeric value providing the lower bound of the threshold in bisection search for the D-optimal criteria. Default value is 0.8
#' @param w_tresh	 the threshold for difference in weights w for the bisection search to converge. The default value is 0.02
#' @param use_w a boolean value indicating whether w_candidate should be used or the bisection search method for finding the compound optimal designs
#' @param w_candidate if use_w is TRUE, a vector of weights w to compute the Uniform optimal design for
#' @param num_cores number of cores for parallel processing. Default value is NA, in which case Sys.getenv("PBS_RESC_TOTAL_PROCS") is used to get the number of cores in a HPC. Otherwise, the user needs to provide a value. Use 1 for serial execution.
#' @param rng random seed for random number generation
#' @return A list of:
#'     \itemize{
#'       \item duag_design-D-optimal desugn for the augmented runs
#'       \item duag_x- factor levels matrix of the D-optimal design for the augmented runs
#'       \item compound_unif_design- A list of compound uniform optimal designs for each weigth w
#'       \item design_eff_results- a matrix siummarizing the design effciencies of the compound optimal design
#'     }
#' @import doRNG
#' @importFrom foreach %dopar%
#' @export
stage2_aug=function(active_me,active_2FI,active_qe,
                    design_matrix,fake_factor_index,naug,se_names,
                    C_dtype,C,heredity="no",n_starts=100,n_starts_compound=10,design_aug_type="compound_unif",
                    unif_crit="MD2", U_opt_tresh=0.8, D_opt_tresh=0.8,w_tresh = 0.02,
                    use_w=FALSE, w_candidate=0, rng=c(),num_cores=NA){
  if(length(rng)==0){rng=sample.int(1e5,size=1)}
  print(paste("Random number seed",rng))
  nfactor=length(C)
  n=nrow(design_matrix)
  me_index=c(1:nfactor)
  qe_index=which(C_dtype=="cont")
  two_fi_index=utils::combn(1:nfactor,2,simplify = FALSE)
  #me_index=which(names(C) %in% active_me_list[[sim_i]])
  if(heredity=="weak"){
    del=two_fi_index[sapply(two_fi_index,FUN=function(x){!(any(x %in% which(se_names[-1]%in%active_me)))})]
    del=sapply(del,function(x){paste0("X",x[1],"*","X",x[2])})
    print("Terms to be deleted under weak heredity")
    print(del)
    #print("Active ME")
    #print(active_me_list[[sim_i]])
    two_fi_index=two_fi_index[sapply(two_fi_index,FUN=function(x){any(x %in% which(se_names[-1]%in%active_me))})]
    print(paste("Length of weak interaction terms",length(two_fi_index)))
    #update se_names with the terms deleted
    if(length(del)>0){se_names=se_names[-which(se_names %in% del)]}
    print("Total number of terms in weak heredity model")
    print(se_names)
  }else if(heredity=="strong"){
    del=two_fi_index[sapply(two_fi_index,FUN=function(x){!(all(x %in% which(se_names[-1]%in%active_me)))})]
    del=sapply(del,function(x){paste0("X",x[1],"*","X",x[2])})
    print("Terms to be deleted under weak heredity")
    print(del)
    #print("Active ME")
    #print(active_me_list[[sim_i]])
    two_fi_index=two_fi_index[sapply(two_fi_index,FUN=function(x){all(x %in% which(se_names[-1]%in%active_me))})]
    #print(paste("Length of weak interaction terms",length(two_fi_index)))
    if(length(del)>0){se_names=se_names[-which(se_names %in% del)]}
      #print("Total number of terms in weak model")
      #print(se_names)
  }else if(heredity=="no"){
    two_fi_index
  }else{
    stop("Error in heredity value")
  }
  K=rep(0.001,length(se_names))
  #intercept
  K[1]=0
  #active me
  K[which(se_names%in%active_me)]=0
  #active 2FI
  K[which(se_names%in%active_2FI)]=0
  #all quadratic terms including ones with no active ME
  K[which(se_names%in%paste0(names(C)[qe_index],"*",names(C)[qe_index]))]=0
  #K[which(se_names%in%active_cubic)]=0
  #K[which(se_names%in%active_quatric)]=0
  #K=K**2
  print(K)
  K=diag(K)
  x_matrix_start=t(apply(design_matrix[,1:nfactor],MARGIN=1,FUN=f_x,model_order="non_standard",
                           me_index=me_index,two_fi_index=two_fi_index,qe_index=qe_index,
                           C_dtype=C_dtype))
  dsd_bayes_D_aug=construct_des_D_Bayes(x_matrix_start=x_matrix_start,
                                          x_start=design_matrix[,1:nfactor],
                                          n=nrow(design_matrix)+naug,C_dtype=C_dtype, C=C,
                                          model_order="non_standard", type= "augment",
                                          freeze_rows=nrow(design_matrix),
                                          me_index = me_index,qe_index = qe_index,
                                          two_fi_index = two_fi_index,
                                          cubic_index=c(),quatric_index=c(),
                                          n_starts=n_starts,
                                          K=K,design = "Bayesian",
                                          rng=rng,num_cores=num_cores)
  daug_mat=dsd_bayes_D_aug[[1]]
  #print(daug_mat)
  colnames(daug_mat)=se_names
  #bayes_daug_list[[sim_i]]=daug_mat
  x_daug_mat=dsd_bayes_D_aug[[2]]
  colnames(x_daug_mat)=se_names[2:(nfactor+1)]
  #bayes_daug_x_mat_list[[sim_i]]=x_daug_mat
  if(design_aug_type=="compound_unif"){
    if(unif_crit=="WD2"){
      unif_start=matrix(stats::runif(nfactor*(nrow(design_matrix)+naug)),ncol=nfactor,nrow=nrow(design_matrix)+naug)
      unif_opt_design=DiceDesign::discrepESE_LHS(unif_start, T0=0.005*DiceDesign::discrepancyCriteria(unif_start,type='W2')[[1]],
                                       inner_it=100, J=50, it=2, criterion="W2")
      x_unifopt=unif_opt_design$design
    }else if (unif_crit=="MD2"){
      unif_opt_design=UniDOE::GenUD(n=nrow(design_matrix)+naug,s=nfactor,q=nrow(design_matrix)+naug)
      x_unifopt=unif_opt_design$final_design
      x_unifopt=(x_unifopt-0.5)/(nrow(design_matrix)+naug)
    }
    compound_unif_designs_cpp=construct_des_compound_unif(x_matrix_start=x_matrix_start,x_start=design_matrix[,1:nfactor],n=nrow(design_matrix)+naug,C_dtype=C_dtype,C=C,
                                                            model_order="non_standard", type= "augment", aug_type="not_foldover",
                                                            D_opt_design=daug_mat,D_opt_x=x_daug_mat,x_unifopt=x_unifopt,unif_crit=unif_crit,
                                                            n_starts=n_starts_compound,D_opt_type="Bayesian",
                                                            w_candidate=w_candidate,use_w=use_w,use_cpp=TRUE,
                                                            K=K,me_index_daug=me_index, qe_index_daug=qe_index,two_fi_index_daug=two_fi_index,
                                                            D_opt_tresh = D_opt_tresh, U_opt_tresh=U_opt_tresh,w_tresh = w_tresh, rng=rng, num_cores=num_cores)
    results_data=data.frame("w"=NA,"w_criteria"=NA,"D_eff"=NA,"U_eff"=NA)
    if(!use_w){
      results_data[1,1]=names(compound_unif_designs_cpp)[1]
      results_data[1,2]=compound_unif_designs_cpp[[1]][[3]][1]
      results_data[1,3]=compound_unif_designs_cpp[[1]][[3]][2]
      results_data[1,4]=compound_unif_designs_cpp[[1]][[3]][3]
      row=2
      for(i in 2:(length(compound_unif_designs_cpp)-1)){
        results_data[row,1]=names(compound_unif_designs_cpp)[i]
        results_data[row,2]=compound_unif_designs_cpp[[i]][[3]][1]
        results_data[row,3]=compound_unif_designs_cpp[[i]][[3]][2]
        results_data[row,4]=compound_unif_designs_cpp[[i]][[3]][3]
        row=row+1
      }
      results_data[row,1]=names(compound_unif_designs_cpp)[length(compound_unif_designs_cpp)]
      results_data[row,2]=compound_unif_designs_cpp[[length(compound_unif_designs_cpp)]][[2]][1]
      results_data[row,3]=compound_unif_designs_cpp[[length(compound_unif_designs_cpp)]][[2]][2]
      results_data[row,4]=compound_unif_designs_cpp[[length(compound_unif_designs_cpp)]][[2]][3]
      compound_unif_designs_cpp_list=list("duag_design"=daug_mat,"duag_x"=x_daug_mat,"compound_unif_design"=compound_unif_designs_cpp,"design_eff_results"=results_data)
    }else{
      row=1
      for(i in 1:length(compound_unif_designs_cpp)){
        results_data[row,1]=names(compound_unif_designs_cpp)[i]
        results_data[row,2]=compound_unif_designs_cpp[[i]][[3]][1]
        results_data[row,3]=compound_unif_designs_cpp[[i]][[3]][2]
        results_data[row,4]=compound_unif_designs_cpp[[i]][[3]][3]
        row=row+1
      }
      compound_unif_designs_cpp_list=list("duag_design"=daug_mat,"duag_x"=x_daug_mat,"compound_unif_design"=compound_unif_designs_cpp,"design_eff_results"=results_data)
    }
    print(results_data)
  }else{
    compound_unif_designs_cpp_list=list("duag_design"=daug_mat,"duag_x"=x_daug_mat,"compound_unif_design"=c(),"design_eff_results"=c())
  }
  return(compound_unif_designs_cpp_list)
}


#' fits a GAUSS-LASSO regression model (LASSO with AICc is used for model selection and an OLS model is fit for the identified terms), on the stage 2 design and returns the lack of fit test result
#' @param design_matrix design matrix from Stage 2 (the matrix with all the regression terms)
#' @param y a vector of responses for the experiment
#' @param se_names column names of the design matrix
#' @param lof_alpha significance level for the Lack of Fit test
#' @return a list of:
#'  \itemize{
#'     \item beta- active terms from model and the Gauss-LASSO estimates for the active terms
#'     \item Sigma- the variance-covariance matrix for the regression model
#'     \item F_stat_lof- F statistic of lack of fit test
#'     \item F_crit_lof- F critical value of lack of fit test
#'     \item lof_result- lack of fit test result
#'     \item resid- the residuals of the Gauss LASSO model
#'    }
#' @export

stage2_aug_analysis=function (design_matrix, y, se_names, lof_alpha = 0.05)
{
  if (length(y) != nrow(design_matrix)) {
    stop("Error: number of design responses and the number of experiments should be same")
  }
  fitlasso = gamlr::gamlr(x = design_matrix[, -1], y = y, family = "gaussian",
                          nlambda = 100)
  beta = as.matrix(stats::coef(fitlasso))
  colnames(design_matrix) = se_names
  X_ols = design_matrix[, which(beta != 0)]
  X_ols = cbind(X_ols, y)
  colnames(X_ols) = c(se_names[beta != 0], "y")
  lm_ols = stats::lm(y ~ 0 + ., data = as.data.frame(X_ols))
  beta_ols = lm_ols$coefficients
  names(beta_ols) = se_names[beta != 0]
  print("Done with fitting OLS models")
  X_ols = X_ols[, -ncol(X_ols),drop=FALSE]
  Sigma_ols = (sum(lm_ols$residuals^2)/lm_ols$df.residual) *
    solve(t(X_ols) %*% X_ols)
  active_terms = names(beta_ols)
  active_factors = unique(unlist(sapply(active_terms, FUN = function(x) {
    stringr::str_split(x, "[*]", simplify = TRUE)
  }, simplify = TRUE)))
  active_factors = se_names[which(se_names %in% active_factors)]
  print(paste("active_factors",active_factors))
  if(length(active_factors)==1&active_factors=="intercept"){
    print("####################################################")
    print("Intercept only model hence the design matrix is used to identify pure error terms")
    active_factors=se_names[sapply(stringr::str_split(se_names,"[*]"),FUN=function(x){length(x)==1})]
    print("####################################################")
  }
  X_ols_pe = design_matrix[, active_factors,drop=FALSE]
  unique_rows = unique(X_ols_pe)
  print(is.null(unique_rows))
  rep_index = t(apply(unique_rows, MARGIN = 1, FUN = function(x) {
    apply(X_ols_pe, MARGIN = 1, FUN = function(y) {
      all(y == x)
    })
  }))
  rep_index = apply(rep_index, MARGIN = 1, FUN = function(z) {
    which(z == TRUE)
  })
  rep_index = rep_index[sapply(rep_index, FUN = function(x) {
    length(x) > 1
  })]
  sse = sum(summary(lm_ols)$residuals^2)
  print(paste("rep index", rep_index))
  print(paste("ybar", sapply(rep_index, FUN = function(x) {
    mean(y[x])
  })))
  print(paste("y", sapply(rep_index, FUN = function(x) {
    y[x]
  })))
  sspe = sapply(rep_index, FUN = function(x) {
    sum((mean(y[x]) - y[x])^2)
  })
  sspe = sum(sspe)
  num = (sse - sspe)/(nrow(unique_rows) - ncol(X_ols))
  den = sspe/(nrow(design_matrix) - nrow(unique_rows))
  F_stat = num/den
  print(paste("sspe", sspe))
  print(paste("MSLF", num))
  print(paste("MSPE", den))
  num_df = (nrow(unique_rows) - ncol(X_ols))
  den_df = (nrow(design_matrix) - nrow(unique_rows))
  print(paste("Numerator df", num_df))
  print(paste("Denominator df", den_df))
  print(paste("F_stat", F_stat))
  F_crit = stats::qf(1 - lof_alpha, num_df, den_df)
  print(paste("F_crit", F_crit))
  test_result = (F_stat > F_crit)
  graphics::plot(lm_ols, which = 1, main = "Residual Plot")
  residuals = lm_ols$residuals
  me_names = se_names[-1][sapply(sapply(se_names[-1], FUN = function(x) {
    stringr::str_split(x, "[*]")
  }), function(x) {
    length(x) == 1
  })]
  qe_names = se_names[-1][sapply(sapply(se_names[-1], FUN = function(x) {
    stringr::str_split(x, "[*]")[1]
  }), function(x) {
    if (length(x) > 1) {
      x[1] == x[2]
    }
    else {
      FALSE
    }
  })]
  two_fi_names = se_names[-1][sapply(sapply(se_names[-1], FUN = function(x) {
    stringr::str_split(x, "[*]")[1]
  }), function(x) {
    if (length(x) > 1) {
      x[1] != x[2]
    }
    else {
      FALSE
    }
  })]
  active_me = me_names[me_names %in% colnames(X_ols)]
  active_qe = qe_names[qe_names %in% colnames(X_ols)]
  active_2FI = two_fi_names[two_fi_names %in% colnames(X_ols)]
  print(active_me)
  print(active_qe)
  print(active_2FI)
  print(colnames(X_ols))
  return_list = list()
  return_list[["beta"]] = beta_ols
  return_list[["Sigma"]] = Sigma_ols
  return_list[["F_stat_lof"]] = F_stat
  return_list[["F_crit_lof"]] = F_crit
  return_list[["lof_result"]] = test_result
  return_list[["resid"]] = residuals
  return(return_list)
}

#deprecated code for stage2_aug_analysis
# stage2_aug_analysis=function(design_matrix,y,se_names,lof_alpha=0.05){
#   #Sigma_ols_list=vector(mode="list",length=nsim)
#   #beta_ols_list=vector(mode="list",length=nsim)
#   #F_stat_list=vector(mode="numeric",length=nsim)
#   #F_val_list=vector(mode="numeric",length=nsim)
#   #test_result_list=vector(mode="logical",length=nsim)
#   #residuals_list=vector(mode="list",length=nsim)
#   #overall_unique_dsd=vector(mode="list",length=nsim)
#   #power_me_dsd=vector(mode="numeric",length=nsim)
#   #power_qe_dsd=vector(mode="numeric",length=nsim)
#   #power_2FI_dsd=vector(mode="numeric",length=nsim)
#
#   #fdr_me_dsd=vector(mode="numeric",length=nsim)
#   #fdr_qe_dsd=vector(mode="numeric",length=nsim)
#   #fdr_2FI_dsd=vector(mode="numeric",length=nsim)
#
#   #se_names=se_names_list[[sim_i]]
#   #design_matrix=design_aug_list[[sim_i]]
#   #y_aug=apply(design_x[[sim_i]][(freeze_rows+1):(freeze_rows+naug),,drop=FALSE],MARGIN=1,y_func)
#   #y=c(y_[[sim_i]],y_aug)
#   #y_[[sim_i]]=y
#   #free_index=which(se_names %in% active_me_list[[sim_i]])-1
#   #fitlasso_me_nopenalty=gamlr(x=design_matrix[,-1],
#   #                        y=y_[[sim_i]],family = "gaussian",
#   #                        nlambda = 100,free=free_index)
#   if(length(y)!=nrow(design_matrix)){
#     stop("Error: number of design responses and the number of experiments should be same")
#   }
#   fitlasso=gamlr::gamlr(x=design_matrix[,-1],
#                    y=y,family = "gaussian",
#                    nlambda = 100)
#   beta=as.matrix(stats::coef(fitlasso))
#   #beta_me_nopenalty=as.matrix(coef(fitlasso_me_nopenalty,corrected = TRUE))
#   colnames(design_matrix)=se_names
#   #print(which(beta!=0))
#   X_ols=design_matrix[,which(beta!=0)]
#   #X_ols_me_nopenalty=design_matrix[,which(beta_me_nopenalty!=0)]
#   #print(dim(X_ols))
#   #print(length(c(se_names[beta!=0],"y")))
#   X_ols=cbind(X_ols,y)
#   #X_ols_me_nopenalty=cbind(X_ols_me_nopenalty,y_[[sim_i]])
#   colnames(X_ols)=c(se_names[beta!=0],"y")
#   #colnames(X_ols_me_nopenalty)=c(se_names[beta_me_nopenalty!=0],"y")
#   #print("Starting ols models")
#   lm_ols=stats::lm(y~0+.,data=as.data.frame(X_ols))
#   #lm_ols_me_nopenalty=lm(y~0+.,data=as.data.frame(X_ols_me_nopenalty))
#   beta_ols=lm_ols$coefficients
#   #beta_ols_me_nopenalty=lm_ols_me_nopenalty$coefficients
#   names(beta_ols)=se_names[beta!=0]
#   #names(beta_ols_me_nopenalty)=se_names[beta_me_nopenalty!=0]
#   print("Done with fitting OLS models")
#   #beta_ols_list[[sim_i]]=beta_ols
#   #beta_ols_freeindex_list[[sim_i]]=beta_ols_me_nopenalty
#   X_ols=X_ols[,-ncol(X_ols)]
#   #print(paste("Active terms in stage 2 analysis",colnames(X_ols)))
#   #X_ols_me_nopenalty=X_ols_me_nopenalty[,-ncol(X_ols_me_nopenalty)]
#   Sigma_ols=(sum(lm_ols$residuals**2)/lm_ols$df.residual)*solve(t(X_ols)%*%X_ols)
#   #Sigma_ols_freeindex_list[[sim_i]]=(sum(lm_ols_me_nopenalty$residuals**2)/lm_ols_me_nopenalty$df.residual)*solve(t(X_ols_me_nopenalty)%*%X_ols_me_nopenalty)
#   active_terms=names(beta_ols)
#   active_factors=unique(unlist(sapply(active_terms,FUN=function(x){stringr::str_split(x,"[*]",simplify = TRUE)},simplify = TRUE)))
#   #print(active_factors)
#   #print(colnames(X_ols))
#   active_factors=se_names[which(se_names%in%active_factors)]
#   X_ols_pe=design_matrix[,active_factors]
#   #X_ols_pe=X_ols_dsd[,active_me_dsd]
#   unique_rows=unique(X_ols_pe)
#   print(is.null(unique_rows))
#   rep_index=t(apply(unique_rows,MARGIN=1,FUN = function(x){apply(X_ols_pe,MARGIN=1,FUN=function(y){all(y==x)})}))
#   rep_index=apply(rep_index,MARGIN=1,FUN=function(z){which(z==TRUE)})
#   rep_index=rep_index[sapply(rep_index, FUN=function(x){length(x)>1})]
#   #X_ols_lm=as.data.frame(X_ols)
#   #y=return_list_dsd[[4]][[i]]
#   #X_ols_lm["y"]=y
#   #lm_ols=lm(y~.,data=X_ols_lm)
#   sse=sum(summary(lm_ols)$residuals**2)
#   print(paste("rep index",rep_index))
#   print(paste("ybar",sapply(rep_index,FUN=function(x){mean(y[x])})))
#   print(paste("y",sapply(rep_index,FUN=function(x){y[x]})))
#   sspe=sapply(rep_index,FUN=function(x){sum((mean(y[x])-y[x])**2)})
#   sspe=sum(sspe)
#   num=(sse-sspe)/(nrow(unique_rows)-ncol(X_ols))
#   den=sspe/(nrow(design_matrix)-nrow(unique_rows))
#   F_stat=num/den
#   print(paste("sspe",sspe))
#   print(paste("MSLF",num))
#   print(paste("MSPE",den))
#   num_df=(nrow(unique_rows)-ncol(X_ols))
#   den_df=(nrow(design_matrix)-nrow(unique_rows))
#   print(paste("Numerator df", num_df))
#   print(paste("Denominator df", den_df))
#   print(paste("F_stat",F_stat))
#   #F_stat_list[[sim_i]]=F_stat
#   F_crit=stats::qf(1-lof_alpha,num_df,den_df)
#   print(paste("F_crit",F_crit))
#   #F_val_list[[sim_i]]=F_crit
#   test_result=(F_stat>F_crit)
#   #y_hat=predict.lm(lm_ols,as.data.frame(X_ols[,-ncol(X_ols)]))
#   #graphics::plot(x=y_hat,y=lm_ols$residuals,main="Residual Plot",xlab="y predicted",ylab="residuals")
#   graphics::plot(lm_ols,which=1,main="Residual Plot")#,xlab="y predicted",ylab="residuals")
#   residuals=lm_ols$residuals
#   #overall_unique_dsd[[sim_i]]=active_factors
#   me_names=se_names[-1][sapply(sapply(se_names[-1],FUN=function(x){stringr::str_split(x,"[*]")}),function(x){length(x)==1})]
#   qe_names=se_names[-1][sapply(sapply(se_names[-1],FUN=function(x){stringr::str_split(x,"[*]")[1]}),function(x){if(length(x)>1){x[1]==x[2]}else{FALSE}})]
#   two_fi_names=se_names[-1][sapply(sapply(se_names[-1],FUN=function(x){stringr::str_split(x,"[*]")[1]}),function(x){if(length(x)>1){x[1]!=x[2]}else{FALSE}})]
#   active_me=me_names[me_names%in%colnames(X_ols)]
#   active_qe=qe_names[qe_names%in%colnames(X_ols)]
#   active_2FI=two_fi_names[two_fi_names%in%colnames(X_ols)]
#   print(active_me)
#   print(active_qe)
#   print(active_2FI)
#   print(colnames(X_ols))
#   #return the two beta lists and design_list
#   return_list=list()
#   return_list[["beta"]]=beta_ols
#   return_list[["Sigma"]]=Sigma_ols
#   return_list[["F_stat_lof"]]=F_stat
#   return_list[["F_crit_lof"]]=F_crit
#   return_list[["lof_result"]]=test_result
#   return_list[["resid"]]=residuals
#   return(return_list)
# }


#' constructs the Adaptive-RSO designs in Stage 3 design
#' @param beta a vector of regression parameter estimates from stage 2 analysis
#' @param design_matrix design matrix from Stage 2
#' @param design_x settings matrix from Stage 2
#' @param Sigma variance covariance matrix of the regression model from stage 2
#' @param se_names design matrix column names from stage 2
#' @param naug number of runs to augment design by in stage 3
#' @param n_starts number of random starts for constructing D-optimal and I-optimal designs
#' @param n_starts_compound number of random starts for design constructing compound optimal designs
#' @param B number of bootstrap samples for x* confidence region
#' @param C_dtype data type of the factors in the design either "cont"(continuous) or "cat"(categorical)
#' @param C  list of of candidate points for each factor in the experiment
#' @param me_index_daug vector of main effects index in the D-optimal design
#' @param qe_index_daug vector of quadratic index in the D-optimal design
#' @param two_fi_index_daug list of two factor interaction index in the D-optimal design
#' @param use_cpp boolean indicating whether to use CPP code or not
#' @param D_opt_tresh a numeric value providing the lower bound of the threshold in bisection search for the D-optimal criteria. Default value is 0.8
#' @param I_opt_tresh a numeric value providing the lower bound of the threshold in bisection search for the I-optimal criteria. Default value is 0.8
#' @param w_tresh	 the threshold for difference in weights w for the bisection search to converge. The default value is 0.02
#' @param minimize a boolean value indicating whether to minimize or maximize. Default value is minimize
#' @param rng seed for random number generation for the random design construction starts
#' @param num_cores number of cores for parallel processing. Default value is NA, in which case Sys.getenv("PBS_RESC_TOTAL_PROCS") is used to get the number of cores in a HPC. Otherwise, the user needs to provide a value. Use 1 for serial execution.
#' @return returns a list of:
#'       \itemize{
#'          \item D_aug- list of Bayesian D-optimal augmentation design matrix and factor levels matrix of the stage 2 design for the specified number of augmentation runs
#'          \item Boostrap_results- matrix containing the subregion of the design space idetified for exploitation
#'          \item I_aug- list of I-optimal augmentation design matrix and factor levels matrix of the stage 2 design for the specified number of augmentation runs
#'          \item compound_weff- list of the best exploration exploitation designs for a set of weights w evaluated by bisection search
#'          \item K- prior variance matrix for Bayesian D-optimal design
#'          \item W_mat- moment matrix for I-lambda designs
#'          \item x_star_beta- optimal location x* for the current estimate beta
#'          }
#' @export
stage3_aug=function(beta,design_matrix,design_x,Sigma,
                             se_names,naug,n_starts=100,n_starts_compound=100,B=5000,
                             C_dtype,C,
                             me_index_daug,qe_index_daug,two_fi_index_daug, use_cpp=TRUE,
                             D_opt_tresh = 0.8, I_opt_tresh = 0.8,w_tresh = 0.02,minimize=TRUE,rng = c(),num_cores=NA){
  if(length(rng)==0){rng=sample.int(1e5,size=1)}
  print(paste("Random number seed",rng))
  print(paste("use_cpp",use_cpp))
  nfactor=length(C)
  #se_names=se_names_list[[sim_i]]
  print(se_names)
  #active=se_names %in% names(beta_list[[sim_i]])
  #active=se_names[active]
  active=names(beta)
  active=active[-1]
  print(active)
  if(length(active)>0){
    me_index=sapply(active,FUN=function(x){length(stringr::str_split(x,"[*]",simplify = TRUE))==1},simplify = TRUE)
    active_me=active[me_index]
    active=active[!me_index]
    if(length(active)>0){
      qe_index=sapply(active,FUN=function(x){stringr::str_split(x,"[*]",simplify = TRUE)[1]==stringr::str_split(x,"[*]",simplify = TRUE)[2]},simplify = TRUE)
      active_qe=active[qe_index]
      active_2FI=active[!qe_index]
    }else{
      active_qe=c()
      active_2FI=c()
    }
  }else{
    active_me=c()
    active_qe=c()
    active_2FI=c()
  }
  #print(beta_list)
  #print(active_qe)
  #print(active_2FI)
  #print(which(se_names %in% active_me)-1)
  parmat=sapply(C,FUN=sample,size=20,replace=TRUE)
  #qe_index=as.numeric(regmatches(active_qe,regexpr("[0-9]+",active_qe)))
  #two_fi_index=lapply(regmatches(active_2FI,gregexpr("[0-9]+",active_2FI)),FUN=as.numeric)
  two_fi_index=lapply(active_2FI,FUN=function(x){which(names(C)%in%as.vector(stringr::str_split(x,"[*]",simplify=TRUE)))})
  qe_index=sapply(active_qe,FUN=function(x){which(names(C)%in%stringr::str_split(x,"[*]",simplify=TRUE)[1])})
  if(length(two_fi_index)==0){two_fi_index=c()}
  if(length(qe_index)==0){qe_index=c()}
  #print(two_fi_index)
  #print(qe_index)
  if(minimize){
    opt_multiplier=1
  }else{
    opt_multiplier=-1
  }
  opt_vals=optimr::multistart(parmat = parmat, fn=fx,method = "L-BFGS-B",
                        b=opt_multiplier*matrix(beta),model_order="non_standard",C_dtype=C_dtype,
                        me_index=which(se_names %in% active_me)-1,
                        qe_index=qe_index,
                        two_fi_index=two_fi_index,
                        lower = rep(-1,nfactor),upper=rep(1,nfactor))

  opt_vals=opt_vals[opt_vals[,"convergence"]==0,]

  x_star=as.matrix(opt_vals[which.min(opt_vals$value),1:nfactor])
  x_star_mat=x_star[1,]
  x_start=design_x
  #x_matrix_start=t(apply(x_start,MARGIN=1,FUN=f_x,me_index=1:nfactor,qe_index=qe_index_daug,two_fi_index=two_fi_index_daug,model_order="non_standard",C_dtype=C_dtype))
  x_matrix_start=design_matrix
  print(dim(x_matrix_start))
  K=rep(0.001,ncol(x_matrix_start))
  K[1]=0
  K[which(se_names%in%active_me)]=0
  K[which(se_names%in%active_2FI)]=0
  #change made on 8 July so that all the pure quadratic terms are active in stage 3
  qe_index_=which(C_dtype == "cont")
  K[which(se_names %in% paste0(names(C)[qe_index_], "*", names(C)[qe_index_]))] = 0
  #K[(length(K)-length(C)+1):length(K)]=0
  #K[which(se_names%in%paste0(active_me,"*",active_me))]=0
  #K=K**2
  #K_list[[sim_i]]=K
  K=diag(K)
  alpha=matrix(0,nrow=(nfactor+1),
               ncol=ncol(x_matrix_start)-(nfactor+1))
  maxmin=matrix(0,nrow=ncol(alpha),ncol=2)
  maxmin[,1]=1
  #C_dtype=rep("cont",ncol(x_start))
  #if(all(w_candidate!=0) | length(w_candidate)==0){
  start=proc.time()[3]
  dsd_bayes_D_aug=construct_des_D_Bayes(x_matrix_start=x_matrix_start,
                                          x_start=x_start,
                                          n=nrow(x_matrix_start)+naug,C_dtype=C_dtype, C=C,
                                          model_order="non_standard", type= "augment",
                                          freeze_rows=nrow(x_matrix_start),
                                          me_index = me_index_daug,qe_index =qe_index_daug,
                                          two_fi_index = two_fi_index_daug, n_starts=n_starts,
                                          K=K,design = "Bayesian",rng=rng,num_cores=num_cores)
  daug_mat=dsd_bayes_D_aug[[1]]
  end=proc.time()[3]
  print(paste("Time to execute Stage 3 D-aug",end-start))
  colnames(daug_mat)=se_names
  #bayes_d_compound[[sim_i]]=daug_mat
  x_daug_mat=dsd_bayes_D_aug[[2]]
  colnames(x_daug_mat)=colnames(x_start)
  #bayes_d_compound_x[[sim_i]]=x_daug_mat
  print("#############Done with Bayes D aug##################")

  start=proc.time()[3]
  boot_results=bootstrap_confidence_set(beta=beta,Sigma=Sigma,B=B,
                                          me_index=which(se_names %in% active_me)-1,
                                          qe_index=qe_index,two_fi_index=two_fi_index,
                                          C_dtype=C_dtype,nfactor=nfactor,
                                          C=C, model_order = "non_standard", lower=rep(-1,nfactor),upper=rep(1,nfactor),fx=fx,minimize=minimize,num_cores=num_cores)
  boot_results=boot_results[[1]]
  end=proc.time()[3]
  print(paste("Time to run bootstrap samples",end-start))
  boot_list=rbind(apply(boot_results,MARGIN=2,min),
                           apply(boot_results,MARGIN=2,max))
  W_mat=E_W(low=apply(boot_results,MARGIN=2,min),
              high=apply(boot_results,MARGIN=2,max),
              shape1=rep(1,length(C)),shape2=rep(1,length(C)),
              C_dtype = C_dtype,
              model_order = "non_standard",
              me_index=which(se_names %in% active_me)-1,
              qe_index=qe_index,
              two_fi_index=two_fi_index)
  #W_mat_list[[sim_i]]=W_mat

  start=proc.time()[3]
  I_lambda_daug=construct_des_I_lambda(x_matrix_start=design_matrix[,which(colnames(design_matrix) %in% names(beta))],
                                         x_start=design_x,W_mat=W_mat,
                                         n=nrow(design_matrix)+naug,
                                         C_dtype=C_dtype, C=C,
                                         model_order="non_standard",
                                         type= "augment",
                                         me_index=which(se_names %in% active_me)-1,
                                         qe_index=qe_index,
                                         two_fi_index=two_fi_index,
                                         n_starts=n_starts,rng=rng,num_cores=num_cores)
  end=proc.time()[3]
  print(paste("Time to execute I optimal designs",end-start))
  iaug_mat=I_lambda_daug[[1]]
  active=se_names %in% names(beta)
  active=se_names[active]
  colnames(iaug_mat)=active
  x_iaug_mat=I_lambda_daug[[2]]
  colnames(x_iaug_mat)=colnames(x_daug_mat)
  i_lambda_compound=list("I_opt_aug_x"=x_iaug_mat,"I_opt_design"=iaug_mat)
  print("#############Done with I lambda aug##################")

  start=proc.time()[3]
  w_eff=construct_des_compound(x_matrix_start=design_matrix,
                                 x_start=design_x,W_mat=W_mat,
                                 n=nrow(design_matrix)+naug,C_dtype=C_dtype,C=C,
                                 model_order="non_standard", type= "augment",
                                 D_opt_design=dsd_bayes_D_aug[[1]],
                                 D_opt_x=dsd_bayes_D_aug[[2]],
                                 I_opt_design=I_lambda_daug[[1]],
                                 me_index=which(se_names %in% active_me)-1,
                                 qe_index=qe_index,
                                 two_fi_index=two_fi_index, n_starts=n_starts_compound,
                                 D_opt_type="Bayesian",
                                 K=K,
                                 me_index_daug=me_index_daug,
                                 qe_index_daug=qe_index_daug,
                                 two_fi_index_daug=two_fi_index_daug, use_cpp=use_cpp,
                                 D_opt_tresh = D_opt_tresh,
                                 I_opt_tresh = I_opt_tresh,
                                 w_tresh = w_tresh,
                                 rng=rng,num_cores=num_cores)
  end=proc.time()[3]
  print(paste("Time to execute compound optimal designs",end-start))
  return_list=list()
  return_list[["D_aug"]]=list("D_aug_Design"=daug_mat,"D_aug_x"=x_daug_mat)
  return_list[["Boostrap_results"]]=boot_list
  return_list[["I_aug"]]=i_lambda_compound
  return_list[["compound_weff"]]=w_eff
  return_list[["K"]]=K
  return_list[["W_mat"]]=W_mat
  return_list[["x_star_beta"]]=x_star_mat
  return(return_list)
}

#'fits a GAUSS LASSO regression model with AICc for model selection to the stage 3 design and returns the regression parameters and optimal response.
#' @param design_matrix stage 3 design matrix of a given order
#' @param design_x stage 3 settings matrix
#' @param y a vector of response for the experiments
#' @param se_names column names of the stage 3 design
#' @param C_dtype data type of the factors in the design either "cont"(continuous) or "cat"(categorical)
#' @param C  list of of candidate points for each factor in the experiment
#' @param minimize a boolean value indicating whether to minimize or maximize. Default value is minimize
#' @return
#'   \itemize{
#'     \item beta- the regression parameters from the Gauss Lasso model
#'     \item x_star- the predicted optimal location for the response surface
#'     \item y_min- the predicted optimal response for the response surface
#'   }
#' @export
stage3_analysis=function(design_matrix,design_x,y,se_names,C,C_dtype,minimize=TRUE){
  #nfactor=length(C)
  n_start=length(y)
  #y_aug=apply(design_matrix[(n_start+1):nrow(design_matrix),2:(nfactor+1),drop=FALSE],MARGIN=1,y_func)
  #y_aug=apply(design_x[(n_start+1):nrow(design_x),,drop=FALSE],MARGIN=1,y_func)
  #y_aug=c(y_obs,y_aug)
  if(length(y)!=nrow(design_matrix)){
    stop("Number of rows in the design matrix and number of responses should match")
  }
  fitlasso=gamlr::gamlr(x=design_matrix[,-1],
                 y=y,family = "gaussian",
                 nlambda = 100)
  beta=as.matrix(stats::coef(fitlasso))
  X_ols=design_matrix[,beta!=0]
  beta_ols=solve(t(X_ols)%*%X_ols)%*%t(X_ols)%*%y
  active=rownames(beta_ols)[-1]
  me_index=sapply(active,FUN=function(x){length(stringr::str_split(x,"[*]",simplify = TRUE))==1},simplify = TRUE)
  active_me=active[me_index]
  active=active[!me_index]
  if(length(active)>0){
    qe_index=sapply(active,FUN=function(x){stringr::str_split(x,"[*]",simplify = TRUE)[1]==stringr::str_split(x,"[*]",simplify = TRUE)[2]},simplify = TRUE)
    active_qe=active[qe_index]
    active_2FI=active[!qe_index]
  }else{
    active_qe=c()
    active_2FI=c()
  }
  print(beta_ols)
  print(dim(beta_ols))
  parmat=sapply(C,FUN=sample,size=20,replace=TRUE)
  print(active_qe)
  print(active_2FI)
  qe_index=as.numeric(regmatches(active_qe,regexpr("[0-9]+",active_qe)))
  two_fi_index=lapply(regmatches(active_2FI,gregexpr("[0-9]+",active_2FI)),FUN=as.numeric)
  print(1+length(which(se_names %in% active_me)-1)+length(qe_index)+length(two_fi_index))
  print(two_fi_index)
  print(qe_index)
  if(minimize){
    opt_multiplier=1
  }else{
    opt_multiplier=-1
  }
  opt_vals=optimr::multistart(parmat = parmat, fn=fx,method = "L-BFGS-B",
                      b=opt_multiplier*matrix(beta_ols),model_order="non_standard",C_dtype=C_dtype,
                      me_index=which(se_names %in% active_me)-1,
                      qe_index=qe_index,
                      two_fi_index=two_fi_index,
                      lower = rep(-1,length(C)),upper=rep(1,length(C)))

  opt_vals=opt_vals[opt_vals[,"convergence"]==0,]

  x_star=as.matrix(opt_vals[which.min(opt_vals$value),1:length(C)])
  y_star=opt_multiplier*min(opt_vals$value)
  return(list("beta"=beta_ols,"x_star"=x_star[1,],"y_min"=y_star))
}


#'drops inactive factors from the design for a given heredity
#' @param drop_inactive_factors boolean value indicating whether to drop inactive factors
#' @param heredity string indicating the heredity of the model "strong", "weak" or "no". Determines which inactive effect should be dropped.
#' @param design_matrix design matrix from stage 2 augmentation
#' @param design_x settings matrix from stage 2 augmentation
#' @param beta the regression parameters from stage 2 analysis
#' @param se_names  vector of column names of the design matrix
#' @param candidate_factors vector of names of experimental factors
#' @param C list of candidate space for each experimental factor in the design
#' @param C_dtype a string vector indicating the data type of the experimental factors in the design either "cont"(continuous) or "cat"(categorical)
#' @return a list of
#'       \itemize{
#'          \item design_matrix- updated design matrix with the specified heredity
#'          \item design_x- updated factor level matrix
#'          \item me_index- updated main effects index vector for the design matrix
#'          \item qe_index- updated quadratic effects index vector for the design matrix
#'          \item two_fi_index- updated list of two factor interaction index
#'          \item C- updated candidate points list
#'          \item C_dtype- updated vector of datatypes of the experimental factors
#'          \item se_names- updated vector of column names of the design matrix
#'       }
#' @export
update_terms_stage3=function(drop_inactive_factors,heredity,design_matrix,design_x,beta,se_names,candidate_factors,C,C_dtype){
  #drop inactive terms
  if(drop_inactive_factors){
    se_names_=se_names
    active_factors=unique(unlist(sapply(names(beta),FUN=function(x){stringr::str_split(x,"[*]",simplify = TRUE)},simplify = TRUE)))
    #define a new candidate factor list
    candidate_factors_=candidate_factors[candidate_factors %in% active_factors]
    print("Updated candiate factors is")
    print(candidate_factors_)
    #define factors to delete
    del_factors=candidate_factors[!(candidate_factors %in% active_factors)]
    print("Factors to be deleted is")
    print(del_factors)
    if(length(del_factors)>0){
      for(i in del_factors){
        se_names_=se_names_[sapply(stringr::str_split(se_names_,"[*]"),FUN=function(x){all(!(x %in% i))})]
      }
    }
    #remove terms which does not meet the weak/strong heredity heirarchy
    C_=C[names(C)%in%candidate_factors_]
    C_dtype_=C_dtype[names(C)%in%candidate_factors_]
  }else{
    C_=C
    C_dtype_=C_dtype
    se_names_=se_names
    candidate_factors_=candidate_factors
  }
  nfactor=length(C_)
  two_fi_index=utils::combn(1:nfactor,2,simplify = FALSE)
  me_index=sapply(names(beta)[-1],FUN=function(x){length(stringr::str_split(x,"[*]",simplify=TRUE))==1})
  me_index=names(beta)[-1][me_index]
  active_me=me_index
  print(active_me)
  me_index_drop=names(C_)[!(names(C_)%in%me_index)]
  me_index=which(names(C_)%in%me_index)
  qe_index=1:nfactor
  #active=se_names_ %in% names(beta)
  #active=se_names_[active]
  active=names(beta)
  active=active[-1]
  print(active)
  if(heredity=="weak"){
    del_=two_fi_index[sapply(two_fi_index,FUN=function(x){!(any(x %in% which(se_names_[-1]%in%active_me)))})]
    del=sapply(del_,function(x){paste0(names(C_)[x[1]],"*",names(C_)[x[2]])})
    if(sum(del %in% active)>0){
      del_=del_[-which(del %in% active)]
      del=del[-which(del %in% active)]
    }
    print("Terms to be deleted under weak heredity")
    print(del)
    #print("Active ME")
    #print(active_me_list[[sim_i]])
    #update se_names with the terms deleted
    if(length(del)>0){
      se_names_=se_names_[-which(se_names_ %in% del)]
      two_fi_index=two_fi_index[which(!(two_fi_index %in% del_))]
    }
    print(paste("Length of weak interaction terms",length(two_fi_index)))
    print(length(se_names_))
    if(length(me_index)<length(C_)){se_names_=se_names_[-which(se_names_ %in% me_index_drop)]}
    print("Terms in weak heredity model")
    print(se_names_)
    print(length(se_names_))
  }else if(heredity=="strong"){
    del=two_fi_index[sapply(two_fi_index,FUN=function(x){!(all(x %in% which(se_names_[-1]%in%active_me)))})]
    del=sapply(del,function(x){paste0(names(C_)[x[1]],"*",names(C_)[x[2]])})
    del=del[!(del %in% active)]
    print("Terms to be deleted under weak heredity")
    print(del)
    #print("Active ME")
    #print(active_me_list[[sim_i]])
    two_fi_index=two_fi_index[sapply(two_fi_index,FUN=function(x){all(x %in% which(se_names_[-1]%in%active_me))})]
    print(paste("Length of weak interaction terms",length(two_fi_index)))
    #update se_names with the terms deleted
    if(length(del)>0){se_names_=se_names_[-which(se_names_ %in% del)]}
    if(length(me_index)<length(C_)){se_names_=se_names_[-which(se_names_ %in% me_index_drop)]}
    print("Terms in strong heredity model")
    print(se_names_)
    print(length(se_names_))
  }
  #se_names_list[[sim_i]]=c("intercept",candidate_factors_,se_names_,
  #           as.vector(sapply(candidate_factors_,FUN=function(x){paste0(x,"*",x)},simplify = TRUE)))
  design_x_=design_x[,candidate_factors_]
  design_matrix_=t(apply(design_x_,MARGIN=1,FUN=f_x,model_order="non_standard",me_index=me_index,qe_index=qe_index,two_fi_index=two_fi_index,C_dtype=C_dtype_))
  print(dim(design_matrix_))
  print(length(se_names_))
  colnames(design_matrix_)=se_names_
  #print(which(colnames(bayes_d_aug_list[[sim_i]]) %in% se_names_list[[sim_i]]))
  #bayes_d_aug_list[[sim_i]]=bayes_d_mat_upd
  #me_index_daug_list[[sim_i]]=me_index
  #qe_index_daug_list[[sim_i]]=qe_index
  #two_fi_index_daug_list[[sim_i]]=two_fi_index
  return(list("design_matrix"=design_matrix_,"design_x"=design_x_,"me_index"=me_index,"qe_index"=qe_index,"two_fi_index"=two_fi_index,
              "C"=C_,"C_dtype"=C_dtype_,"se_names"=se_names_))
}

#' function to transform the design space
#' @param design_x a settings matrix
#' @param std a boolean value, when TRUE transforms design from (a,b) to (-1,1) and when FALSE transforms (-1,1) to (a,b)
#' @param lower a vector of lower levels for the factors in the experiment
#' @param upper a vector of upper levels for the factors in the experiment
#' @description a function that transforms the design space from (a,b) to (-1,1) or from (-1,1) to (a,b)
#' @return a transformed matrix of design factor levels
#' @example
#' @export
transform_designspace=function(design_x,std=TRUE, lower,upper){
  if(length(lower)!=length(upper)){
    stop("The length of lower and upper vector should be the same")
  }else if(length(lower)!=ncol(design_x)){
    stop("The length of lower and upper vectors should be the same as number of factors in the design")
  }
  #if std=TRUE transform to [-1,1] design space
  if(std){
    transformed_std_x=t(apply(design_x,MARGIN=1,FUN=function(x){(x-((upper+lower)/2))/((upper-lower)/2)}))
    return(transformed_std_x)
  }else{
    #if std=FALSE transform to [a,b] design space
    transformed_x=t(apply(design_x,MARGIN=1,FUN=function(x){x*((upper-lower)/2)+((upper+lower)/2)}))
    return(transformed_x)
  }
}

#' function to transform the factor names to a standard notation and backtransform
#' @param design a design matrix
#' @param backtransform a boolean value, when TRUE transforms the column names from standard notation "X1"..."Xn" to factor names
#' @param colnames_vector a vector mapping the factor names to a standard notation "X1"..."Xn". To be provided when backtransform=TRUE
#' @param factor_levels a boolean value, when TRUE indicates that the input matrix is the settings matrix and when FALSE indicates that the input matrix is the design matrix. Design matrix is assumed to contain intercept as the first column.
#' @description function to transform the factor names to a standard notation "X1"..."Xn" and backtransform
#' @return a list of:
#'      \itemize{
#'        \item colnames_vector- a vector mapping the factor names to a standard notation "X1"..."Xn"
#'        \item design- the matrix with coded column names
#'      }
#' @example
#' @export
factornames_standardize=function(design,backtransform=FALSE,colnames_vector=c(),factor_levels=FALSE){
  if(backtransform){
    if(factor_levels==FALSE){
      #assuming the first column is "intercept"
      colnames_vector_=c("intercept",colnames_vector)
      names(colnames_vector_)[1]="intercept"
      x_names=stringr::str_split(colnames(design),"[*]")
      x_names=lapply(x_names,FUN=function(x){sapply(x,FUN=function(y){colnames_vector_[which(names(colnames_vector_)%in%y)]})})
      x_names=sapply(x_names,FUN=function(x){if(length(x)>1){paste0(x[1],"*",x[2])}else{x}})
      colnames(design)=x_names
    }else{
      indx=sapply(colnames(design),FUN=function(x){colnames_vector[which(names(colnames_vector)%in%x)]})
      colnames(design)=indx
    }
  }else{
    colnames_vector=colnames(design)
    print(colnames_vector)
    names(colnames_vector)=paste0("X",1:ncol(design))
    colnames(design)=names(colnames_vector)
  }
  return(list("colnames_vector"=colnames_vector,"design"=design))
}

