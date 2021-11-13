#' an internal function to create interaction matrix for model selection
#' @description constructs the interaction matrix for model selection
#' @param C_active_dtype a string vector containing the datatype of the active ME
#' @param C_active_names a string vector names of active experimental factor
#' @param C_non_active_names a string vector names of inactive experimental factors. It is not NA only when a weak heredity interaction matrix is desired.
#' @param x_matrix_me a matrix of ME for the design
#' @return a matrix of interaction terms
#' @export
create_interaction_matrix=function(C_active_dtype, C_active_names,
                                   C_non_active_names,
                                   x_matrix_me){
  #creates an interaction matrix for a given me and active and non_active terms
  #for active terms a full second order model is created
  #for non-active terms 2fi with active ME is created
  #returns a matrix
  active_me=x_matrix_me[,colnames(x_matrix_me)%in%C_active_names,drop=FALSE]
  #print(C_active_names)
  #print(colnames(x_matrix_me))
  #print(active_me)

  interaction_colnames=c()
  interaction_matrix=c()

  if(class(active_me)[1]!="matrix"){
    print("Only one ME is active")
    active_me=matrix(active_me)
  }

  p=ncol(active_me)
  for(i in 1:p){
    interaction_matrix=cbind(interaction_matrix,active_me[,i]*active_me[,i:p])
    interaction_colnames=c(interaction_colnames,paste0(colnames(active_me)[i],
                                                       "*",colnames(active_me)[i:p]))
  }

  del=paste0(colnames(active_me)[C_active_dtype=="cat"],"*",
             colnames(active_me)[C_active_dtype=="cat"])
  interaction_matrix=interaction_matrix[,!(interaction_colnames %in% del),drop=FALSE]
  interaction_colnames=interaction_colnames[!(interaction_colnames %in% del)]

  if(length(C_non_active_names)>0){
    non_active_me=x_matrix_me[,colnames(x_matrix_me)%in%C_non_active_names,drop=FALSE]
    #print(class(non_active_me))
    for(i in 1:ncol(active_me)){
      interaction_matrix=cbind(interaction_matrix,active_me[,i]*non_active_me)
      interaction_colnames=c(interaction_colnames,paste0(colnames(active_me)[i],
                                                         "*",colnames(non_active_me)))
    }
  }
  colnames(interaction_matrix)=interaction_colnames
  #print("colnames of interactio matrix")
  #print(colnames(interaction_matrix))
  return(interaction_matrix)
}

#' function to fit all possible regression, not recommended for large experiments
#' @description function to fit all possible regression, not recommended for large experiments
#' @param y a vector of response
#' @param XME a matrix of ctive Main Effects from stage 1 JN method
#' @param XSE a mtraix of second order terms
#' @param maxterms a numeric value indicating the max terms in the design
#' @return a list of:
#'      \itemize{
#'        \item bestModel- the matrix of regressors for the model chosen from all possible regression
#'        \item bestSSE- the SSE of the chosen model
#'        \item minAICc- AICCc of the chosen model
#'        \item best_beta- the parameter estimates of the chosen model
#'      }
#' @export
allpossibleregression = function(y, XME, XSE, maxterms){
  #  Best subsets using AICc
  # y, y2nd, XME, fq, maxterms, mseME, errorDF, pME
  #save number of rows and cols
  n = dim(XSE)[1]
  nc = dim(XSE)[2]
  # browser
  #create matrix of ME with intercept
  pME=ncol(XME)
  XME = cbind(matrix(1,n,1),XME)
  #number of columns of ME
  nTermsME = dim(XME)[2]
  #compute betaME
  betaME = solve(t(XME) %*% XME) %*% t(XME) %*% y
  betaME = matrix(betaME,nTermsME,1)
  #compute residual ME
  rME =  y - XME %*% betaME
  #compute SSE me
  sseME = as.numeric(t(rME) %*% rME)
  #p is adjusted for the AIC formula (p+1)
  p = 2 + pME
  # compute aicc for ME model
  aicc0 = n*(log(2*pi*sseME/n) + 1) + 2*p + 2*p*(p+1)/(n-p-1)
  minAICc = aicc0
  bestModel = c()
  bestSSE = sseME
  best_beta=betaME
  for (nterms in 1:maxterms){
    factors = 1:nc
    mods = t(utils::combn(factors, nterms))
    nmods = dim(mods)[1]
    for (midx in 1:nmods){
      x2nd = XSE[, mods[midx,]]
      x = cbind(XME,x2nd)
      nColsInx = dim(x)[2]
      #if (rankMatrix(x) == nColsInx) {
      if (rcond(t(x)%*%x) > 1e-8){ # make sure the x matrix is not singular
        b = solve( t(x)  %*%  x ) %*% (t(x)  %*%  y)
        r = y - x  %*%  b
        sse = t(r)  %*%  r
        if (sse > 1e-8) {  # do not include saturated models
          p = nterms + pME + 2
          AICc = n*(log(2*pi*sse/n) + 1) + 2*p + 2*p*(p+1)/(n-p-1)
          if (AICc < minAICc){
            minAICc = AICc
            bestSSE = sse
            bestModel = mods[midx,]
            best_beta=b
          }
        }
      }
    }
  }
  return(list(bestModel, bestSSE, minAICc,best_beta))
}

#' function for stepwise AICc model selection
#' @description function for stepwise AICc model selection on the second order terms in the JN model selection method
#' @param y a vector of response
#' @param XME a matrix of active Main Effects from stage 1 JN method
#' @param XSE a mtraix of second order terms
#' @param maxterms a numeric value indicating the max terms in the design
#' @return a list of:
#'      \itemize{
#'        \item x_se- the matrix of second order terms for the model chosen from stepwise AICc
#'        \item bestSSE- the SSE of the chosen model
#'        \item minAICc- AICCc of the chosen model
#'        \item best_beta- the parameter estimates of the chosen model
#'      }
#' @export
stepAICc = function(y, XME, XSE, maxterms){
  #  Best subsets using AICc
  # y, y2nd, XME, fq, maxterms, mseME, errorDF, pME
  #save number of rows and cols
  n = dim(XSE)[1]
  nc = dim(XSE)[2]
  # browser
  #create matrix of ME with intercept
  pME=ncol(XME)
  XME = cbind(matrix(1,n,1),XME)
  #print(XME)
  #number of columns of ME
  nTermsME = dim(XME)[2]
  #compute betaME
  betaME = solve(t(XME) %*% XME) %*% t(XME) %*% y
  betaME = matrix(betaME,nTermsME,1)
  #print(betaME)
  #compute residual ME
  rME =  y - XME %*% betaME
  #compute SSE me
  sseME = as.numeric(t(rME) %*% rME)
  #p is adjusted for the AIC formula (p+1)
  p = 2 + pME
  #print(paste("p",p))
  # compute aicc for ME model
  aicc0 = n*(log(2*pi*sseME/n) + 1) + 2*p + 2*p*(p+1)/(n-p-1)
  minAICc = aicc0
  #print(minAICc)
  #print(n*(log(2*pi*sseME/n)))
  #print(2*p + 2*p*(p+1)/(n-p-1))
  x_se_local = c()
  bestSSE = sseME
  best_beta=betaME
  x_se=c()
  xse_index=colnames(XSE)
  epsilon=1e4
  nse=0
  #print("xse_index")
  #print(xse_index)
  while(nse<=maxterms & epsilon>0){
    #print("Starting while loop")
    prev=minAICc
    aicc_list=c()
    for(term in xse_index){
      #print(term)
      x_se_local = cbind(x_se,XSE[, term,drop=FALSE])
      #print(x_se_local)
      x = cbind(XME,x_se_local)
      nColsInx = dim(x)[2]
      #if (rankMatrix(x) == nColsInx) {
      if (rcond(t(x)%*%x) > 1e-8){ # make sure the x matrix is not singular
        b = solve( t(x)  %*%  x ) %*% (t(x)  %*%  y)
        r = y - x  %*%  b
        sse = t(r)  %*%  r
        if (sse > 1e-8) {  # do not include saturated models
          p = (nse+1) + pME + 2
          AICc = n*(log(2*pi*sse/n) + 1) + 2*p + 2*p*(p+1)/(n-p-1)
          #print(AICc)
          aicc_list=c(aicc_list,AICc)
        }
      }else{print("Singular, hence model not fit")}
    }
    #print(aicc_list)
    #print(min(aicc_list))
    #print(which.min(aicc_list))
    #print(xse_index[which.min(aicc_list)])
    #print(XSE[,xse_index[which.min(aicc_list)],drop=FALSE])
    #print(aicc_list)
    #print(minAICc)
    if (min(aicc_list) < minAICc){
      minAICc = min(aicc_list)
      bestSSE = sse
      x_se = cbind(x_se,XSE[,xse_index[which.min(aicc_list)],drop=FALSE])
      x = cbind(XME,x_se)
      best_beta=solve( t(x)  %*%  x ) %*% (t(x)  %*%  y)
    }
    #print(x_se)
    #track change in minAICc
    epsilon=prev-minAICc
    #print(paste("epsilon",epsilon))
    #remove se_best from index
    if(epsilon>0){
      rm_index=which(xse_index==colnames(x_se)[ncol(x_se)])
      xse_index=xse_index[-rm_index]
      nse=ncol(x_se)
      #print(xse_index)
    }
  }
  return(list(x_se, bestSSE, minAICc,best_beta))
}

#' a function for Jones and Nachtsheim (JN) model selection method
#' @description a function for Jones and Nachtsheim model selection method. Works only for foldover designs.
#' @param dsd_matrix a DSD matrix for model selection
#' @param y a vector of response
#' @param alpha a numeric value the signficance level for model selection, set at the value recommended in the JN 2017 paper
#' @param heredity a string indicating the heredity for model selection. One of "Strong", "weak", and "no". Heredity determines the number of second order terms considered in model selection.
#' @param fake_factor_index a vector indeicating the index of facle factor terms in the design matrix
#' @param C_dtype a string vector indicating the data type of the experimental factors in the design either "cont"(continuous) or "cat"(categorical)
#' @param model_selection_type a string indicating the type of model selection method adopted for second order terms. One of "JN_org"- method proposed in the paper, "LASSO"- LASSO with AICc, "APR"- all possible regression, not recommended for large designs and "StepAICc"- stepwise AICc
#' @return a list of:
#'    \itemize{
#'      \item ME Index- index of active Main Effects
#'      \item Beta_me- parameter estimates of the main effects
#'      \item Tstat- t-statistic of the hypothesis test of regression parameters
#'      \item p-value- p-value of the significance test of regression parameters
#'      \item active_se- regression matrix of active second order terms
#'      \item beta_se- parameter estimates of the second order terms
#'    }
#' @export
JN_model_selection=function(dsd_matrix,y,alpha=0.2,heredity="strong",
                            fake_factor_index,
                            C_dtype,model_selection_type="JN_org"){
  y=matrix(y)
  ff_matrix=dsd_matrix[,fake_factor_index]
  x_matrix_me=dsd_matrix[,-fake_factor_index]
  n=nrow(dsd_matrix)
  #compute se and dof_se
  if(nrow(dsd_matrix)>(2*ncol(x_matrix_me)+1)){#ensure that there are fake factors, ff_df for computing sse for ME screening
    #center_runs=apply(x_matrix_me,MARGIN=1,FUN = function(x){all(x==0)})
    unique_rows=unique(x_matrix_me)
    rep_index=t(apply(unique_rows,MARGIN=1,FUN = function(x){apply(x_matrix_me,MARGIN=1,FUN=function(y){all(y==x)})}))
    rep_index=apply(rep_index,MARGIN=1,FUN=function(z){which(z==TRUE)})
    rep_index=rep_index[sapply(rep_index, FUN=function(x){length(x)>1})]
    if(class(rep_index)=="list"){
      center_runs=rep_index[[1]]
    }else{center_runs=rep_index}
    #print(center_runs)
    #print(length(center_runs))
    sf=t(y)%*%ff_matrix%*%solve(t(ff_matrix)%*%ff_matrix)%*%t(ff_matrix)%*%y
    sf=sf/ncol(ff_matrix)
    #if(sum(center_runs)>1 && !("cat"%in%C_dtype)){
    if(length(center_runs)>1){
      print("updating the se estimate to include center runs")
      print("Center Run index")
      print(center_runs)
      y_bar_sc=mean(y[center_runs,])
      #sc=sum((y[center_runs]-y_bar_sc)**2)/((sum(center_runs)-1))
      sc=sum((y[center_runs]-y_bar_sc)**2)/((length(center_runs)-1))
      #se=((sum(center_runs)-1)*sc+ncol(ff_matrix)*sf)/(sum(center_runs)-1+ncol(ff_matrix))
      se=((length(center_runs)-1)*sc+ncol(ff_matrix)*sf)/(length(center_runs)-1+ncol(ff_matrix))
    }else{
      se=sf
    }

    #identify active ME
    #dsd_matrix=cbind(rep(1,nrow(dsd_matrix)),dsd_matrix)
    intercept=(t(matrix(rep(1,nrow(dsd_matrix))))%*%y/n)[1]
    #print(intercept)
    beta_me=solve(t(dsd_matrix)%*%dsd_matrix)%*%t(dsd_matrix)%*%(y-intercept)
    y_me=dsd_matrix%*%beta_me
    y_se=(y-intercept)-y_me
    beta_me=beta_me[-fake_factor_index]
    s_beta_me=se*diag(solve(t(dsd_matrix)%*%dsd_matrix))
    s_beta_me=s_beta_me[-fake_factor_index]
    t_stat=abs(beta_me)/sqrt(s_beta_me)
    p_value=stats::pt(t_stat,df=n-length(beta_me),lower.tail = FALSE)
    active_me_index=which(p_value<=0.05)
    #active_me_index=active_me_index-1
    #update se and dof_se
    ff_matrix=cbind(ff_matrix,x_matrix_me[,-active_me_index])
    sf=t(y)%*%ff_matrix%*%solve(t(ff_matrix)%*%ff_matrix)%*%t(ff_matrix)%*%y
    sf=sf/ncol(ff_matrix)
    #if(sum(center_runs)>1 & !("cat"%in%C_dtype)){
    if(length(center_runs)>1){
      print("updating the se estimate to include center runs")
      print("Center Run index")
      print(center_runs)
      y_bar_sc=mean(y[center_runs,])
      #sc=sum((y[center_runs]-y_bar_sc)**2)/((sum(center_runs)-1))
      sc=sum((y[center_runs]-y_bar_sc)**2)/((length(center_runs)-1))
      #se=((sum(center_runs)-1)*sc+ncol(ff_matrix)*sf)/(sum(center_runs)-1+ncol(ff_matrix))
      se=((length(center_runs)-1)*sc+ncol(ff_matrix)*sf)/(length(center_runs)-1+ncol(ff_matrix))
      #sspe_df=sum(center_runs)-1+ncol(ff_matrix)
      sspe_df=length(center_runs)-1+ncol(ff_matrix)
    }else{
      se=sf
      sspe_df=ncol(ff_matrix)
    }
    #print(p_value)
    #print(beta_me)
    #print(t_stat)
    #print(active_me_index)
    #account for heredity and non-heredity models
    if(length(active_me_index)<=1){
      print(paste("number of active ME is", length(active_me_index),"model is automatically changed to no heredity"))
      heredity="no"
    }
    if(heredity=="strong"){
      #create interaction matrix
      C_active_names=colnames(x_matrix_me)[active_me_index]
      #C_non_active_names=names(x_matrix_me)[-active_me_index]
      interaction_matrix=create_interaction_matrix(C_active_dtype = C_dtype,
                                                   C_active_names = C_active_names,
                                                   C_non_active_names = c(),
                                                   x_matrix_me = x_matrix_me)
      #center interaction matrix
      interaction_matrix_local=apply(interaction_matrix,MARGIN=2,FUN = function(x){x-mean(x)})
      if(class(interaction_matrix_local)[1]!="matrix"){
        interaction_matrix_local=as.matrix(interaction_matrix_local)
      }
    }else if(heredity=="weak"){
      #create interaction matrix
      #print("active_me_index")
      #print(active_me_index)
      C_active_names=colnames(x_matrix_me)[active_me_index]
      #print(C_active_names)
      #print(colnames(x_matrix_me))
      C_non_active_names=colnames(x_matrix_me)[-active_me_index]
      interaction_matrix=create_interaction_matrix(C_active_dtype = C_dtype,
                                                   C_active_names = C_active_names,
                                                   C_non_active_names = C_non_active_names,
                                                   x_matrix_me = x_matrix_me)
      #center interaction matrix
      interaction_matrix_local=apply(interaction_matrix,MARGIN=2,FUN = function(x){x-mean(x)})
    }else if(heredity=="no"){
      #create interaction matrix
      C_active_names=colnames(x_matrix_me)
      interaction_matrix=create_interaction_matrix(C_active_dtype = C_dtype,
                                                   C_active_names = C_active_names,
                                                   C_non_active_names = c(),
                                                   x_matrix_me = x_matrix_me)
      #center interaction matrix
      interaction_matrix_local=apply(interaction_matrix,MARGIN=2,FUN = function(x){x-mean(x)})
    }else if(heredity=="noheredity-quadratic"){
      C_active_names=colnames(x_matrix_me)[active_me_index]
      #print(C_active_names)
      interaction_matrix=create_interaction_matrix(C_active_dtype = C_dtype,
                                                   C_active_names = C_active_names,
                                                   C_non_active_names = c(),
                                                   x_matrix_me = x_matrix_me)
      if(sum(!(colnames(x_matrix_me)%in%C_active_names))>0){
        nonactive_qe=x_matrix_me[,!(colnames(x_matrix_me)%in%C_active_names),drop=FALSE]**2
        #print(nonactive_qe)
        nonactive_names=colnames(x_matrix_me)[!(colnames(x_matrix_me)%in%C_active_names)]
        colnames(nonactive_qe)=paste0(nonactive_names,"*",nonactive_names)
        interaction_matrix=cbind(interaction_matrix,nonactive_qe)
        #print(interaction_matrix)
      }else{
        print("All me are active")
      }
      #center interaction matrix
      interaction_matrix_local=apply(interaction_matrix,MARGIN=2,FUN = function(x){x-mean(x)})
    }
    #interaction_matrix_local=interaction_matrix
    #identify active 2FI
    c=ncol(dsd_matrix)
    if(model_selection_type=="JN_org"){
      #y_se=y_se-matrix(rep(1,nrow(dsd_matrix)))%*%intercept
      TSS=t(y_se)%*%y_se
      F_stat=(TSS/c)/se
      if(F_stat<stats::qf(alpha,c,ncol(ff_matrix))){
        return(list(active_me_index,c(intercept,beta_me),t_stat,p_value,c(),c()))
      }else{
        i=0
        active_se=c()
        beta_se=c()
        #print(c)
        while(F_stat>=stats::qf(alpha,c-i,sspe_df)&&(i+1)<=c/2&&
              ncol(interaction_matrix_local)>0){
          i=i+1
          R=c()
          for(j in 1:ncol(interaction_matrix_local)){
            x_se=cbind(active_se,interaction_matrix_local[,j,drop=FALSE])
            #print(dim(x_se))
            #print(interaction_matrix_local)
            #print(x_se)
            beta_se=solve(t(x_se)%*%x_se)%*%t(x_se)%*%y_se
            #print(dim(t(y_se-(x_se%*%beta_se))))
            #print(dim(y_se-(x_se%*%beta_se)))
            #print(t(y_se-(x_se%*%beta_se))%*%(y_se-(x_se%*%beta_se)))
            R=c(R,t(y_se-(x_se%*%beta_se))%*%(y_se-(x_se%*%beta_se)))
          }
          #print(R)
          active_se=cbind(active_se,interaction_matrix_local[,which.min(R)[1],drop=FALSE])
          interaction_matrix_local=interaction_matrix_local[,-which.min(R)[1],drop=FALSE]
          F_stat=(R[which.min(R)]/(c-i))/se
          beta_se=solve(t(active_se)%*%active_se)%*%t(active_se)%*%y_se
          active_name=colnames(active_se)[ncol(active_se)]
          #print(interaction_matrix_local)
          #print(ncol(interaction_matrix_local))
          #print(acive_se)
          #print(beta_se)
          #print(F_stat)
        }
        print("Completed model selection")
        #print(active_se)
        #print(colnames(active_se))
        #print(interaction_matrix[,colnames(active_se),drop=FALSE])
        #print(apply(interaction_matrix[,colnames(active_se)],MARGIN=2,FUN=mean))
        #print(mean(interaction_matrix[,colnames(active_se)])*beta_se)
        if(!is.null(active_se)){
          intercept=intercept-sum(apply(interaction_matrix[,colnames(active_se),drop=FALSE],MARGIN=2,FUN=mean)*beta_se)
        }
        return(list("ME Index"=active_me_index,"Beta_me"=c(intercept,beta_me),"Tstat"=t_stat,"p-value"=p_value,
                    "active_se"=active_se,"beta_se"=beta_se))
      }
    }else if(model_selection_type=="APR"){
      apr_results=allpossibleregression(y=y, XME=x_matrix_me[,active_me_index,drop=FALSE], XSE=interaction_matrix, maxterms=c/2)
      #print(apr_results[[4]])
      return(list("ME Index"=active_me_index,"Beta_me"=apr_results[[4]][1:(length(active_me_index)+1),],"Tstat"=t_stat,"p-value"=p_value,
                  "active_se"=interaction_matrix[,apr_results[[1]]],"beta_se"=apr_results[[4]][(length(active_me_index)+2):nrow(apr_results[[4]]),]))
    }
    else if(model_selection_type=="LASSO"){
      if(length(active_me_index)!=0){
        x=cbind(x_matrix_me[,active_me_index],interaction_matrix)
        lasso_model=gamlr::gamlr(x=x, y=y,family = "gaussian", nlambda = 100, free=1:length(active_me_index))
        b=as.matrix(stats::coef(lasso_model,corrected = TRUE))
        #print(b)
      }else{
        x=interaction_matrix
        lasso_model=gamlr::gamlr(x=x, y=y,family = "gaussian", nlambda = 100)
        b=as.matrix(stats::coef(lasso_model,corrected = TRUE))
      }
      x_ols=x[,which(b[-1,]!=0)]
      x_ols=cbind(matrix(1,nrow=nrow(x_ols),ncol=1),x_ols)
      if(rcond(x_ols)<1e-15){
        print("Lasso solution is over specified, hence switching to setpwise AIC")
        intermediate=JN_model_selection(dsd_matrix=dsd_matrix,y=y,alpha=alpha,heredity=heredity,fake_factor_index=fake_factor_index,
                                        C_dtype=C_dtype,model_selection_type="stepAICc")
        #print(intermediate)
        return(intermediate)
      }
      b_ols=solve(t(x_ols)%*%x_ols)%*%t(x_ols)%*%y
      #print(b_ols)
      if(nrow(b_ols)==length(active_me_index)+1){
        active_se=c(NA)
        beta_se=c(NA)
      }else{
        active_se=x_ols[,(length(active_me_index)+2):nrow(b_ols)]
        beta_se=b_ols[(length(active_me_index)+2):nrow(b_ols),]
      }
      return(list("ME Index"=active_me_index,"Beta_me"=b_ols[1:(length(active_me_index)+1),],"Tstat"=t_stat,"p-value"=p_value,
                  "active_se"=active_se,"beta_se"=beta_se))
    }
    else if(model_selection_type=="stepAICc"){
      apr_results=stepAICc(y=y, XME=x_matrix_me[,active_me_index,drop=FALSE], XSE=interaction_matrix, maxterms=c/2)
      #print(apr_results[[4]])
      if(!is.null(apr_results[[1]])){
        beta_se=apr_results[[4]][(length(active_me_index)+2):nrow(apr_results[[4]]),]
        active_se=apr_results[[1]]
      }else{
        beta_se=c(NA)
        active_se=c(NA)
      }
      return(list("ME Index"=active_me_index,"Beta_me"=apr_results[[4]][1:(length(active_me_index)+1),],"Tstat"=t_stat,"p-value"=p_value,
                  "active_se"=active_se,"beta_se"=beta_se))
    }
  }
}

# dsd_six_factor=matrix(c(0,1,1,1,1,1,1,1,
#                         0,-1,-1,-1,-1,-1,-1,-1,
#                         1,0,1,1,-1,1,-1,-1,
#                         -1,0,-1,-1,1,-1,1,1,
#                         1,-1,0,1,1,-1,1,-1,
#                         -1,1,0,-1,-1,1,-1,1,
#                         1,-1,-1,0,1,1,-1,1,
#                         -1,1,1,0,-1,-1,1,-1,
#                         1,1,-1,-1,0,1,1,-1,
#                         -1,-1,1,1,0,-1,-1,1,
#                         1,-1,1,-1,-1,0,1,1,
#                         -1,1,-1,1,1,0,-1,-1,
#                         1,1,-1,1,-1,-1,0,1,
#                         -1,-1,1,-1,1,1,0,-1,
#                         1,1,1,-1,1,-1,-1,0,
#                         -1,-1,-1,1,-1,1,1,0,
#                         0,0,0,0,0,0,0,0),ncol=8,byrow=TRUE)
#
# colnames(dsd_six_factor)=c("X1","X2","X3","X4","X5","X6","X7","X8")
# fake_factor_index = c(7,8)

# y_func=function(x){
#   return(9+2.3*x[1]+2.7*x[3]-2.6*x[4]+2.3*x[5]+2.9*x[6]+
#            -3.1*x[4]*x[5]+2.9*x[3]*x[1]-2.3*x[1]*x[1]-
#            2.2*x[5]*x[6]-2.7*x[6]*x[6]+rnorm(1))
# }

# y=apply(dsd_six_factor,MARGIN=1,y_func)

# JN_result=JN_model_selection(dsd_matrix = dsd_six_factor, y=y, alpha=0.2,
#                              heredity = "strong",fake_factor_index = fake_factor_index,
#                              C_dtype = rep("cont",8))


# dat=read.csv("C:/Users/sunde/Downloads/designData2020-11-14.csv")
# dsd_matrix=dat[,-ncol(dat)]
# dsd_matrix=as.matrix(dsd_matrix)
# y=dat[,ncol(dat)]
# JN_model_selection(dsd_matrix = dsd_matrix, y=y, fake_factor_index = c(7,8),C_dtype = rep("cont",8))
# JN_model_selection(dsd_matrix = dsd_matrix, y=y, fake_factor_index = c(7,8),C_dtype = rep("cont",8),model_selection_type="APR")
# JN_model_selection(dsd_matrix = dsd_matrix, y=y, fake_factor_index = c(7,8),C_dtype = rep("cont",8),model_selection_type="LASSO")
# JN_model_selection(dsd_matrix = dsd_matrix, y=y, fake_factor_index = c(7,8),C_dtype = rep("cont",8),model_selection_type="stepAICc")
