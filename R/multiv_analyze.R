#'Multivariate analysis
#'@description Perform multivariate analysis.
#'@usage multiv_analyze(METBObj, method="pca", scale="center")
#'@param METBObj METBObj object contains list of data.
#'@param method name of multivariate analysis method. Choose one from the list: pca, pls, opls. Default is pca.
#'@param scale text indicating a scaling method for parameter \code{scaleC} of \code{ropls::opls()}: no centering nor scaling ('none'),
#'mean-centering only ('center') [default], mean-centering and pareto scaling ('pareto'), or mean-centering and auto scaling ('standard').
#'@return a list of the following components:
#'
#'model_summary = a data frame with the model overview.
#'
#'score_val = a numerical matrix of x scores.
#'
#'loading_val = a numerical matrix of x loadings.
#'
#'oscore_val = a numerical matrix of orthogonal scores.
#'
#'oloading_val = a numerical matrix of orthogonal loadings.
#'
#'vip_val = a numerical vector of VIP.
#'
#'ovip_val = a numerical vector of variable importance for orthogonal modeling.
#'
#'details = a list of analysis details: itestMethod, scale.
#'
#'@author Kwanjeera W \email{kwanjeera.wan@@mahidol.ac.th}
#'@references Thevenot EA, et al. (2015) ropls.
#'@seealso \code{\link[ropls:opls]{ropls::opls()}}
#'@examples
#'#sugar_dt = set_input_obj(sugar, 1,2,5)
#'#out=multiv_analyze(sugar_dt)
#'@export
multiv_analyze <- function(METBObj, method="pca", scale="center"){
  #Check argument
  tmparg <- try(method <- match.arg(method, c("pca","pls","opls"), several.ok = FALSE), silent = TRUE)
  if (class(tmparg) == "try-error") {
    cat("\nERROR! Argument 'method' is not valid, choose one from the list: pca,pls,opls.\nData was not analyzed.\n")
    return(FALSE)
  }
  tmparg <- try(scale <- match.arg(scale, c("none","center","pareto","standard"), several.ok = FALSE), silent = TRUE)
  if (class(tmparg) == "try-error") {
    cat("\nERROR! Argument 'scale' is not valid, choose one from the list: none,center,pareto,standard.\nData was not analyzed.\n")
    return(FALSE)
  }
  if(sum(is.na(METBObj$X)) > 0){#Data contains missing values
    cat("\nThe data contains missing values. Data was not analyzed.\n")
    return(FALSE)
  }

  #Initialize parameters
  multiv_result = list(); methodls = list(); #Working data
  dat = METBObj$X; #Working data
  cat("\nChecking essential parameters ...")
  inst_pkg = NULL
  if(!requireNamespace("ropls", quietly = TRUE)){#check and install required package
    cat("\nMissing the required package 'ropls', trying to install the package ...\n")
    inst_pkg = install_pkgs('ropls')
  }
  if(length(unlist(inst_pkg))){
    cat("\nERROR! Could not install the required package 'ropls'. Data was not analyzed.\n")
    return(FALSE)
  }
  #Check type of category/factor column
  if (is.factor(METBObj$Y)) {
    Y = METBObj$Y #working data
    cat('\nFind a factor of ',nlevels(Y),' levels.',sep='')
  }else if(is.numeric(METBObj$Y)){
    Y = METBObj$Y #working data
    cat('\nFind a numeric vector. Regression analysis will be performed.')
  }else{
    Y = as.factor(METBObj$Y) #working data
    cat('\nThe category/factor column is', class(METBObj$Y), '.\nConverting to a factor of ',
        nlevels(Y),' levels. ','Classification analysis will be performed.',sep='')
  }
  #Multivariate analysis
  cat("\n\nPerforming multivariate analysis ...\n")
  if (method == "pca"){#PCA
    out_data = tryCatch({
      ropls::opls(x=dat, fig.pdfC="none", scaleC=scale, crossvalI=6)
    },
    error=function(e){
      cat(e$message)
      #message(e)
      cat("\nData was not analyzed.\n")
      NULL
    })
  }
  if (method == "pls"){#PLS
    out_data = tryCatch({
      ropls::opls(x=dat, y=Y, fig.pdfC="none", scaleC=scale, crossvalI=6)
    },
    error=function(e){
      cat(e$message)
      #message(e)
      cat("\nData was not analyzed.\n")
      NULL
    })
  }
  if (method == "opls"){#OPLS
    out_data = tryCatch({
      ropls::opls(x=dat, y=Y, predI = 1, orthoI = 2, fig.pdfC="none", scaleC=scale, crossvalI=6)
    },
    error=function(e){
      cat(e$message)
      #message(e)
      cat("\nData was not analyzed.\n")
      NULL
    })
  }
  if(!is.null(out_data)){
    methodls$testMethod = out_data@typeC; methodls$scale = out_data@suppLs$scaleC; multiv_result$model_summary = out_data@modelDF;
    multiv_result$score_val = out_data@scoreMN; multiv_result$loading_val = out_data@loadingMN;
    multiv_result$oscore_val = out_data@orthoScoreMN; multiv_result$oloading_val = out_data@orthoLoadingMN;
    multiv_result$vip_val = out_data@vipVn; multiv_result$ovip_val = out_data@orthoVipVn; multiv_result$details = methodls;
    cat("\n",methodls$testMethod,"was performed.\n")
  }
  return(multiv_result)
}
