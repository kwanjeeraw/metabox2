#'Correlation analysis
#'@description Perform correlation analysis.
#'@usage correlation_analyze(METBObj, y=NULL, method="pearson")
#'@param METBObj METBObj object contains list of data.
#'@param y a numeric vector or matrix or data frame of other variables e.g. genes, proteins or compounds, with samples/subjects as row.names. See details.
#'@param method name of multivariate analysis method. Choose one from the list: pearson, spearman. Default is pearson.
#'@details If y is provided, it will be combined to \code{METBObj$X} and correlation will be calculated for all \code{METBObj$X} and y variables.
#'@return a list of the following components:
#'
#'p_value = a matrix of p-values.
#'
#'p_adj = a matrix of adjusted p-values.
#'
#'corr_coef = a matrix of correlation coefficient.
#'
#'corr_data = a data frame of analysis result containing coeff,p_value,p_adj.
#'
#'details = a list of analysis details: testMethod, pAdjusted.
#'
#'@author Kwanjeera W \email{kwanjeera.wan@@mahidol.ac.th}
#'@seealso \code{\link[Hmisc:rcorr]{Hmisc::rcorr()}}, \code{\link[stats:p.adjust]{stats::p.adjust()}}
#'@examples
#'#sugar_dt = set_input_obj(sugar, 1,2,5)
#'#out=correlation_analyze(sugar_dt)
#'@export
correlation_analyze <- function(METBObj, y=NULL, method="pearson"){
  #Check argument
  tmparg <- try(method <- match.arg(tolower(method), c("pearson","spearman"), several.ok = FALSE), silent = TRUE)
  if (class(tmparg) == "try-error") {
    cat("\nERROR! Argument 'method' is not valid, choose one from the list: pearson,spearman.\nData was not analyzed.\n")
    return(FALSE)
  }
  #Check argument
  if(sum(is.na(METBObj$X)) > 0){#Data contains missing values
    cat("\nWarning! The data contains missing values.\n")
  }

  #Initialize parameters
  correlation_result = list(); methodls = list(); #Working data
  dat = as.matrix(METBObj$X); #Working data
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
  #Correltation analysis
  cat("\n\nPerforming correltation analysis ...\n")
  if(is.null(y)){
    out_data = tryCatch({
      Hmisc::rcorr(x=dat, y=y, type=method)
    },
    error=function(e){
      cat(e$message)
      #message(e)
      cat("\nData was not analyzed.\n")
      NULL
    })

  }else{
    if(!is.matrix(y)){
      y = as.matrix(y)
      cat("\nConverting a data y to a numeric matrix.\n")
    }
    out_data = tryCatch({
      Hmisc::rcorr(x=dat, y=y, type=method)
    },
    error=function(e){
      cat(e$message)
      #message(e)
      cat("\nData was not analyzed.\n")
      NULL
    })
  }
  if(!is.null(out_data)){
    corr_coef = out_data$r; p_value = out_data$P
    #do p-adjustment
    p_adj = out_data$P
    plw = p_adj[lower.tri(p_adj)]
    adj.plw = p.adjust(plw, method = "fdr")
    p_adj[lower.tri(p_adj)] = adj.plw
    p_adj[upper.tri(p_adj)] = 0
    p_adj = p_adj + t(p_adj)
    correlation_result$p_value = p_value; correlation_result$p_adj = p_adj; correlation_result$corr_coef = corr_coef;
    #melt matrix
    out_data$P[upper.tri(out_data$P)] = NA
    pv_tb = reshape2::melt(out_data$P)
    pv_tb = na.omit(pv_tb)
    p_adj[upper.tri(p_adj)] = NA
    pa_tb = reshape2::melt(p_adj)
    pa_tb = na.omit(pa_tb)
    out_data$r[upper.tri(out_data$r)] = NA
    r_tb = reshape2::melt(out_data$r)
    r_tb = na.omit(r_tb)
    p_tb = merge(pv_tb, pa_tb, by=c("Var1","Var2"), sort = F)
    #create dataframe output
    corr_data = merge(r_tb , p_tb, by=c("Var1","Var2"), sort = F)
    colnames(corr_data) = c("Var1","Var2","coeff","p_value","p_adj")
    methodls$testMethod = method; methodls$pAdjusted = "Benjamini & Hochberg or FDR";
    correlation_result$corr_data = corr_data; correlation_result$details = methodls;
    cat("\n",methodls$testMethod,"was performed.\n")
  }
  return(correlation_result)
}
