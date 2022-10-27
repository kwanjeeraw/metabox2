#'Impute missing values
#'@description Impute missing values in the input data.
#'@usage impute_missing_data(METBObj, method="mean", removeall=FALSE, cutoff=FALSE)
#'@param METBObj METBObj object contains list of data.
#'@param method name of missing value imputation method. Choose one from the list: zero, halfmin, min, mean, median, knn, bpca, ppca, svd, rf. Default is median.
#'@param removeall logical indicating all variables/features with missing values are removed (TRUE) or not removed (FALSE; Default).
#'@param cutoff number indicating percent missing value cutoff (from 1 to 99; Default = FALSE). If cutoff = FALSE, it will impute any missing values with a chosen method.
#'Or else variables with missing values > percent cutoff will be removed.
#'@return METBObj object.
#'@author Kwanjeera W \email{kwanjeera.wan@@mahidol.ac.th}
#'@references Wolfram Stacklies, et al. pcaMethods (2007) \url{https://doi.org/10.1093/bioinformatics/btm069}.
#'@references rfImpute \url{https://www.rdocumentation.org/packages/randomForest/versions/4.6-14/topics/rfImpute}.
#'@references Marietta Kokla, et al. Random forest-based imputation (2019) \url{https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-019-3110-0}.
#'@examples
#'#out = impute_missing_data(METBObj)
#'@export
impute_missing_data <- function(METBObj, method="mean", removeall=FALSE, cutoff=FALSE){
  #Check argument
  tmparg <- try(method <- match.arg(method, c("zero","halfmin","min","mean","median","knn","bpca","ppca","svd","rf"), several.ok = FALSE), silent = TRUE)
  if (class(tmparg) == "try-error") {
    cat("\nERROR! Argument 'method' is not valid, choose one from the list: zero,halfmin,min,mean,median,knn,bpca,ppca,svd,rf.\nMissing values were not imputed.\n")
    return(METBObj)
  }
  if(sum(is.na(METBObj$X)) == 0){#No missing value to impute
    cat("\nNo missing value to impute.\n")
    return(METBObj)
  }
  #Initialize parameters
  methodls = list(); #Working data
  cat("\nExecuting function ...")
  if(sum(is.na(METBObj$X))>0){#Missing values to impute
    dat = METBObj$X #working data
    if(removeall){#remove all variables/features with missing values
      incl = apply(is.na(dat), 2, sum) == 0
      new_dat = dat[,incl]
      removedls = colnames(dat)[!incl] #keep track of removed variables
      imputedls = NULL
      cat("\nRemove all variables with missing values.\n")
    }else{#impute variables/features with missing values
      if(cutoff < 0 | cutoff >= 100) {#Check argument
        cat("\nERROR! Argument 'cutoff' is not valid, the cutoff value is from 1 to 99.\nMissing values were not imputed.\n"); new_dat = METBObj$X
      }else{
        if(!cutoff) {#impute any missing values with chosen method
          incl_dat = METBObj$X
          removedls = NULL
          cat("\nImpute missing values with a chosen method.\n")
        }else{#remove variables/features with defined % missing values
          incl = apply(is.na(dat), 2, sum)/nrow(dat) <= cutoff/100
          incl_dat = dat[,incl]
          if(!is.data.frame(incl_dat)){
            incl_dat = as.data.frame(incl_dat)
            colnames(incl_dat) = colnames(dat)[incl]
          }
          removedls = colnames(dat)[!incl] #keep track of removed variables
          cat("\nRemove variables with >",cutoff,"% missing values.\n")
        }

        #Impute missing values with chosen method
        if(method=="zero"){#impute variables/features with min values
          imp = apply(is.na(incl_dat), 2, sum) > 0
          imputedls = colnames(incl_dat)[imp] #keep track of imputed variables
          new_dat = data.frame(apply(incl_dat, 2, function(x){
            if(sum(is.na(x)) > 0){
              x[is.na(x)] = 0}
            x
          }), check.names = FALSE)
          cat("\nImpute missing values with '0' value.\n")
        }
        if(method=="halfmin"){#impute variables/features with min values
          imp = apply(is.na(incl_dat), 2, sum) > 0
          imputedls = colnames(incl_dat)[imp] #keep track of imputed variables
          new_dat = data.frame(apply(incl_dat, 2, function(x){
            if(sum(is.na(x)) > 0){
              x[is.na(x)] = min(x, na.rm=TRUE)/2}
            x
          }), check.names = FALSE)
          cat("\nImpute missing values with 'half of min' value of each column.\n")
        }
        if(method=="min"){#impute variables/features with min values
          imp = apply(is.na(incl_dat), 2, sum) > 0
          imputedls = colnames(incl_dat)[imp] #keep track of imputed variables
          new_dat = data.frame(apply(incl_dat, 2, function(x){
            if(sum(is.na(x)) > 0){
              x[is.na(x)] = min(x, na.rm=TRUE)}
            x
          }), check.names = FALSE)
          cat("\nImpute missing values with 'min' value of each column.\n")
        }
        if (method=="mean"){#impute variables/features with mean values
          imp = apply(is.na(incl_dat), 2, sum) > 0
          imputedls = colnames(incl_dat)[imp] #keep track of imputed variables
          new_dat = data.frame(apply(incl_dat, 2, function(x){
            if(sum(is.na(x)) > 0){
              x[is.na(x)] = mean(x, na.rm=TRUE)}
            x
          }), check.names = FALSE)
          cat("\nImpute missing values with 'mean' of each column.\n")
        }
        if (method == "median"){#impute variables/features with median values
          imp = apply(is.na(incl_dat), 2, sum) > 0
          imputedls = colnames(incl_dat)[imp] #keep track of imputed variables
          new_dat = data.frame(apply(incl_dat, 2, function(x){
            if(sum(is.na(x)) > 0){
              x[is.na(x)] = median(x, na.rm=TRUE)}
            x
          }), check.names = FALSE)
          cat("\nImpute missing values with 'median' of each column.\n")
        }
        if(method == "knn"){#impute variables/features with k-nearest neighbours based on similar samples
          imp = apply(is.na(incl_dat), 2, sum) > 0
          imputedls = colnames(incl_dat)[imp] #keep track of imputed variables
          new_dat = data.frame(impute::impute.knn(data.matrix(incl_dat))$data, check.names = FALSE)
          cat("\nImpute missing values with 'knn' method.\n")
        }
        if(method == "bpca"){#impute variables/features with Bayesian PCA (BPCA)
          imp = apply(is.na(incl_dat), 2, sum) > 0
          imputedls = colnames(incl_dat)[imp] #keep track of imputed variables
          new_dat = data.frame(pcaMethods::pca(incl_dat, nPcs =5, method="bpca", center=T)@completeObs, check.names = FALSE)
          cat("\nImpute missing values with 'BPCA' method.\n")
        }
        if(method == "ppca"){#impute variables/features with probabilistic PCA (PPCA)
          imp = apply(is.na(incl_dat), 2, sum) > 0
          imputedls = colnames(incl_dat)[imp] #keep track of imputed variables
          new_dat = data.frame(pcaMethods::pca(incl_dat, nPcs =5, method="ppca", center=T)@completeObs, check.names = FALSE)
          cat("\nImpute missing values with 'PPCA' method.\n")
        }
        if(method == "svd"){#impute variables/features with Singular Value Decomposition (SVD) impute
          imp = apply(is.na(incl_dat), 2, sum) > 0
          imputedls = colnames(incl_dat)[imp] #keep track of imputed variables
          new_dat = data.frame(pcaMethods::pca(incl_dat, nPcs =5, method="svdImpute", center=T)@completeObs, check.names = FALSE)
          cat("\nImpute missing values with 'SVD' method.\n")
        }
        if(method == "rf"){#impute variables/features with randomForest
          imp = apply(is.na(incl_dat), 2, sum) > 0
          imputedls = colnames(incl_dat)[imp] #keep track of imputed variables
          new_dat = randomForest::rfImpute(incl_dat, METBObj$Y, iter=5, ntree=300)[,-1]
          cat("\nImpute missing values with 'RF' method.\n")
        }
      }
    }
    methodls$method = method; methodls$removeall = removeall; methodls$cutoff = cutoff; methodls$removed = removedls; methodls$imputed = imputedls;
    METBObj$X = new_dat; METBObj$details$impute_data = methodls; #update input data
    return(METBObj)
  }
}
