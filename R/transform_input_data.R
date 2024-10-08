#'Transform data
#'@description Transform variables/features to be closer to normal or Gaussian distribution. The generalized logarithm transformation is used to stabilize variance.
#'@usage transform_input_data(METBObj, method="sqrt")
#'@param METBObj METBObj object contains list of data.
#'@param method name of scaling method. Choose one from the list: cube, glog10, glog2, log10, log2, sqrt. Default is sqrt.
#'@return METBObj object.
#'@author Kwanjeera W \email{kwanjeera.wan@@mahidol.ac.th}
#'@references David M. Rocke, Blythe Durbin. glog (2003) \url{https://doi.org/10.1093/bioinformatics/btg107}.
#'@references Wolfgang Huber, et al. VSN (2002) \url{https://doi.org/10.1093/bioinformatics/18.suppl_1.S96}.
#'@references E. S. Motakis, et al. Variance stabilization and normalization (2006) \url{https://doi.org/10.1093/bioinformatics/btl412}.
#'@examples
#'#out = transform_input_data(METBObj)
#'@export
transform_input_data <- function(METBObj, method="sqrt"){
  #Check argument
  tmparg <- try(method <- match.arg(method, c("cube","glog10","glog2","log10","log2","sqrt"), several.ok = FALSE), silent = TRUE)
  if (class(tmparg) == "try-error") {
    cat("ERROR! Argument 'method' is not valid, choose one from the list: cube,glog10,glog2,log10,log2,sqrt.\nData was not transformed.\n")
    return(METBObj)
  }
  if(sum(is.na(METBObj$X)) > 0){#Data contains missing values
    cat("The data contains missing values. glog10 and glog2 methods will not transform the data. The others will retain NAs in the transformed data.
        Imputing the missing values before data transformation is recommended.\n")
    #return(METBObj)
  }
  #Initialize parameters
  dat = METBObj$X; methodls = list(); printtxt=""; #working data
  cat("Executing function ...")
  #Transform variables so that they are a more normal or Gaussian distribution as opposed to a skewed distribution
  if (method == "cube"){#cube root
    new_dat = sign(dat) * abs(dat)^{1/3}
    printtxt = paste0(printtxt,"\nData transformation with 'cube'.\n")
  }
  if (method == "glog10"){#glog10
    min.val = min(abs(dat[dat!=0]))/10
    new_dat = log10((dat + sqrt(dat^2 + min.val^2))/2)
    printtxt = paste0(printtxt,"\nData transformation with 'glog10'.\n")
  }
  if (method == "glog2"){#glog2^2
    min.val = min(abs(dat[dat!=0]))/2 #dat[datX<=0] = 1
    new_dat = log2((dat + sqrt(dat^2 + min.val^2))/2)
    printtxt = paste0(printtxt,"\nData transformation with 'glog2'.\n")
  }
  if (method == "log10"){#log10
    dat[dat==0] = 1
    new_dat = log10(dat)
    printtxt = paste0(printtxt,"\nData transformation with 'log10'.\n")
  }
  if (method == "log2"){#log2
    dat[dat==0] = 1
    new_dat = log2(dat)
    printtxt = paste0(printtxt,"\nData transformation with 'log2'.\n")
  }
  if (method == "sqrt"){#square root
    dat[dat==0] = 0
    new_dat = sqrt(dat)
    printtxt = paste0(printtxt,"\nData transformation with 'sqrt'.\n")
  }
  # if (method == "log2_2"){#glog2
  #   min.val = min(abs(dat[dat!=0]))/2 #dat[datX<=0] = 1
  #   new_dat = log2((dat + sqrt(dat^2 + min.val))/2)
  #   cat("\nData transformation with 'log2 transformation'.\n")
  # }
  printtxt = paste(printtxt,"\nData summary:\n*",
                   nrow(new_dat), "samples and", ncol(new_dat), "variables.\n**",
                   sum(apply(new_dat, 2, function(x) { sum(x<0) }),na.rm = TRUE)," negative variables.\n***",
                   sum(is.na(new_dat)), "(",round((sum(is.na(new_dat))/(nrow(new_dat)*ncol(new_dat)))*100,2),"% ) missing values.\n")
  cat(printtxt)
  methodls$method = method;
  METBObj$X = new_dat; METBObj$details$transform_data = methodls; METBObj$text = printtxt;#update input data
  return(METBObj)
}
