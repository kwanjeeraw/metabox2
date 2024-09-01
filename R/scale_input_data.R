#'Scale data
#'@description Scale data for between-variable/feature variations
#'@usage scale_input_data(METBObj, method="pareto")
#'@param METBObj METBObj object contains list of data.
#'@param method name of scaling method. Choose one from the list: auto, level, pareto, power, range, vast. Default is pareto.
#'@return METBObj object
#'@author Kwanjeera W \email{kwanjeera.wan@@mahidol.ac.th}
#'@references Qingxia Yang, et al. NOREVA (2020) \url{https://doi.org/10.1093/bib/bbz137}.
#'@examples
#'#out = scale_input_data(METBObj)
#'@export
scale_input_data <- function(METBObj, method="pareto"){
  #Check argument
  tmparg <- try(method <- match.arg(method, c("auto","level","pareto","power","range","vast"), several.ok = FALSE), silent = TRUE)
  if (class(tmparg) == "try-error") {
    cat("ERROR! Argument 'method' is not valid, choose one from the list: auto,level,pareto,power,range,vast.\nData was not scaled.\n")
    return(METBObj)
  }
  if(sum(is.na(METBObj$X)) > 0){#Data contains missing values
    cat("The data contains missing values, which will be retained in the scaled data.
        Imputing the missing values before scaling is recommended.\n")
    #return(METBObj)
  }
  #Initialize parameters
  dat = METBObj$X; methodls = list(); printtxt=""; #working data
  cat("Executing function ...")
  #Metabolite-based normalization method to adjust metabolite variances
  if (method == "auto"){#auto scaling or unit variance scaling
    new_dat = data.frame(scale(dat), check.names = FALSE)
    printtxt = paste0(printtxt,"\nData scaling with 'auto scaling'.\n")
  }
  if(method == "level") {#level scaling
    new_dat = data.frame(apply(dat, 2, function(x){
      (x - mean(x, na.rm=TRUE))/mean(x, na.rm=TRUE)
    }), check.names=FALSE)
    printtxt = paste0(printtxt,"\nData scaling with 'level scaling'.\n")
  }
  if(method == "pareto") {#pareto scaling
    new_dat = data.frame(apply(dat, 2, function(x){
      (x - mean(x, na.rm=TRUE))/sqrt(sd(x, na.rm=TRUE))
    }), check.names=FALSE)
    printtxt = paste0(printtxt,"\nData scaling with 'pareto scaling'.\n")
  }
  if(method == "power") {#power scaling
    new_dat = data.frame(apply(dat, 2, function(x){
      sqrt(x) - mean(sqrt(x))
    }), check.names=FALSE)
    printtxt = paste0(printtxt,"\nData scaling with 'power scaling'.\n")
  }
  if(method == "range") {#range scaling
    new_dat = data.frame(apply(dat, 2, function(x){
      (x - mean(x, na.rm=TRUE))/(max(x)-min(x))
    }), check.names=FALSE)
    printtxt = paste0(printtxt,"\nData scaling with 'range scaling'.\n")
  }
  if(method == "vast") {#vast scaling
    new_dat = data.frame(apply(dat, 2, function(x){
      mean(x, na.rm=TRUE) *
        (x - mean(x, na.rm=TRUE)) / (sd(x, na.rm=TRUE)**2)
    }), check.names=FALSE)
    printtxt = paste0(printtxt,"\nData scaling with 'vast scaling'.\n")
  }
  printtxt = paste(printtxt,"\nData summary:\n*",
                   nrow(new_dat), "samples and", ncol(new_dat), "variables.\n**",
                   sum(apply(new_dat, 2, function(x) { sum(x<0) }),na.rm = TRUE)," negative variables.\n***",
                   sum(is.na(new_dat)), "(",round((sum(is.na(new_dat))/(nrow(new_dat)*ncol(new_dat)))*100,2),"% ) missing values.\n")
  cat(printtxt)
  methodls$method = method;
  METBObj$X = new_dat; METBObj$details$scale_data = methodls; METBObj$text = printtxt;#update input data
  return(METBObj)
}
