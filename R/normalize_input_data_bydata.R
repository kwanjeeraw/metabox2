#'Data-driven normalizations
#'@description Normalize samples by data-driven normalization methods.
#'@usage normalize_input_data_bydata(METBObj, method="median", ref=NULL)
#'@param METBObj METBObj object contains list of data.
#'@param method name of normalization method. Choose one from the list: contrast, cubic, eigenms, linear, loess, mean, median, pqn, quantile, sum, vsn. Default is median.
#'@param ref a row number/index or a vector of row numbers of reference/control sample(s) in the \code{METBObj$inputdata} data frame.
#'This parameter is required for pqn method.
#'@details The row number or the row index always begins with 1.
#'@return METBObj object.
#'@author Kwanjeera W \email{kwanjeera.wan@@mahidol.ac.th}
#'@references Astrand M. contrast (2003) \url{https://doi.org/10.1089/106652703763255697}.
#'@references Frank Dieterle, et al. pqn (2006) \url{https://pubs.acs.org/doi/10.1021/ac051632c}.
#'@references Bolstad Ben. preprocessCore.quantiles (2021) \url{https://github.com/bmbolstad/preprocessCore}.
#'@references Yuliya V. Karpievitch, et al. EigenMS (2014) \url{https://doi.org/10.1371/journal.pone.0116221}.
#'@references Bethanne M.Warrack, et al. MSTUS (2009) \url{https://doi.org/10.1016/j.jchromb.2009.01.007}.
#'@references Sandrine Dudoit, et al. cyclic loess (2002) \url{https://wwwf.imperial.ac.uk/~das01/BioinformaticsMSc/Papers/sinica.final.pdf}.
#'@references Christopher Workman, et al. cubic spline (2002) \url{https://genomebiology.biomedcentral.com/articles/10.1186/gb-2002-3-9-research0048}.
#'@references Chetwynd, AJ, et al. MSTUS (2016) \url{https://doi.org/10.1016/j.chroma.2015.12.056}.
#'@references Wolfgang Huber, et al. VSN (2002) \url{https://doi.org/10.1093/bioinformatics/18.suppl_1.S96}.
#'@references E. S. Motakis, et al. Variance stabilization and normalization (2006) \url{https://doi.org/10.1093/bioinformatics/btl412}.
#'@seealso \code{\link[preprocessCore:normalize.quantiles]{preprocessCore::normalize.quantiles()}}, \code{\link[affy:normalize.loess]{affy::normalize.loess()}}, \code{\link[affy:normalize.qspline]{affy::normalize.qspline()}}, \code{\link[affy:normalize.invariantset]{affy::normalize.invariantset()}}, \code{\link[vsn:vsn2]{vsn::vsn2()}}
#'@examples
#'#out = normalize_input_data_bydata(METBObj)
#'@export
normalize_input_data_bydata <- function(METBObj, method="median", ref=NULL){
  #Check argument
  tmparg <- try(method <- match.arg(method, c("contrast","cubic","eigenms","linear","loess","mean","median","pqn","quantile","sum","vsn"), several.ok = FALSE), silent = TRUE)
  if (class(tmparg) == "try-error") {
    cat("ERROR! Argument 'method' is not valid, choose one from the list: contrast,cubic,eigenms,linear,loess,mean,median,pqn,quantile,sum,vsn.\nData was not normalized.\n")
    return(METBObj)
  }
  if(sum(is.na(METBObj$X)) > 0){#Data contains missing values
    cat("WARNING! The data contains missing values. contrast, cubic methods will not normalize the data. The others will retain NAs in the normalized data.
        Imputing the missing values before normalization is recommended.\n")
    #return(METBObj)
  }
  #Initialize parameters
  dat = METBObj$X; methodls = list(); printtxt=""; #working data
  cat("Executing function ...\n")
  #Sample-based normalization method to make each sample comparable
  if (method == "mean"){#mean
    new_dat = data.frame(t(apply(dat, 1, function(x){
      x/mean(x, na.rm=TRUE)
    })), check.names=FALSE)
    printtxt = paste0(printtxt,"\nData normalization with 'mean'.\n")
  }
  if (method == "median"){#median
    new_dat = data.frame(t(apply(dat, 1, function(x){
      x/median(x, na.rm=TRUE)
    })), check.names=FALSE)
    printtxt = paste0(printtxt,"\nData normalization with 'median'.\n")
  }
  if (method == "sum"){#sum or MSTUS
    new_dat = data.frame(t(apply(dat, 1, function(x){
      1000*x/sum(x, na.rm=TRUE)
    })), check.names=FALSE)
    printtxt = paste0(printtxt,"\nData normalization with 'sum'.\n")
  }

  if (method == "contrast"){#contrast
    inst_pkg = NULL
    if(!requireNamespace("affy", quietly = TRUE)){#check and install required package
      cat("\nMissing the required package 'affy', trying to install the package ...\n")
      inst_pkg = install_pkgs('affy')
    }
    if(length(unlist(inst_pkg))){
      new_dat = dat
      printtxt = paste0(printtxt,"\nERROR! Could not install the required package 'affy'. Data was not normalized.\nReturning unprocessed data ...\n")
    }else{
      m_dat = data.frame(t(dat))
      threshold = 1e-11
      m_dat[m_dat <= 0] = threshold
      new_dat = tryCatch({
        contrast_data = affy::maffy.normalize(data.matrix(m_dat), subset = 1:nrow(m_dat), span = 0.75, verbose = TRUE, family = "gaussian", log.it = FALSE)
        subtract <- function(x){
          t(t(x)-apply(x,2,quantile,0.1,na.rm = TRUE))
        }
        printtxt = paste0(printtxt,"\nData normalization with 'contrast normalization'.\n")
        ct_dat = data.frame(t(subtract(contrast_data)), row.names = c(1:nrow(dat)))
        colnames(ct_dat) = colnames(dat)
        ct_dat
      },
      error=function(e){
        message(e)
        printtxt = paste0(printtxt,e$message,"\nERROR! Data was not normalized.\nReturning unprocessed data ...\n")
        dat
      })
    }
  }
  if (method == "cubic"){#cubic
    inst_pkg = NULL
    if(!requireNamespace("affy", quietly = TRUE)){#check and install required package
      cat("\nMissing the required package 'affy', trying to install the package ...\n")
      inst_pkg = install_pkgs('affy')
    }
    if(length(unlist(inst_pkg))){
      new_dat = dat
      printtxt = paste0(printtxt,"\nERROR! Could not install the required package 'affy'. Data was not normalized.\nReturning unprocessed data ...\n")
    }else if(ncol(dat) < 35){
      new_dat = dat
      printtxt = paste0(printtxt,"\nERROR! Number of variables/features must be more than 35 variables. Data was not normalized.\nReturning unprocessed data ...\n")
    }else{
      m_dat = data.frame(t(dat))
      new_dat = tryCatch({
        spline_data = affy::normalize.qspline(m_dat, samples=0.1, target=apply(m_dat,1,mean,na.rm = TRUE), verbose = FALSE)
        printtxt = paste0(printtxt,"\nData normalization with 'cubic spline normalization'.\n")
        cb_dat = data.frame(t(spline_data), row.names = c(1:nrow(dat)))
        colnames(cb_dat) = colnames(dat)
        cb_dat
      },
      error=function(e){
        message(e)
        printtxt = paste0(printtxt,e$message,"\nERROR! Data was not normalized.\nReturning unprocessed data ...\n")
        dat
      })
    }
  }
  if (method == "eigenms"){#eigenms
    grps = as.factor(METBObj$inputdata[, METBObj$classCol])
    m_info = cbind(seq(1:ncol(dat)), colnames(dat))
    m_dat = data.frame(t(dat))
    m_dat = log2(m_dat) #log abundances of the metabolites
    new_dat = tryCatch({
      printtxt = paste0(printtxt,"\nData normalization with 'EigenMS'.\n")
      m_ints_eig1 = eig_norm1(m=m_dat, treatment=grps, prot.info=m_info)
      m_ints_norm1 = eig_norm2(rv=m_ints_eig1)
      eg_dat = data.frame(t(m_ints_norm1$norm_m), row.names = c(1:nrow(dat)))
      colnames(eg_dat) = colnames(dat)
      2^(eg_dat) #anti-log
    },
    error=function(e){
      message(e)
      printtxt = paste0(printtxt,e$message,"\nERROR! Data was not normalized.\nReturning unprocessed data ...\n")
      dat
    })
  }
  if (method == "linear"){#linear
    m_dat = data.frame(t(dat))
    new_dat = tryCatch({
      printtxt = paste0(printtxt,"\nData normalization with 'linear baseline normalization'.\n")
      linear_baseline = apply(m_dat, 1, median, na.rm = TRUE)
      baseline_mean = mean(linear_baseline, na.rm = TRUE)
      sample_means = apply(m_dat, 2, mean, na.rm = TRUE)
      linear_scaling = baseline_mean/sample_means
      linear_data = t(t(m_dat) * linear_scaling)
      ln_dat = data.frame(t(linear_data), row.names = c(1:nrow(dat)))
      colnames(ln_dat) = colnames(dat)
      ln_dat
    },
    error=function(e){
      message(e)
      printtxt = paste0(printtxt,e$message,"\nERROR! Data was not normalized.\nReturning unprocessed data ...\n")
      dat
    })
  }
  # if (method == "liwong"){#liwong - unknown error with negative values
  #   inst_pkg = NULL
  #   if(!requireNamespace("affy", quietly = TRUE)){#check and install required package
  #     cat("\nMissing the required package 'affy', trying to install the package ...\n")
  #     inst_pkg = install_pkgs('affy')
  #   }
  #   if(length(unlist(inst_pkg))){
  #     new_dat = dat
  #     printtxt = paste0(printtxt,"\nERROR! Could not install the required package 'affy'. Data was not normalized.\nReturning unprocessed data ...\n")
  #   }else{
  #     new_dat = tryCatch({
  #       printtxt = paste0(printtxt,"\nNormalizing based on a baseline sample.\n")
  #       m_dat = data.frame(t(dat))
  #       average.intensity = apply(m_dat, 2, mean, na.rm = TRUE)
  #       median.number = round(ncol(m_dat)/2 + 0.1)
  #       ordering = order(average.intensity)
  #       median.sample.number = ordering[median.number]
  #       median.sample = m_dat[,median.sample.number]
  #       liwong_data <- vector()
  #       for(i in 1:ncol(m_dat)){
  #         tryCatch({
  #           liwong.model = affy::normalize.invariantset(data=m_dat[,i], ref=median.sample, prd.td=c(0.003,0.007))
  #           liwong.sample = predict(liwong.model$n.curve$fit, m_dat[,i])
  #         },error = function(er){
  #           liwong.sample = list()
  #           liwong.sample$y = m_dat[,i]
  #         })
  #         liwong_data <- cbind(liwong_data, liwong.sample$y)
  #       }
  #       data.frame(t(liwong_data), row.names = c(1:nrow(dat)))
  #     },
  #     error=function(e){
  #       message(e)
  #       printtxt = paste0(printtxt,e$message,"\nERROR! Data was not normalized.\nReturning unprocessed data ...\n")
  #       dat
  #     })
  #     colnames(new_dat) = colnames(dat)
  #   }
  # }
  if (method == "loess"){#loess
    inst_pkg = NULL
    if(!requireNamespace("affy", quietly = TRUE)){#check and install required package
      cat("\nMissing the required package 'affy', trying to install the package ...\n")
      inst_pkg = install_pkgs('affy')
    }
    if(length(unlist(inst_pkg))){
      new_dat = dat
      printtxt = paste0(printtxt,"\nERROR! Could not install the required package 'affy'. Data was not normalized.\nReturning unprocessed data ...\n")
    }else{
      m_dat = data.frame(t(dat))
      new_dat = tryCatch({
        loess_data = affy::normalize.loess(m_dat, subset = 1:nrow(m_dat), log.it = TRUE, verbose = TRUE, span = 0.75, family.loess = "gaussian")
        printtxt = paste0(printtxt,"\nData normalization with 'cyclic locally weighted regression'.\n")
        lo_dat = data.frame(t(loess_data), row.names = c(1:nrow(dat)))
        colnames(lo_dat) = colnames(dat)
        lo_dat
      },
      error=function(e){
        message(e)
        printtxt = paste0(printtxt,e$message,"\nERROR! Data was not normalized.\nReturning unprocessed data ...\n")
        dat
      })
    }
  }
  if (method == "pqn"){#pqn
    #Probabilistic quotient normalization (PQN)
    PQN_fn <- function(x, ref_samp){
      x/median(as.numeric(x/ref_samp), na.rm=TRUE)
    }
    new_dat = tryCatch({
      if(!is.null(ref)){#norm by a reference group, control samples recommended
        printtxt = paste0(printtxt,"\nProbabilistic quotient normalization (PQN) by the median of reference sample(s).\n")
        ref_samp = apply(dat[ref, , drop=FALSE], 2, median)
        data.frame(t(apply(dat, 1, PQN_fn, ref_samp)), check.names=FALSE)
      }else{#norm by all samples
        printtxt = paste0(printtxt,"\nProbabilistic quotient normalization (PQN) by the median of all samples.\n")
        ref_samp = apply(dat, 2, median)
        data.frame(t(apply(dat, 1, PQN_fn, ref_samp)), check.names=FALSE)
      }
    },
    error=function(e){
      message(e)
      printtxt = paste0(printtxt,e$message,"\nERROR! Data was not normalized.\nReturning unprocessed data ...\n")
      dat
    })
  }
  if (method == "quantile"){#quantile
    inst_pkg = NULL
    if(!requireNamespace("preprocessCore", quietly = TRUE)){#check and install required package
      cat("\nMissing the required package 'preprocessCore', trying to install the package ...\n")
      inst_pkg = install_pkgs('preprocessCore')
    }
    if(length(unlist(inst_pkg))){
      new_dat = dat
      printtxt = paste0(printtxt,"\nERROR! Could not install the required package 'preprocessCore'. Data was not normalized.\nReturning unprocessed data ...\n")
    }else{
      new_dat = tryCatch({
        printtxt = paste0(printtxt,"\nData normalization with 'quantile normalization'.\n")
        qn_data = data.frame(t(preprocessCore::normalize.quantiles(t(dat), copy=FALSE)), check.names=FALSE)
        qn_data
      },
      error=function(e){
        message(e)
        printtxt = paste0(printtxt,e$message,"\nERROR! Data was not normalized.\nReturning unprocessed data ...\n")
        dat
      })
    }
  }
  if (method == "vsn"){#vsn
    inst_pkg = NULL
    if(!requireNamespace("vsn", quietly = TRUE)){#check and install required package
      cat("\nMissing the required package 'vsn', trying to install the package ...\n")
      inst_pkg = install_pkgs('vsn')
    }
    if(length(unlist(inst_pkg))){
      new_dat = dat
      printtxt = paste0(printtxt,"\nERROR! Could not install the required package 'vsn'. Data was not normalized.\nReturning unprocessed data ...\n")
    }else{
      new_dat = tryCatch({
        printtxt = paste0(printtxt,"\nData normalization and transformation with 'vsn'.\n")
        vsn_model = vsn::vsn2(t(dat),minDataPointsPerStratum=5L)
        vsn_data = vsn::predict(vsn_model, t(dat))
        data.frame(t(vsn_data), check.names=FALSE)
      },
      error=function(e){
        message(e)
        printtxt = paste0(printtxt,e$message,"\nERROR! Data was not normalized.\nReturning unprocessed data ...\n")
        dat
      })
    }
  }
  printtxt = paste(printtxt,"\nData summary:\n*",
                   nrow(new_dat), "samples and", ncol(new_dat), "variables.\n**",
                   sum(apply(new_dat, 2, function(x) { sum(x<0) }),na.rm = TRUE)," negative variables.\n***",
                   sum(is.na(new_dat)), "(",round((sum(is.na(new_dat))/(nrow(new_dat)*ncol(new_dat)))*100,2),"% ) missing values.\n")
  cat(printtxt)
  methodls$method = method; methodls$ref = ref;
  METBObj$X = new_dat; METBObj$details$normalize_bydata = methodls; METBObj$text = printtxt;#update input data
  return(METBObj)
}
