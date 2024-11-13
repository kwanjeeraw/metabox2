#'Internal standards (IS)-based and quality control samples (QC)-based normalizations
#'@description Normalize samples by IS-based and QC-based methods.
#'@usage normalize_input_data_byqc(METBObj, method="nomis", istd=NULL, factorCol=NULL,
#'batch=FALSE, sampleType=FALSE, injectionOrder=FALSE)
#'@param METBObj METBObj object contains list of data.
#'@param method name of scaling method. Choose one from the list: ccmn, nomis, ruvrand, serrf, loess. Default is nomis.
#'@param istd a column number/index or a vector of column numbers of internal standard(s) (IS) in the \code{METBObj$X} data frame.
#'This parameter is not required for serrf and loess methods.
#'@param factorCol a column number/index of biological factors of interest in the \code{METBObj$inputdata} data frame.
#'This parameter is required for ccmn and ruv2 methods.
#'@param batch a column number/index of the batch information in the \code{METBObj$inputdata} data frame.
#'This parameter is required for serrf and loess methods. Each batch must contain at least 5 qc samples.
#'@param sampleType a column number/index indicating sample types in the \code{METBObj$inputdata} data frame.
#'This parameter is required for serrf and loess methods.
#'
#'Use "qc" to indicate pooled/quality control (QC) samples.
#'Use "sample" to indicate other samples.
#'sampleType allows NA, which will not be normalized by serrf or loess, e.g. blank samples.
#'sampleType allows other type of samples like "validate1", "validate2", e.g. validate QCs.
#'
#'@param injectionOrder a column number/index of the injection order in the \code{METBObj$inputdata} data frame.
#'This parameter is required for serrf and loess methods. The injection order must be numeric, and unique.
#'@details The column number or the column index always begins with 1.
#'@return METBObj object
#'@author Kwanjeera W \email{kwanjeera.wan@@mahidol.ac.th}
#'@references Henning Redestig, et al. CCMN (2009) \url{https://pubs.acs.org/doi/10.1021/ac901143w}.
#'@references Marko Sysi-Aho, et al. NOMIS (2007) \url{https://doi.org/10.1186/1471-2105-8-93}.
#'@references Alysha M. De Livera, et al. RUV (2012) \url{https://doi.org/10.1021/ac302748b}.
#'@references Alysha M. De Livera, et al. RUV (2015) \url{https://doi.org/10.1021/ac502439y}.
#'@references Ramyar Molania, et al. RUV-III (2019) \url{https://doi.org/10.1093/nar/gkz433}.
#'@references Sili Fan, et al. SERRF (2019) \url{https://doi.org/10.1021/acs.analchem.8b05592}.
#'@references LOESS.
#'@seealso \code{\link[crmn:normalize]{crmn::normalize()}}, \code{\link[MetNorm:NormalizeRUVRand]{MetNorm::NormalizeRUVRand()}}
#'@examples
#'#load('GitHub/metabox-ml/data/qc_batch_example.RData')
#'#out = normalize_input_data_byqc(qc_batch_example, method="serrf",
#'#batch=2, sampleType=3, injectionOrder=5)
#'@export
normalize_input_data_byqc <- function(METBObj, method="nomis", istd=NULL, factorCol=NULL, batch=FALSE, sampleType=FALSE, injectionOrder=FALSE){
  #Check argument
  tmparg <- try(method <- match.arg(method, c("ccmn","nomis","ruvrand","serrf","loess"), several.ok = FALSE), silent = TRUE)
  if (class(tmparg) == "try-error") {
    cat("ERROR! Argument 'method' is not valid, choose one from the list: ccmn,nomis,ruvrand,serrf,loess.\nData was not normalized.\n")
    return(METBObj)
  }
  if(sum(is.na(METBObj$X)) > 0){#Data contains missing values
    cat("WARNING! The data contains missing values, which will be retained in the normalized data.
        Imputing the missing values before normalization is recommended.\n")
    #return(METBObj)
  }

  #Initialize parameters
  dat = METBObj$X; methodls = list(); istdls=""; printtxt=""; #Working data
  cat("Executing function ...")
  #IS-based normalization methods
  if (method == "ccmn"){#ccmn
    if(!is.null(factorCol) && !is.null(istd)){
      inputdat = METBObj$inputdata #working data
      G = model.matrix(~-1+. , data=data.frame(inputdat[,factorCol])) #biological factors, G
      m_dat = t(dat)
      colnames(m_dat) = rownames(dat)
      isIS = logical(ncol(dat))
      isIS[istd] = TRUE #IS list
      istdls = colnames(dat)[isIS] = paste0("STD_",colnames(dat)[isIS]) #rename IS
      new_dat = tryCatch({
        ccmn_data = crmn::normalize(m_dat, "crmn", factors=G, standards=isIS, ncomp=2, lg=TRUE)
        printtxt = paste0(printtxt,"\nData normalization with 'ccmn'.\n")
        cbind(dat[isIS],data.frame(t(ccmn_data), check.names = FALSE)) #output includes unadjusted IS columns
      },
      error=function(e){
        message(e)
        printtxt = paste0(printtxt,e$message,"\nERROR! Data was not normalized.\nReturning unprocessed data ...\n")
        dat
      })
    }else{new_dat = dat; printtxt = paste0(printtxt,"\nERROR! Argument 'factorCol' and 'istd' are required for ccmn method.\nData was not normalized.\nReturning unprocessed data ...\n")}
  }
  if (method == "nomis"){#nomis
    if(!is.null(istd)){
      m_dat = t(dat)
      colnames(m_dat) = rownames(dat)
      isIS = logical(ncol(dat))
      isIS[istd] = TRUE #IS list
      istdls = colnames(dat)[isIS] = paste0("STD_",colnames(dat)[isIS]) #rename IS
      new_dat = tryCatch({
        nomis_data = crmn::normalize(m_dat, "nomis", standards=isIS, lg=TRUE)
        printtxt = paste0(printtxt,"\nData normalization with 'nomis'.\n")
        cbind(dat[isIS],data.frame(t(nomis_data), check.names = FALSE)) #output includes unadjusted IS columns
      },
      error=function(e){
        message(e)
        printtxt = paste0(printtxt,e$message,"\nERROR! Data was not normalized.\nReturning unprocessed data ...\n")
        dat
      })
    }else{new_dat = dat; printtxt = paste0(printtxt,"\nERROR! Argument 'istd' are required for nomis method.\nData was not normalized.\nReturning unprocessed data ...\n")}
  }
  if (method == "ruvrand"){#RUVRand
    cat("\nWARNING! RUVRand function is in log2 scale.\n")
    if(!is.null(istd)){
      m_dat = log2(data.matrix(dat, rownames.force = T)) #log abundances of the metabolites
      isIS = log2(dat[,istd])
      istdls = colnames(m_dat)[istd] = paste0("STD_",colnames(m_dat)[istd]) #rename IS
      if(length(istd)>1){
        printtxt = paste0(printtxt,"\nUse the average of the internal standards.\n")
        #if multiple internal standards are available,
        #simply find those which correlate highly with the internal standard or the average of the internal standards
        isIS = rowMeans(isIS)
      }
      r = numeric(dim(m_dat)[2])
      for(j in 1:length(r)){
        r[j] = cor(isIS, m_dat[,j])
      }
      ctl = logical(length(r))
      ctl[which(r>round(quantile(r,0.7,na.rm = TRUE),2))] = TRUE
      new_dat = tryCatch({
        ruv_data = MetNorm::NormalizeRUVRand(Y=m_dat, ctl=ctl, k=1, lambda=0.03, plotk = FALSE)
        printtxt = paste0(printtxt,"\nData normalization with 'RUV-random'.\n")
        data.frame((ruv_data$newY), check.names = FALSE) #output includes adjusted IS columns, it is in log-scale
      },
      error=function(e){
        message(e)
        printtxt = paste0(printtxt,e$message,"\nERROR! Data was not normalized.\nReturning unprocessed data ...\n")
        dat
      })
    }else{new_dat = dat; printtxt = paste0(printtxt,"\nERROR! Argument 'istd' are required for nomis method.\nData was not normalized.\nReturning unprocessed data ...\n")}
  }
  # if (method == "ruv2"){#RUV2
  #   if(!is.null(factorCol) && !is.null(istd)){
  #     inst_pkg = NULL
  #     if(!requireNamespace("ruv", quietly = TRUE)){#check and install required package
  #       cat("\nMissing the required package 'ruv', trying to install the package ...\n")
  #       inst_pkg = install_pkgs('ruv')
  #     }
  #     if(length(unlist(inst_pkg))){
  #       new_dat = dat
  #       cat("\nERROR! Could not install the required package 'ruv'. Data was not normalized.\n")
  #     }else{
  #       inputdat = METBObj$inputdata #working data
  #       G = model.matrix(~-1+. , data=data.frame(inputdat[,factorCol])) #biological factors, G, error if > 2 groups
  #       m_dat = log2(data.matrix(dat, rownames.force = T)) #log abundances of the metabolites
  #       isIS = log2(dat[,istd])
  #       istdls = colnames(m_dat)[istd] = paste0("STD_",colnames(m_dat)[istd]) #rename IS
  #       if(length(istd)>1){
  #         cat("\nUse the average of the internal standards.\n")
  #         #if multiple internal standards are available,
  #         #simply find those which correlate highly with the internal standard or the average of the internal standards
  #         isIS = rowMeans(isIS)
  #       }
  #       r = numeric(dim(m_dat)[2])
  #       for(j in 1:length(r)){
  #         r[j] = cor(isIS, m_dat[,j])
  #       }
  #       ctl = logical(length(r))
  #       ctl[which(r>round(quantile(r,0.7,na.rm = TRUE),2))] = TRUE
  #       new_dat = tryCatch({
  #         ruv_data = 2^(ruv::RUV2(Y=m_dat, X=G, ctl=ctl, k=1)) #inverse log
  #         cat("\nData normalization with 'RUV-2'.\n")
  #         data.frame(ruv_data, check.names = FALSE) #output includes adjusted IS columns
  #       },
  #       error=function(e){
  #         message(e)
  #         cat("\nERROR! Data was not normalized.\n")
  #         dat
  #       })
  #     }
  #   }else{new_dat = dat; cat("\nERROR! Argument 'factorCol' and 'istd' are required for ruv-2 method.\nData was not normalized.\n")}
  # }
  # if (method == "ruviii"){#RUV-III
  #   if(!is.null(factorCol) && !is.null(istd)){
  #     inst_pkg = NULL
  #     if(!requireNamespace("ruv", quietly = TRUE)){#check and install required package
  #       cat("\nMissing the required package 'ruv', trying to install the package ...\n")
  #       inst_pkg = install_pkgs('ruv')
  #     }
  #     if(length(unlist(inst_pkg))){
  #       new_dat = dat
  #       cat("\nERROR! Could not install the required package 'ruv'. Data was not normalized.\n")
  #     }else{
  #       inputdat = METBObj$inputdata #working data
  #       #e.g. use ruv::replicate.matrix(data.frame(sugar_dt$inputdata$subjectID, sugar_dt$inputdata$time)) #replicate structure, M
  #       M = ruv::replicate.matrix(data.frame(METBObj$orgID,METBObj$Y))
  #       m_dat = log2(data.matrix(dat, rownames.force = T)) #log abundances of the metabolites
  #       isIS = log2(dat[,istd])
  #       istdls = colnames(m_dat)[istd] = paste0("STD_",colnames(m_dat)[istd]) #rename IS
  #       if(length(istd)>1){
  #         cat("\nUse the average of the internal standards.\n")
  #         #if multiple internal standards are available,
  #         #simply find those which correlate highly with the internal standard or the average of the internal standards
  #         isIS = rowMeans(isIS)
  #       }
  #       r = numeric(dim(m_dat)[2])
  #       for(j in 1:length(r)){
  #         r[j] = cor(isIS, m_dat[,j])
  #       }
  #       ctl = logical(length(r))
  #       ctl[which(r>round(quantile(r,0.7,na.rm = TRUE),2))] = TRUE
  #       new_dat = tryCatch({
  #         ruv_data = 2^(ruv::RUVIII(Y=m_dat, M=M, ctl=ctl, k=1)) #inverse log
  #         cat("\nData normalization with 'RUV-III'.\n")
  #         data.frame(ruv_data, check.names = FALSE) #output includes adjusted IS columns
  #       },
  #       error=function(e){
  #         message(e)
  #         cat("\nERROR! Data was not normalized.\n")
  #         dat
  #       })
  #     }
  #   }else{new_dat = dat; cat("\nERROR! Argument 'factorCol' and 'istd' are required for ruviii method.\nData was not normalized.\n")}
  # }
  #QC-based normalization methods
  if (method == "serrf"){#SERRF
    istdls = FALSE
    if((sampleType) && (injectionOrder) && (batch)){
      #--convert input data for SERRF --
      m_dat = t(dat)
      colnames(m_dat) = rownames(dat) #matrix of metabolites, e and e_matrix
      met_info = data.table::data.table(data.frame(label=as.character(colnames(dat)), No=as.character(1:ncol(dat)))) #metabolite info, f
      inputdat = METBObj$inputdata #working data
      inputdat[,sampleType]=tolower(inputdat[,sampleType]) #use lower case
      colnames(inputdat)[sampleType] = 'sampleType' #sampleType column
      colnames(inputdat)[injectionOrder] = 'time' #time column
      colnames(inputdat)[batch] = 'batch' #batch column
      samp_info = data.table::data.table(data.frame(lapply(inputdat[,c("sampleType","time","batch")], as.character))) #sample info, p
      #--convert input data for SERRF --
      checkinp = check_serrf(samp_info)
      if(checkinp==""){
        serrf_data = run_serrf(p=samp_info, f=met_info, e=m_dat, e_matrix=m_dat) #call SERRF function
        printtxt = paste0(printtxt,"\nData normalization with 'SERRF'.\n")
        new_dat = data.frame(t(serrf_data), check.names = FALSE) #output will not exclude QC samples and IS columns
      }else{
        printtxt = paste0(checkinp,printtxt,"\nERROR! Data was not normalized.\nReturning unprocessed data ...\n")
        new_dat = dat
      }
    }else{new_dat = dat; printtxt = paste0(printtxt,"\nERROR! Argument 'sampleType', 'injectionOrder' and 'batch' are required for SERRF method.\nData was not normalized.\nReturning unprocessed data ...\n")}
  }
  if (method == "loess"){#LOESS
    cat("\nWARNING! LOESS function could not normalize data in many cases.\n")
    istdls = FALSE
    if((sampleType) && (injectionOrder) && (batch)){
      #--convert input data for LOESS --
      m_dat = t(dat)
      colnames(m_dat) = rownames(dat) #matrix of metabolites, e and e_matrix
      met_info = data.table::data.table(data.frame(label=as.character(colnames(dat)), No=as.character(1:ncol(dat)))) #metabolite info, f
      inputdat = METBObj$inputdata #working data
      inputdat[,sampleType]=tolower(inputdat[,sampleType]) #use lower case
      colnames(inputdat)[sampleType] = 'sampleType' #sampleType column
      colnames(inputdat)[injectionOrder] = 'time' #time column
      colnames(inputdat)[batch] = 'batch' #batch column
      samp_info = data.table::data.table(data.frame(lapply(inputdat[,c("sampleType","time","batch")], as.character))) #sample info, p
      #--convert input data for LOESS --
      checkinp = check_serrf(samp_info)
      if(checkinp==""){
        loess_data = run_loess(p=samp_info, f=met_info, e=m_dat, e_matrix=m_dat) #call LOESS function
        printtxt = paste0(printtxt,"\nData normalization with 'LOESS'.\n")
        new_dat = data.frame(t(loess_data), check.names = FALSE) #output will not exclude QC samples and IS columns
      }else{
        printtxt = paste0(checkinp,printtxt,"\nERROR! Data was not normalized.\nReturning unprocessed data ...\n")
        new_dat = dat
      }
    }else{new_dat = dat; printtxt = paste0(printtxt,"\nERROR! Argument 'sampleType', 'injectionOrder' and 'batch' are required for LOESS method.\nData was not normalized.\nReturning unprocessed data ...\n")}
  }
  printtxt = paste(printtxt,"\nData summary:\n*",
                   nrow(new_dat), "samples and", ncol(new_dat), "variables.\n**",
                   sum(apply(new_dat, 2, function(x) { sum(x<0) }),na.rm = TRUE)," negative variables.\n***",
                   sum(is.na(new_dat)), "(",round((sum(is.na(new_dat))/(nrow(new_dat)*ncol(new_dat)))*100,2),"% ) missing values.\n")
  cat(printtxt)
  methodls$method = method; methodls$factorCol = factorCol; methodls$istd = istdls; methodls$sampleType=sampleType;
  methodls$injectionOrder = injectionOrder; methodls$batch = batch;
  METBObj$X = new_dat; METBObj$details$normalize_byqc = methodls; METBObj$text = printtxt;#update input data
  return(METBObj)
}
