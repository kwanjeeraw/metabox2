#'LOESS normalization
#'@description Normalize samples by QC-based methods, LOESS
#'@usage run_loess(p, f, e, e_matrix)
#'@param p is phenotypes. p is a data.frame containing following columns:
#'sampleType must include 'qc', 'sample'; It can include NA, which will not be normalized by LOESS, e.g. blank samples; It can include other type of samples like 'validate1', 'validate2', e.g. validate QCs.
#'time is the injection order of each sample. Must be numeric, and unique.
#'batch is batch information of samples. Each batch must contain at least 3 'qc's, although 5 qcs are recommended.
#'@param f is phenotypes.
#'@param e is the data matrix. Each row is a compound and each column is a sample.
#'@param e_matrix is the \code{e} data matrix. Each row is a compound and each column is a sample.
#'@return matrix with samples in columns and compounds in rows.
#'@author Sili Fan \email{}
#'@references Sili Fan, et al. SERRF (2019) \url{https://doi.org/10.1021/acs.analchem.8b05592}.
#'@examples
#'#loess_dat = run_loess(p, f, e, e_matrix)
#'@export
run_loess = function(p, f, e, e_matrix){

  remove_outlier = function(v){
    out = boxplot.stats(v)$out
    return(list(value = v[!v%in%out],index = which(v%in%out)))
  }
  loess_wrapper_extrapolate <- function (x, y, span.vals = seq(0.25, 1, by = 0.05), folds = 5){
    # Do model selection using mean absolute error, which is more robust than squared error.
    mean.abs.error <- numeric(length(span.vals))

    # Quantify error for each span, using CV
    loess.model <- function(x, y, span){
      loess(y ~ x, span = span, control=loess.control(surface="interpolate",statistics='exact'),family = "gaussian")
    }

    loess.predict <- function(fit, newdata) {
      predict(fit, newdata = newdata)
    }

    span.index <- 0

    for (each.span in span.vals) {
      span.index <- span.index + 1
      mean.abs.error[span.index] = tryCatch({
        y.hat.cv <- bootstrap::crossval(x, y, theta.fit = loess.model, theta.predict = loess.predict, span = each.span, ngroup = folds)$cv.fit
        non.empty.indices <- !is.na(y.hat.cv)
        diff = (y[non.empty.indices] / y.hat.cv[non.empty.indices]) * mean(y[non.empty.indices])
        sd(diff)/mean(diff)
      },error = function(er){
        NA
      })
    }
    best.span <- span.vals[which.min(mean.abs.error)]
    if(length(best.span)==0){
      best.span = 0.75
    }

    best.model <- loess(y ~ x, span = best.span, control=loess.control(surface="interpolate",statistics='exact'),family = "gaussian")

    return(list(best.model, min(mean.abs.error, na.rm = TRUE),best.span))
  }
  shiftData = function(ori,norm){
    ori.min = apply(ori,1,min,na.rm=T)
    norm.min = apply(norm,1,min,na.rm=T)
    return(norm - c(norm.min - ori.min))
  }
  # check if 'qc' or 'sample' is not in sampleType
  if(any(!c('qc','sample') %in% p$sampleType)){
    stop("The 'sampleType' must contain at least 'qc' and 'sample'. Please see example data for more information.")
  }
  # check if time is in the datasheet
  if(!"time" %in% colnames(p)){
    stop("Your data must have 'time'. Please see example data for more information.")
  }
  # check if time has duplicated value.
  if(any(duplicated(p$time))){
    time_duplicated = duplicated(p$time)
    stop("Your dataset has ",sum(time_duplicated)," duplicated 'time' value. 'time' of each sample should be unique. The duplicated 'time' values are ", paste0(p$time[time_duplicated],collapse = ', '),'.')
  }
  # check if batch is in the datasheet
  if(!"batch" %in% colnames(p)){
    stop("Your data must have 'batch'. Please see example data for more information.")
  }
  # check if any batch has too little qc
  if(any(table(p$batch, p$sampleType)[,'qc']<3)){
    stop("Some batches has a small number of QC that is not enough for training model. Each batch should have at least 3 QCs.")
  }

  # check with missing value. All missing values are by default replaced by half minimum.
  num_miss = sum(is.na(e))
  if(num_miss>0){
    cat(paste0("Your dataset has ",num_miss," missing values (i.e. empty cell). These missing values will be replaced by half of the minimum of non-missing value for each compound.\n"))
    for(i in 1:nrow(e)){
      e[i, is.na(e[i,])] = 0.5 * min(e[i, !is.na(e[i,])])
    }
  }
  # check with zero values.
  num_zero = sum(e == 0)
  if(num_zero>0){
    cat(paste0("Your dataset has ",num_zero," zeros. These zeros will be kept zeros in the final normalized data. \n"))
  }
  p$sample_index = paste0('p',1:nrow(p))
  original_colnames_e = colnames(e)
  colnames(e) = p$sample_index

  # Empty sampleType will not be used for normalization, but will be put back in the final sheet.
  if(any(is.na(p$sampleType))){
    cat(paste0("There are ",sum(is.na(p$sampleType))," empty cells in the 'sampleType'. They will not be used for normalization, but will be put back in the final sheet. \n"))
  }

  p_empty_sampleType = p[is.na(p$sampleType), ]
  e_empty_sampleType = e[,is.na(p$sampleType)]
  if('numeric' %in% class(e_empty_sampleType)){
    e_empty_sampleType = matrix(e_empty_sampleType, ncol = 1)
  }
  colnames(e_empty_sampleType) = p_empty_sampleType$sample_index
  e = e[, !is.na(p$sampleType)]
  p = p[!is.na(p$sampleType), ]

  # check if there is any infinite value.
  infinite_index = which(apply(e, 1, function(x){
    sum(is.infinite(x)) == length(x)
  }))
  if(!length(infinite_index)==0){
    e_infinite = e[infinite_index,]
    f_infinite = f[infinite_index,]

    e = e[-infinite_index,]
    f = f[-infinite_index,]
  }
  # split e, and p to different sample type.
  e_qc = e[, p$sampleType == 'qc']
  e_sample = e[, p$sampleType == 'sample']

  p_qc = p[p$sampleType == 'qc',]
  p_sample = p[p$sampleType == 'sample',]

  e_validates = list()
  p_validates = list()
  with_validate = any(!p$sampleType %in% c('qc','sample'))
  if(with_validate){
    val_RSDs = list()
    validate_types = unique(p$sampleType[!p$sampleType %in% c('qc','sample')])
    for(validate_type in validate_types){
      e_validates[[validate_type]] = e[, p$sampleType %in% validate_type]
      p_validates[[validate_type]] = p[p$sampleType %in% validate_type, ]
      val_RSDs[[validate_type]] = list()
    }
  }else{
    validate_types = NULL
  }
  aggregate_e = function(e_qc,e_sample,e_validates){
    e = do.call('cbind',c(list(e_qc, e_sample), e_validates))
    e = e[,order(as.numeric(gsub("p","",colnames(e))))]
    return(e)
  }

  #################### up to here. Identical with serrf as the above checks the input dataset format.

  e_norm = matrix(,nrow=nrow(e),ncol=ncol(e))
  QC.index = p[["sampleType"]]
  batch = p[["batch"]]
  time = as.numeric(p[["time"]])

  loess_fun_cv = function(e,train.index = QC.index,test.index=NULL,batch,time){
    cl = parallel::makeCluster(1)
    e_norm = parallel::parSapply(cl, X=1:nrow(e), function(i,e,train.index,batch,time,remove_outlier,loess_wrapper_extrapolate){

      e_norm = tryCatch({
        # for(i in 1:nrow(e)){
        line = e[i,]
        for(b in 1:length(unique(batch))){
          outlier_remove = remove_outlier(e[i,(batch %in% unique(batch)[b]) & (train.index=='qc')])
          if(length(outlier_remove$index) == 0){
            lm = loess_wrapper_extrapolate(x=time[(batch %in% unique(batch)[b]) & (train.index=='qc')], y = e[i,(batch %in% unique(batch)[b]) & (train.index=='qc')])[[1]]
          }else{
            lm = loess_wrapper_extrapolate(x=time[(batch %in% unique(batch)[b]) & (train.index=='qc')][-outlier_remove$index], y = e[i,(batch %in% unique(batch)[b]) & (train.index=='qc')][-outlier_remove$index])[[1]]
          }
          line[batch %in% unique(batch)[b]] = predict(lm,newdata  = time[batch %in% unique(batch)[b]])


          if(length(which(is.na(line[batch %in% unique(batch)[b]])))>0){
            for(j in which(is.na(line[batch %in% unique(batch)[b]]))){
              time_notNA = time[batch %in% unique(batch)[b]][-which(is.na(line[batch %in% unique(batch)[b]]))]
              closest_time = time_notNA[which.min(abs(time_notNA - time[batch %in% unique(batch)[b]][j]))]
              line[batch %in% unique(batch)[b]][j] = line[batch %in% unique(batch)[b]][which(time[batch %in% unique(batch)[b]]==closest_time)]
            }
          }

        }

        if(sum(line<0)>(length(line)/5)){
          stop("too many negative value. LOESS failed.")
        }else{
          line[line<0] = runif(sum(line<0), min = max(c(median(e[i,], na.rm = TRUE) -  0.1 * sd(e[i,], na.rm = TRUE),0)), max = max(c(median(e[i,], na.rm = TRUE) +  0.1 * sd(e[i,], na.rm = TRUE),1)))
        }

        # if(length(which(is.na(line)))>0){
        #   for(j in which(is.na(line))){
        #     time_notNA = time[-which(is.na(line))]
        #     closest_time = time_notNA[which.min(abs(time_notNA - time[j]))]
        #     line[j] = line[which(time==closest_time)]
        #   }
        # }
        e[i,] / (line / median(e[i,], na.rm = TRUE))
        # if(sum((e[i,] / (line / median(e[i,], na.rm = TRUE)))<0)){
        #   stop(i)
        # }
        # }


      },error = function(er){
        e[i,]
      })
      return(e_norm)
    },e,train.index,batch,time,remove_outlier,loess_wrapper_extrapolate)
    parallel::stopCluster(cl)
    e_norm = t(e_norm)
    if(!is.null(test.index)){
      rsd = RSD(e_norm[,test.index=='qc'])
    }else(
      rsd = NULL
    )
    return(list(data = e_norm,rsd = rsd))
  }

  cl = parallel::makeCluster(1)
  e_norm = loess_fun_cv(e,train.index = QC.index,test.index=NULL,batch=batch,time=time)[[1]]
  rownames(e_norm) = rownames(e)
  colnames(e_norm) = colnames(e)
  normalized_dataset = e_norm


  comb_p = rbind(p_empty_sampleType, p)
  comb_e = cbind(e_empty_sampleType, e)
  comb_p$sample_index = as.numeric(substring(comb_p$sample_index,2))
  comb_e = comb_e[,order(comb_p$sample_index)]
  if(!length(infinite_index)==0){
    e_full = rbind(e_infinite, normalized_dataset)
    f_full = rbind(f_infinite, f)
    e_full[infinite_index,] = NA
    e_full[-infinite_index,] = normalized_dataset
    normalized_dataset = e_full
  }
  normalized_dataset = cbind(e_empty_sampleType, normalized_dataset)
  normalized_dataset = normalized_dataset[,order(comb_p$sample_index)]
  comb_p = comb_p[order(comb_p$sample_index),]
  colnames(normalized_dataset) = comb_p$label
  normalized_dataset[e_matrix==0] = 0
  normalized_dataset[is.na(e_matrix)] = NA

  return(normalized_dataset) # the returned dataset: row compound; column samples.
}
