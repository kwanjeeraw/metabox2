#'SERRF normalization
#'@description Normalize samples by QC-based methods, SERRF.
#'@usage run_serrf(p, f, e, e_matrix)
#'@param p is phenotypes. p is a data.frame containing following columns:
#'sampleType must include 'qc', 'sample'; It can include NA, which will not be normalized by SERRF, e.g. blank samples; It can include other type of samples like 'validate1', 'validate2', e.g. validate QCs.
#'time is the injection order of each sample. Must be numeric, and unique.
#'batch is batch information of samples. Each batch must contain at least 2 'qc's, although 5 qcs are recommended.
#'@param f is phenotypes.
#'@param e is the data matrix. Each row is a compound and each column is a sample.
#'@param e_matrix is the \code{e} data matrix. Each row is a compound and each column is a sample.
#'@return matrix with samples in columns and compounds in rows.
#'@author Sili Fan \email{}
#'@references Sili Fan, et al. SERRF (2019) \url{https://doi.org/10.1021/acs.analchem.8b05592}.
#'@examples
#'#serrf_dat = run_serrf(p, f, e, e_matrix)
#'@export
run_serrf = function(p, f, e, e_matrix){

  pacman::p_load(ranger) # ranger is required package.

  shiftData = function(ori,norm){
    ori.min = apply(ori,1,min,na.rm=T)
    norm.min = apply(norm,1,min,na.rm=T)
    return(norm - c(norm.min - ori.min))
  }
  # check if 'qc' or 'sample' is not in sampleType
  if(any(!c('qc','sample') %in% p$sampleType)){
    stop("The 'sampleType' must contain at least 'qc' and 'sample'. Please see example data for more information. Data was not normalized.")
  }
  # check if time is in the datasheet
  if(!"time" %in% colnames(p)){
    stop("Your data must have 'time'. Please see example data for more information. Data was not normalized.")
  }
  # check if time has duplicated value.
  if(any(duplicated(p$time))){
    time_duplicated = duplicated(p$time)
    stop("Your dataset has ",sum(time_duplicated)," duplicated 'time' value. 'time' of each sample should be unique. The duplicated 'time' values are ", paste0(p$time[time_duplicated],collapse = ', '),'. Data was not normalized.')
  }
  # check if batch is in the datasheet
  if(!"batch" %in% colnames(p)){
    stop("Your data must have 'batch'. Please see example data for more information. Data was not normalized.")
  }
  # check if any batch has too little qc
  if(any(table(p$batch, p$sampleType)[,'qc']<2)){
    stop("Some batches has a small number of QC that is not enough for training model. Each batch should have at least 2 QCs. Data was not normalized.")
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
  cl = parallel::makeCluster(1)
  num = 10
  minus = FALSE
  serrfR = function(train = e[,p$sampleType == 'qc'],
                    target = e[,p$sampleType == 'sample'],
                    num = 10,
                    batch. = factor(c(batch[p$sampleType=='qc'],batch[p$sampleType=='sample'])),
                    time. = c(time[p$sampleType=='qc'],time[p$sampleType=='sample']),
                    sampleType. = c(p$sampleType[p$sampleType=='qc'],p$sampleType[p$sampleType=='sample']),minus=minus,cl){


    all = cbind(train, target)
    normalized = rep(0, ncol(all))

    for(j in 1:nrow(all)){
      for(b in 1:length(unique(batch.))){
        current_batch = levels(batch.)[b]
        # fixing zeros.
        all[j,batch.%in%current_batch][all[j,batch.%in%current_batch] == 0] = rnorm(length(all[j,batch.%in%current_batch][all[j,batch.%in%current_batch] == 0]),mean = min(all[j,batch.%in%current_batch][!is.na(all[j,batch.%in%current_batch])])+1,sd = 0.1*(min(all[j,batch.%in%current_batch][!is.na(all[j,batch.%in%current_batch])])+.1))
        # fixing NAs.
        all[j,batch.%in%current_batch][is.na(all[j,batch.%in%current_batch])] = rnorm(length(all[j,batch.%in%current_batch][is.na(all[j,batch.%in%current_batch])]),mean = 0.5*min(all[j,batch.%in%current_batch][!is.na(all[j,batch.%in%current_batch])])+1,sd = 0.1*(min(all[j,batch.%in%current_batch][!is.na(all[j,batch.%in%current_batch])])+.1))
      }
    }

    corrs_train = list()
    corrs_target = list()
    for(b in 1:length(unique(batch.))){

      current_batch = levels(batch.)[b]

      train_scale = t(apply(train[,batch.[sampleType.=='qc']%in%current_batch],1,scale))
      if(is.null(dim(target[,batch.[!sampleType.=='qc']%in%current_batch]))){ # !!!
        target_scale = scale(target[,batch.[!sampleType.=='qc']%in%current_batch])#!!!
      }else{
        # target_scale = scale(target[,batch.[!sampleType.=='qc']%in%current_batch])#!!!
        target_scale = t(apply(target[,batch.[!sampleType.=='qc']%in%current_batch],1,scale))
      }

      # if(ncol(target_scale)==1){
      #   target_scale = cbind(target_scale,train_scale[,1:2])
      #   train_scale = train_scale[,-c(1:2)]
      # }

      # all_scale = cbind(train_scale, target_scale)

      # e_current_batch = all_scale
      corrs_train[[current_batch]] = cor(t(train_scale), method = "spearman")
      corrs_target[[current_batch]] = cor(t(target_scale), method = "spearman")
      # corrs[[current_batch]][is.na(corrs[[current_batch]])] = 0
    }

    pred = parallel::parSapply(cl, X = 1:nrow(all), function(j,all,batch.,ranger, sampleType., time., num,corrs_train,corrs_target,minus){
      # cs = c()
      # for(j in 1:nrow(all)){
      # print(j)
      # j = j+1
      # print(j)
      normalized  = rep(0, ncol(all))
      qc_train_value = list()
      qc_predict_value = list()
      sample_value = list()
      sample_predict_value = list()

      for(b in 1:length(levels(batch.))){
        current_batch = levels(batch.)[b]
        current_time = time.[batch. %in% current_batch]
        e_current_batch = all[,batch.%in%current_batch]
        corr_train = corrs_train[[current_batch]]
        corr_target = corrs_target[[current_batch]]


        corr_train_order = order(abs(corr_train[,j]),decreasing = TRUE) #!!! try using corr_target_order[1:num]
        corr_target_order = order(abs(corr_target[,j]),decreasing = TRUE)

        sel_var = c()
        l = num
        # print(intersect(corr_train_order[1:l], corr_target_order[1:l]))
        while(length(sel_var)<(num)){
          sel_var = intersect(corr_train_order[1:l], corr_target_order[1:l])
          sel_var = sel_var[!sel_var == j]
          l = l+1
        }

        train.index_current_batch = sampleType.[batch.%in%current_batch]
        # train_data_y = scale(e_current_batch[j, train.index_current_batch=='qc'],scale=F) #!!! trying to use different scale. Scale to the test_y

        # remove_outlier(e_current_batch[j, train.index_current_batch=='qc'])

        factor = sd(e_current_batch[j, train.index_current_batch=='qc'])/sd(e_current_batch[j, !train.index_current_batch=='qc'])
        if(factor==0 | is.nan(factor) | factor<1 | is.na(factor)){#!!!
          train_data_y = scale(e_current_batch[j, train.index_current_batch=='qc'],scale=F)
        }else{
          # print(j)
          # print("!!")
          if(sum(train.index_current_batch=='qc')*2>=sum(!train.index_current_batch=='qc')){
            train_data_y = (e_current_batch[j, train.index_current_batch=='qc'] - mean(e_current_batch[j, train.index_current_batch=='qc']))/factor ### need to be careful with outlier!
          }else{
            train_data_y = scale(e_current_batch[j, train.index_current_batch=='qc'],scale=F)
          }
        }

        train_data_x = apply(e_current_batch[sel_var, train.index_current_batch=='qc'],1,scale)

        if(is.null(dim(e_current_batch[sel_var, !train.index_current_batch=='qc']))){
          test_data_x = t(scale(e_current_batch[sel_var, !train.index_current_batch=='qc']))
        }else{
          test_data_x = apply(e_current_batch[sel_var, !train.index_current_batch=='qc'],1,scale)
        }

        train_NA_index  = apply(train_data_x,2,function(x){
          sum(is.na(x))>0
        })

        train_data_x = train_data_x[,!train_NA_index]
        test_data_x = test_data_x[,!train_NA_index]

        if(!"matrix" %in% class(test_data_x)){ # !!!
          test_data_x = t(test_data_x)
        }

        good_column = apply(train_data_x,2,function(x){sum(is.na(x))==0}) & apply(test_data_x,2,function(x){sum(is.na(x))==0})
        train_data_x = train_data_x[,good_column]
        test_data_x = test_data_x[,good_column]
        if(!"matrix" %in% class(test_data_x)){ # !!!
          test_data_x = t(test_data_x)
        }
        train_data = data.frame(y = train_data_y,train_data_x )

        if(ncol(train_data)==1){# some samples have all QC constent.
          norm = e_current_batch[j,]
          normalized[batch.%in%current_batch] = norm
        }else{
          colnames(train_data) = c("y", paste0("V",1:(ncol(train_data)-1)))

          set.seed(1)
          model = ranger(y~., data = train_data)

          test_data = data.frame(test_data_x)
          colnames(test_data) = colnames(train_data)[-1]

          norm = e_current_batch[j,]

          # plot(e_current_batch[j, train.index_current_batch=='qc'], ylim = c(min(c(e_current_batch[j, train.index_current_batch=='qc'],e_current_batch[j, !train.index_current_batch=='qc'])),max(c(e_current_batch[j, train.index_current_batch=='qc'],e_current_batch[j, !train.index_current_batch=='qc']))))
          # points(predict(model, data = train_data)$prediction+mean(e_current_batch[j,train.index_current_batch=='qc'],na.rm=TRUE), col = 'red')
          # plot(e_current_batch[j,!train.index_current_batch=='qc'], ylim = c(min(c(e_current_batch[j, train.index_current_batch=='qc'],e_current_batch[j, !train.index_current_batch=='qc'])),max(c(e_current_batch[j, train.index_current_batch=='qc'],e_current_batch[j, !train.index_current_batch=='qc']))))
          # points(predict(model,data = test_data)$predictions  + mean(e_current_batch[j, !train.index_current_batch=='qc'],na.rm=TRUE) - mean(predict(model,data = test_data)$predictions), col = 'red')
          #
          # plot(y =  e_current_batch[j,], x = current_time, col = factor(train.index_current_batch), ylim = range(c(e_current_batch[j,], norm)))
          # plot(y = norm, x = current_time, col = factor(train.index_current_batch), ylim = range(c(e_current_batch[j,], norm)))

          if(minus){

            norm[train.index_current_batch=='qc'] = e_current_batch[j, train.index_current_batch=='qc']-((predict(model, data = train_data)$prediction+mean(e_current_batch[j,train.index_current_batch=='qc'],na.rm=TRUE))-mean(all[j,sampleType.=='qc'],na.rm=TRUE))

            norm[!train.index_current_batch=='qc'] = e_current_batch[j,!train.index_current_batch=='qc']-((predict(model,data = test_data)$predictions  + mean(e_current_batch[j, !train.index_current_batch=='qc'],na.rm=TRUE))-(median(all[j,!sampleType.=='qc'],na.rm = TRUE)))

          }else{
            norm[train.index_current_batch=='qc'] = e_current_batch[j, train.index_current_batch=='qc']/((predict(model, data = train_data)$prediction+mean(e_current_batch[j,train.index_current_batch=='qc'],na.rm=TRUE))/mean(all[j,sampleType.=='qc'],na.rm=TRUE))

            # norm[!train.index_current_batch=='qc'] = e_current_batch[j,!train.index_current_batch=='qc']/((predict(model,data = test_data)$predictions  + mean(e_current_batch[j, !train.index_current_batch=='qc'],na.rm=TRUE))/(median(all[j,!sampleType.=='qc'],na.rm = TRUE)))

            norm[!train.index_current_batch=='qc'] = e_current_batch[j,!train.index_current_batch=='qc']/((predict(model,data = test_data)$predictions  + mean(e_current_batch[j, !train.index_current_batch=='qc'],na.rm=TRUE)- mean(predict(model,data = test_data)$predictions))/(median(all[j,!sampleType.=='qc'],na.rm = TRUE)))

            # norm[!train.index_current_batch=='qc'] = e_current_batch[j,!train.index_current_batch=='qc']/((predict(model,data = test_data)$predictions  + mean(e_current_batch[j, !train.index_current_batch=='qc'],na.rm=TRUE))/(median(e_current_batch[j,!train.index_current_batch=='qc'],na.rm = TRUE)))

          }

          # norm[!train.index_current_batch=='qc'] =(e_current_batch[j,!train.index_current_batch=='qc'])/((predict(model, data = test_data)$prediction + mean(e_current_batch[j,!train.index_current_batch=='qc'],na.rm=TRUE))/mean(e_current_batch[j,!train.index_current_batch=='qc'],na.rm=TRUE))

          norm[!train.index_current_batch=='qc'][norm[!train.index_current_batch=='qc']<0]=e_current_batch[j,!train.index_current_batch=='qc'][norm[!train.index_current_batch=='qc']<0] # fix negative value

          # plot(p$time[batch.%in%b][!train.index_current_batch=='qc'], (e_current_batch[j,!train.index_current_batch=='qc'])/((predict(model,data = test_data)$predictions  + mean(e_current_batch[j, train.index_current_batch=='qc'],na.rm=TRUE))/(median(e_current_batch[j,!train.index_current_batch=='qc'],na.rm = TRUE))))

          norm[train.index_current_batch=='qc'] = norm[train.index_current_batch=='qc']/(median(norm[train.index_current_batch=='qc'],na.rm=TRUE)/median(all[j,sampleType.=='qc'],na.rm=TRUE)) #!!! putting all to the same batch level.

          norm[!train.index_current_batch=='qc'] = norm[!train.index_current_batch=='qc']/(median(norm[!train.index_current_batch=='qc'],na.rm=TRUE)/median(all[j,!sampleType.=='qc'],na.rm=TRUE))

          norm[!is.finite(norm)] = rnorm(length(norm[!is.finite(norm)]),sd = sd(norm[is.finite(norm)],na.rm=TRUE)*0.01)# fix infinite

          out = boxplot.stats(norm, coef = 3)$out
          attempt = ((e_current_batch[j,!train.index_current_batch=='qc'])-((predict(model,data = test_data)$predictions  + mean(e_current_batch[j, !train.index_current_batch=='qc'],na.rm=TRUE))-(median(all[j,!sampleType.=='qc'],na.rm = TRUE))))[norm[!train.index_current_batch=='qc']%in%out];
          if(length(out)>0 & length(attempt)>0){
            if(mean(out)>mean(norm)){
              if(mean(attempt)<mean(out)){
                norm[!train.index_current_batch=='qc'][norm[!train.index_current_batch=='qc']%in%out] =  attempt# !!! this may not help deal with outlier effect..
              }
            }else{
              if(mean(attempt)>mean(out)){
                norm[!train.index_current_batch=='qc'][norm[!train.index_current_batch=='qc']%in%out] =  attempt# !!! this may not help deal with outlier effect..
              }
            }
          }

          norm[!train.index_current_batch=='qc'][norm[!train.index_current_batch=='qc']<0]=e_current_batch[j,!train.index_current_batch=='qc'][norm[!train.index_current_batch=='qc']<0]

          normalized[batch.%in%current_batch] = norm

          # points(current_time, norm, pch = (as.numeric(factor(train.index_current_batch))-1)*19, col = "blue", cex = 0.7)

          # qc_train_value[[b]] = train_data_y + mean(e_current_batch[j, train.index_current_batch=='qc'])
          # qc_predict_value[[b]] = predict(model,data = train_data)$predictions + mean(e_current_batch[j, train.index_current_batch=='qc'])
          # sample_value[[b]] = e_current_batch[j,!train.index_current_batch=='qc']
          # sample_predict_value[[b]] = predict(model,data = test_data)$predictions  + mean(e_current_batch[j, !train.index_current_batch=='qc'])
        }

      }

      # par(mfrow=c(1,2))
      # ylim = c(min(e[j,],norm), max(e[j,],norm))
      # plot(time.[sampleType.=='qc'], unlist(qc_train_value),col = "red",ylim = ylim,main=j)
      # points(time.[sampleType.=='qc'],unlist(qc_predict_value),col = "yellow")
      #
      # points(time.[!sampleType.=='qc'],unlist(sample_value),col = "blue")
      # points(time.[!sampleType.=='qc'],unlist(sample_predict_value),col = "green")
      #
      # plot(time.,normalized, col = factor(sampleType.), ylim = ylim,main=f$label[j])
      #
      # j = j + 1
      #
      # plot(x = time., y = normalized, col = factor(sampleType.), ylim = range(c(normalized, all[j,])))
      # plot(x = time., y = all[j,], col = factor(sampleType.), ylim = range(c(normalized, all[j,])))
      # o = normalized
      c = (median(normalized[sampleType.=="sample"])+(median(all[j,sampleType.=="qc"])-median(all[j,!sampleType.=="qc"]))/sd(all[j,!sampleType.=="qc"]) * sd(normalized[sampleType.=="sample"]))/median(normalized[!sampleType.=="sample"])
      normalized[sampleType.=="qc"] = normalized[sampleType.=="qc"] * ifelse(c>0,c,1)
      # cs[j] = c

      # print(j)
      # }#!!!!

      return(normalized)
    },all,batch.,ranger, sampleType., time., num,corrs_train,corrs_target,minus)


    normed = t(pred)

    normed_target = normed[,!sampleType.=='qc']

    for(i in 1:nrow(normed_target)){ # fix NA
      normed_target[i,is.na(normed_target[i,])] = rnorm(sum(is.na(normed_target[i,])), mean = min(normed_target[i,!is.na(normed_target[i,])], na.rm = TRUE), sd = sd(normed_target[i,!is.na(normed_target[i,])])*0.1)
    }
    for(i in 1:nrow(normed_target)){ # fix negative value
      normed_target[i,normed_target[i,]<0] = runif(1) * min(normed_target[i,normed_target[i,]>0], na.rm = TRUE)
    }

    normed_train = normed[,sampleType.=='qc']

    for(i in 1:nrow(normed_train)){ # fix NA
      normed_train[i,is.na(normed_train[i,])] = rnorm(sum(is.na(normed_train[i,])), mean = min(normed_train[i,!is.na(normed_train[i,])], na.rm = TRUE), sd = sd(normed_train[i,!is.na(normed_train[i,])])*0.1)
    }
    for(i in 1:nrow(normed_train)){ # fix negative value
      normed_train[i,normed_train[i,]<0] = runif(1) * min(normed_train[i,normed_train[i,]>0], na.rm = TRUE)
    }
    return(list(normed_train=normed_train,normed_target=normed_target))
  }

  cat("\nRunning 'SERRF' ...\n")
  serrf_normalized = e
  serrf_normalized_modeled = serrfR(train = e_qc, target = e_sample, num = num,batch. = factor(c(p_qc$batch, p_sample$batch)),time. = c(p_qc$time, p_sample$time),sampleType. = c(p_qc$sampleType, p_sample$sampleType),minus,cl)
  serrf_qc = serrf_normalized_modeled$normed_train
  colnames(serrf_qc) = colnames(e_qc)
  serrf_sample = serrf_normalized_modeled$normed_target
  colnames(serrf_sample) = colnames(e_sample)
  serrf_cross_validated_qc = e_qc
  serrf_validates = list()
  if(with_validate){
    for(validate_type in validate_types){

      serrf_validates[[validate_type]] = serrfR(train = e_qc, target = e_validates[[validate_type]], num = num,batch. = factor(c(p_qc$batch, p_validates[[validate_type]]$batch)),time. = c(p_qc$time, p_validates[[validate_type]]$time),sampleType. = rep(c("qc","sample"),c(nrow(p_qc),nrow(p_validates[[validate_type]]))),minus,cl)$normed_target
      colnames(serrf_validates[[validate_type]]) = colnames(e_validates[[validate_type]])
    }
    normalized_dataset = aggregate_e(serrf_qc,serrf_sample,serrf_validates)
  }else{
    normalized_dataset = aggregate_e(serrf_qc,serrf_sample,NULL)
  }
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
