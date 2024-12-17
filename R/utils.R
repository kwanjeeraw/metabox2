## Utility functions ##

#'Check serrf input
#'@description The function check input data for serrf normalization.
#'@usage check_serrf(p)
#'@param p data.table of sample information
#'@details The function is used inside normalize_input_data_byqc.
#'@return text message
#'@examples
#'#check_serrf(samp_info)
#'@export
check_serrf <- function(p){
  outtxt = ""
  # check if 'qc' or 'sample' is not in sampleType
  if(any(!c('qc','sample') %in% p$sampleType)){
    outtxt = paste0(outtxt,"\nThe 'sampleType' must contain at least 'qc' and 'sample'. See example data for more information.")
  }
  # check if time is in the datasheet
  if(!"time" %in% colnames(p)){
    outtxt = paste0(outtxt,"\nYour data must have 'injectionOrder'. See example data for more information.")
  }
  # check if time has duplicated value.
  if(any(duplicated(p$time))){
    outtxt = paste0(outtxt,"\nData has duplicated 'injectionOrder' values. 'injectionOrder' of each sample should be unique.")
  }
  # check if batch is in the datasheet
  if(!"batch" %in% colnames(p)){
    outtxt = paste0(outtxt,"\nYour data must have 'batch'. See example data for more information.")
  }
  # check if any batch has too little qc
  if(any(table(p$batch, p$sampleType)[,'qc']<2)){
    outtxt = paste0(outtxt,"\nSome batches has a small number of QC that is not enough for training model. Each batch should have at least 2 QCs.")
  }
  return(outtxt)
}

#'Combine p-values
#'@description The function can combine p-values using Fisher's method.
#'@usage combine_pvals(x)
#'@param x numeric vectors of p-values
#'@return numeric number of combined p-values
#'@author Kwanjeera W \email{kwanich@@ucdavis.edu}
#'@references Fisher R. (1932) Statistical methods for research workers. Oliver and Boyd, Edinburgh.
#'@seealso \code{\link{pchisq}}
#'@examples
#'#result <- combine_pvals(c(0.01,0.005,0.1))
#'@export
combine_pvals <- function(x){
  x = x[!x<=0]
  lt <- log(x) #log transform
  chisqu <- (-2) * sum(lt)
  degreeof <- 2 * length(lt)
  pchisq(chisqu, degreeof, lower.tail = FALSE)
}

#'Install package
#'@description Check and install required packages.
#'@usage install_pkgs(pkg_ls)
#'@param pkg_ls text or list indicating the name of a package.
#'@return list of unable to install packages, if any.
#'@examples
#'#install_pkgs("MUVR")
#'@export
install_pkgs <- function(pkg_ls){
  need_pkg = pkg_ls[!(pkg_ls %in% installed.packages()[,"Package"])]
  if(length(need_pkg)){
    temp_ls = sapply(need_pkg, install.packages,dependencies = TRUE) #install from CRAN
    temp_pkg = pkg_ls[!(pkg_ls %in% installed.packages()[,"Package"])]
    if(length(temp_pkg)){
      cat("\nMissing the required package, trying Bioconductor ...\n")
      temp_ls = sapply(temp_pkg, function(x){
        tryCatch(BiocManager::install(x,ask=TRUE,update=FALSE), error=function(e){x}) #install from Bioconductor
      })
      unin_pkg = pkg_ls[!(pkg_ls %in% installed.packages()[,"Package"])]
    }else{
      unin_pkg = NULL
    }
  }else{
    unin_pkg = NULL
  }
  if(length(unin_pkg)) cat("\nCould not find package: ",paste(as.character(unlist(unin_pkg)),collapse=", "),"\n")
  return(unin_pkg)
}

#'Generate report
#'@description Generate report after data processing or data analysis.
#'@usage generate_report(METBObj=NULL, datsummary1=NULL, datsummary2=NULL,
#'note=NULL, reportfile, savedir=getwd())
#'@param METBObj METBObj object returned from data processing functions, including \code{impute_missing_data},
#'\code{normalize_input_data_bydata} and \code{normalize_input_data_byqc}, \code{scale_input_data}, \code{set_input_obj}, \code{transform_input_data}.
#'@param datsummary1 list of data returned from data analysis functions, including \code{comb_overrep_analyze}, \code{correlation_analyze},
#'\code{enrichment_analyze}, \code{lme_analyze}, \code{multiv_analyze}, \code{overrep_analyze}, \code{run_muvr}, \code{test_multinormality},
#'\code{univ_analyze}, \code{mbplsda_analyze}.
#'@param datsummary2 list of data returned from data analysis function, including \code{test_uninormality}.
#'@param note text of additional information.
#'@param reportfile name of a report file specifically for each data processing or analysis function,
#'including correlation_report, impute_data_report, lme_report, mbplsda_report, multivariate_report, biomarker_report,
#'normalize_bydata_report, normalize_byqc_report, scale_data_report, transform_data_report and univariate_report.
#'@param savedir path to a directory for saving the report. Default is the current working directory.
#'@return pdf file.
#'@examples
#'#sugar_dt = set_input_obj(sugar, 1,2,5)
#'#out=multiv_analyze(sugar_dt)
#'#generate_report(sugar_dt,datsummary1=out,reportfile="multivariate_report")
#'@export
generate_report = function(METBObj=NULL, datsummary1=NULL, datsummary2=NULL, note=NULL, reportfile, savedir=getwd()) {
  rmarkdown::render(
    paste0(system.file("shiny/reports", package = "metabox2"), "/", reportfile, ".Rmd"), params = list(
      METBObj = METBObj,datsummary1 = datsummary1, datsummary2=datsummary2,note=note
    ),
    output_format = "pdf_document",output_file = paste0(reportfile,".pdf"), output_dir = savedir)
}

#'Get posthoc summary
#'@description Get posthoc summary from ANOVA.
#'@usage get_posthoc_summary(aov_res, cutoff=0.05)
#'@param aov_res list of ANOVA results containing posthoc tests from \code{univ_analyze()}.
#'@param cutoff a number indicating cutoff for statistical significance (Default = 0.05).
#'@return data frame of test results.
#'@examples
#'#sugar_dt = set_input_obj(sugar, 1,2,5)
#'#out=univ_analyze(sugar_dt, factor2Col=3, doposthoc=TRUE)
#'#test_out=get_posthoc_summary(out)
#'@export
get_posthoc_summary <- function(aov_res, cutoff=0.05){
  #Initialize parameters
  dat = aov_res$posthoc_table; #working data
  if(nrow(dat)==0){
    cat('\nPosthoc results are missing. Returning no data.\n')
    return(data.frame())
  }else{
    cat('\nFormating posthoc output ...\n')
    significant = apply(dat, 1, function(x){
      inx = x < cutoff
      paste(colnames(dat)[inx], collapse = " | ")
    })
    out_dat = cbind(as.data.frame(significant),dat)
    cat('\n',sum(out_dat$significant!=""),'variables are significant at least in one condition.\n')
    return(out_dat)
  }
}

#'Get percent RSD of QC samples
#'@description Calculate the relative standard deviation (RSD) of QC samples.
#'@usage get_rsd(METBObj,sampleType=FALSE)
#'@param METBObj METBObj object contains list of data.
#'@param sampleType a column number/index indicating sample types in the \code{METBObj$inputdata} data frame.
#'This parameter is required for serrf and loess methods.
#'
#'Require "qc" to indicate pooled/quality control (QC) samples.
#'
#'@details The function calculates median RSD, mean RSD and
#'RSD of each variable in the QC samples using RSD = sd/mean.
#'@return data frame.
#'@examples
#'#load('GitHub/metabox-ml/data/qc_batch_example.RData')
#'#out=get_rsd(qc_batch_example,sampleType=3)
#'@export
get_rsd <- function(METBObj, sampleType=FALSE){
  #Initialize parameters
  dat = METBObj$X; rsd_table = data.frame(); #working data
  inputdat = METBObj$inputdata #working data
  inputdat[,sampleType]=tolower(inputdat[,sampleType]) #use lower case
  if(any(!c('qc') %in% inputdat[,sampleType])){
    cat("\nERROR! The 'sampleType' must contain at 'qc'. Please see example data for more information. Data was not computed.\n")
    return(rsd_table)
  }
  #Check sampleType column
  if (is.factor(inputdat[,sampleType])) {
    F1 = inputdat[,sampleType] #working data
  }else{
    F1 = as.factor(inputdat[,sampleType]) #working data
  }
  calc_grmean = aggregate(dat, list(F1), FUN = mean, na.rm=TRUE) #calculate mean for each level of F1
  calc_grsd = aggregate(dat, list(F1), FUN = sd, na.rm=TRUE) #calculate sd for each level of F1
  qc_ind = which(calc_grmean[,1] == "qc") #qc index
  calc_rsd = (calc_grsd[qc_ind,-1]/abs(calc_grmean[qc_ind,-1]))*100 #calculate percent rsd for qc
  mean_rsd = mean(unlist(calc_rsd[1,]),na.rm=TRUE) #calculate mean rsd
  median_rsd = median(unlist(calc_rsd[1,]),na.rm=TRUE) #calculate median rsd
  rsd_table = cbind(mean_rsd,median_rsd,calc_rsd)
  return(rsd_table)
}

#'Get statistical summary of METBObj
#'@description Get statistical summary of METBObj.
#'@usage get_stat_summary(METBObj)
#'@param METBObj METBObj object contains list of data.
#'@details The function calculates median, mean, SD, coefficient of variation (cv) and normality.
#'@return data frame of statistical summary.
#'@examples
#'#sugar_dt = set_input_obj(sugar, 1,2,5)
#'#out=get_stat_summary(sugar_dt)
#'@export
get_stat_summary <- function(METBObj){
  #Initialize parameters
  dat = METBObj$X; stat_summ = data.frame(); #working data
  normshaptest = function(x){ #shapiro.test by group
    tryCatch(shapiro.test(x)$p.value, error=function(e){0})
  }
  if(is.numeric(METBObj$Y)){#numeric response variables
    cat('\nFind a numeric vector. Assuming response variables are continuous.')
    calc_median = apply(dat, 2, function(x) median(x, na.rm=TRUE)) #calculate median across levels
    calc_mean = apply(dat, 2, function(x) mean(x, na.rm=TRUE)) #calculate mean across levels
    calc_sd = apply(dat, 2, function(x) sd(x, na.rm=TRUE)) #calculate sd across levels
    stat_summ = data.frame(all_median=calc_median, all_mean=calc_mean, all_sd=calc_sd, check.names = FALSE)
  }else{#categorical response variables
    #Check category/factor column
    if (is.factor(METBObj$Y)) {
      F1 = METBObj$Y #working data
      cat('\nFind a factor F1 of ',nlevels(F1),' levels.',sep='')
    }else{
      F1 = as.factor(METBObj$Y) #working data
      cat('\nThe category/factor column is', class(METBObj$Y), '.\nConverting to a factor F1 of ',nlevels(F1),' levels.',sep='')
    }
    calc_median = apply(dat, 2, function(x) median(x, na.rm=TRUE)) #calculate median across levels
    calc_mean = apply(dat, 2, function(x) mean(x, na.rm=TRUE)) #calculate mean across levels
    calc_sd = apply(dat, 2, function(x) sd(x, na.rm=TRUE)) #calculate sd across levels
    calc_grmedian = aggregate(dat, list(F1), FUN = median, na.rm=TRUE) #calculate median for each level of F1
    calc_grmedian=t(calc_grmedian[-1])
    colnames(calc_grmedian) = paste0(levels(F1),"_median")
    calc_grmean = aggregate(dat, list(F1), FUN = mean, na.rm=TRUE) #calculate mean for each level of F1
    calc_grmean=t(calc_grmean[-1])
    colnames(calc_grmean) = paste0(levels(F1),"_mean")
    calc_grsd = aggregate(dat, list(F1), FUN = sd, na.rm=TRUE) #calculate sd for each level of F1
    calc_grsd=t(calc_grsd[-1])
    colnames(calc_grsd) = paste0(levels(F1),"_sd")
    calc_grcv = abs(data.frame((matrix(mapply('/',calc_grsd,calc_grmean),ncol = nlevels(F1),byrow = F)),row.names = row.names(calc_grmean))) #calculate abs coefficient of variation for each level of F1
    colnames(calc_grcv) = paste0(levels(F1),"_cv")
    calc_grnorm = aggregate(dat, list(F1), FUN = normshaptest) #calculate normality for each level of F1
    calc_grnorm=t(calc_grnorm[-1])
    colnames(calc_grnorm) = paste0(levels(F1),"_pnormal")
    # calc_fc = apply(calc_grmean,2, function(x){x/calc_mean}) #foldChange = each_level/all_mean
    # rm <- colMeans(dat);rm_mean <- sweep(dat, 2, rm);sx <- apply(rm_mean, 2, sd);
    # calc_fc <- sweep(rm_mean, 2, sx, "/"); #scale z-score
    # calc_fcgrmean = aggregate(calc_fc, list(F1), FUN = mean, na.rm=TRUE) #calculate z-score mean for each level of F1
    # calc_fcgrmean=t(calc_fcgrmean[-1])
    # colnames(calc_fcgrmean) = paste0(levels(F1),"_mean_zscore")
    # calc_fcgrsd = aggregate(calc_fc, list(F1), FUN = sd, na.rm=TRUE) #calculate z-score sd for each level of F1
    # calc_fcgrsd=t(calc_fcgrsd[-1])
    # colnames(calc_fcgrsd) = paste0(levels(F1),"_sd_zscore")
    stat_summ = data.frame(all_median=calc_median, all_mean=calc_mean, all_sd=calc_sd,
                           calc_grmedian, calc_grmean, calc_grsd, calc_grcv, calc_grnorm, check.names = FALSE)
  }
  return(stat_summ)
}

#'Shapiro-Wilk test of univariate normality
#'@description Perform Shapiro-Wilk test on each variable.
#'@usage test_uninormality(METBObj)
#'@param METBObj METBObj object contains list of data.
#'@return list of test results.
#'@examples
#'#test_uninormality(METBObj)
#'@export
test_uninormality <- function(METBObj){
  dat = METBObj$X #working data
  test_res = apply(dat, 2, function(x){
    tryCatch(shapiro.test(x)$p.value, error=function(e){0})
    })
  norm_count = sum(test_res > 0.05)
  norm_mean = mean(test_res)
  cat("\nTotal normally distributed variables:",norm_count,"from",ncol(dat),"\n")
  return(list(p_value=test_res,norm_count=norm_count, norm_mean=norm_mean))
}

#'Shapiro-Wilk test of multivariate normality
#'@description Perform Shapiro-Wilk test on all variables.
#'@usage test_multinormality(METBObj)
#'@param METBObj METBObj object contains list of data.
#'@return list of test results.
#'@examples
#'#test_multinormality(METBObj) #Not exported function
test_multinormality <- function(METBObj){
  dat = cbind(METBObj$inputdata[,METBObj$classCol], METBObj$X) #working data
  mvn_all = tryCatch({
    mvn_test = MVN::mvn(dat[,-1], mvnTest = 'royston', univariateTest = 'SW') #mvn test all
    norm_count = dplyr::count(mvn_test$univariateNormality, Normality)
    cat("\nMultivariate normality:",mvn_test$multivariateNormality$MVN,"; p-value =",mvn_test$multivariateNormality$`p value`,
        "\nTotal normally distributed variables:",norm_count[order(norm_count$Normality),"n"][2],"from",ncol(dat[,-1]),"\n")
    mvn_test
  },
  error=function(e){
    message(e)
    cat("\nERROR! Cannot perform multivariate normality test.\n")
    list()
  })
  mvn_class = tryCatch({
    MVN::mvn(dat, mvnTest = 'royston', subset=colnames(dat)[1], univariateTest = 'SW') #mvn test by class
  },
  error=function(e){
    message(e)
    cat("\nERROR! Cannot perform normality test on subgroups.\n")
    list()
  })
  return(list(mvn_all=mvn_all, mvn_class=mvn_class))
}

#'Levene's test for homogeneity of variance across groups
#'@description Perform Levene's test on each variable.
#'@usage test_uniequalvar(METBObj)
#'@param METBObj METBObj object contains list of data.
#'@return list of test results.
#'@examples
#'#test_uniequalvar(METBObj)
#'@export
test_uniequalvar <- function(METBObj){
  dat = METBObj$X #working data
  test_res = apply(dat, 2, function(x){car::leveneTest(x ~ METBObj$Y)[1,3]})
  norm_count = sum(test_res > 0.05)
  norm_mean = mean(test_res)
  cat("\nTotal equal variances:",norm_count,"from",ncol(dat),"variables\n")
  return(list(p_value=test_res,norm_count=norm_count, norm_mean=norm_mean))
}

#'Box's M-test for homogeneity of multivariate data
#'@description Perform Box's M-test on all variables.
#'@usage test_multiequalvar(METBObj)
#'@param METBObj METBObj object contains list of data.
#'@return list of test results.
#'@examples
#'#test_multiequalvar(METBObj) #Not exported function
test_multiequalvar <- function(METBObj){
  dat = METBObj$X #Working data
  F1 = factor(METBObj$Y) #Working data
  test_res = biotools::boxM(dat,F1)
  return(boxm_all=test_res)
}

## Plot functions ##

#'Box plot overview
#'@description Provide box plot.
#'@usage boxplot_overview(METBObj, plotvar=TRUE, plot_title="")
#'@param METBObj METBObj object contains list of data.
#'@param plotvar a logical variable indicating box plot of variables/features (TRUE) or samples (FALSE).
#'@param plot_title text indicating plot title.
#'@return ggplot object.
#'@examples
#'#boxplot_overview(METBObj)
#'@import ggplot2
#'@import RColorBrewer
#'@export
boxplot_overview <- function(METBObj, plotvar=TRUE, plot_title=""){
  if(plotvar){
    if (ncol(METBObj$X)<2){#need at least 2 points
      cat("\nError! in density.default, need at least 2 data points.\nReturn no plot.\n")
      return(FALSE)
    }
    dat = cbind(Sample_ID=METBObj$ID, METBObj$X) #working data
    plotdata = reshape2::melt(dat, id = c("Sample_ID"))#for ggplot
    ggplot(plotdata, aes(x = variable, y = value, label=Sample_ID)) + geom_boxplot(outlier.shape = NA)+
      geom_jitter(width=0.1, alpha=0.2) + theme_light() + ggtitle(plot_title) +
      theme(axis.title=element_blank(), axis.text=element_text(size=10)) + coord_flip() #by variables
  }else{
    if (nrow(METBObj$X)<2){#need at least 2 points
      cat("\nError! in density.default, need at least 2 data points.\nReturn no plot.\n")
      return(FALSE)
    }
    dat = cbind(Variable_ID=colnames(METBObj$X), data.frame(t(METBObj$X))) #working data
    plotdata = reshape2::melt(dat, id = c("Variable_ID"))#for ggplot
    ggplot(plotdata, aes(x = variable, y = value, label=Variable_ID)) + geom_boxplot(outlier.shape = NA)+
      geom_jitter(width=0.1, alpha=0.2) + theme_light() + ggtitle(plot_title) +
      theme(axis.title=element_blank(), axis.text=element_text(size=10)) + coord_flip() #by samples
  }
}

#'Box plot by category/factor F1
#'@description Provide box plot of a variable/feature by category/factor F1.
#'@usage boxplot_byF1(METBObj, xCol=1, factorLv1)
#'@param METBObj METBObj object contains list of data.
#'@param xCol a column number or column name indicating which variable/feature in the \code{x_input} data frame to be plotted.
#'@param factorLv1 a vector of category/factor levels.
#'@return ggplot object.
#'@examples
#'#boxplot_byF1(METBObj, 1, METBObj$Y)
#'@import ggplot2
#'@import RColorBrewer
#'@export
boxplot_byF1 <- function(METBObj, xCol=1, factorLv1){
  #Check category/factor
  if (is.factor(factorLv1)) {
    F1 = factorLv1 #working data
    cat('\nFind a factor of ',nlevels(factorLv1),' levels.',sep='')
  }else{
    F1 = as.factor(factorLv1) #working data
    cat('\nConverting to a factor of ',nlevels(factorLv1),' levels.',sep='')
  }
  dat = METBObj$X
  #Set plot title
  if(is.numeric(xCol)){#xCol is index
    xnames = colnames(dat)[xCol]
  }else{#xCol is name
    xnames = as.character(xCol)
  }
  colnames(dat)[xCol] = 'value' #change column name
  plotdata = cbind(Sample_ID=METBObj$ID, dat, Groups=F1) #working data; add the factor column at the end of the dataframe
  numlevels = nlevels(F1) #no. of levels
  if(numlevels <= 8){
    grcolors = mbcolors
  }else{
    grcolors = colorRampPalette(mbcolors)(numlevels)
  }
  ggplot(plotdata, aes(x=Groups, y = value, color=Groups, label=Sample_ID)) + geom_boxplot(outlier.shape = NA) +
    scale_color_manual(values = grcolors) + geom_jitter(width=0.1, alpha=0.2) + theme_light() + ggtitle(xnames) + labs(color="") +
    theme(plot.title = element_text(hjust = 0.5, size = 12), axis.title=element_blank(),
          axis.text=element_text(size=10), axis.text.x = element_text(angle=90,hjust=1,vjust=0.2))
}

#'Box plot by category/factor F1 and F2
#'@description Provide box plot of a variable/feature by category/factor F1 and F2.
#'@usage boxplot_byF1F2(METBObj, xCol=1, factorLv1, factorLv2)
#'@param METBObj METBObj object contains list of data.
#'@param xCol a column number or column name indicating which variable/feature in the \code{x_input} data frame to be plotted.
#'@param factorLv1 a vector of F1 category/factor levels.
#'@param factorLv2 a vector of F2 category/factor levels.
#'@return ggplot object.
#'@examples
#'#sugar_dt = set_input_obj(sugar, 1,2,5)
#'#boxplot_byF1F2(sugar_dt, 1, sugar_dt$Y, sugar_dt$inputdata$time)
#'@import ggplot2
#'@import RColorBrewer
#'@export
boxplot_byF1F2 <- function(METBObj, xCol=1, factorLv1, factorLv2){
  #Check category/factor
  if (is.factor(factorLv1)) {
    F1 = factorLv1 #working data
    cat('\nFind a factor F1 of ',nlevels(factorLv1),' levels.',sep='')
  }else{
    F1 = as.factor(factorLv1) #working data
    cat('\nConverting to a factor F1 of ',nlevels(factorLv1),' levels.',sep='')
  }
  if (is.factor(factorLv2)) {
    F2 = factorLv2 #working data
    cat('\nFind a factor F2 of ',nlevels(factorLv2),' levels.',sep='')
  }else{
    F2 = as.factor(factorLv2) #working data
    cat('\nConverting to a factor F2 of ',nlevels(factorLv2),' levels.',sep='')
  }
  dat = METBObj$X
  #Set plot title
  if(is.numeric(xCol)){#xCol is index
    xnames = colnames(dat)[xCol]
  }else{#xCol is name
    xnames = as.character(xCol)
  }
  colnames(dat)[xCol] = 'value' #change column name
  plotdata = cbind(Sample_ID=METBObj$ID, dat, F1=F1, Groups=F2) #working data; add the F1 and F2 columns at the end of the dataframe
  numlevels = max(nlevels(F1), nlevels(F2)) #no. of levels for color
  if(numlevels <= 8){
    grcolors = mbcolors
  }else{
    grcolors = colorRampPalette(mbcolors)(numlevels)
  }
  ggplot(plotdata, aes(x=F1, y = value, color=Groups, label=Sample_ID)) + geom_boxplot(outlier.shape = NA) + geom_point(position=position_jitterdodge(), alpha=0.2) +
    scale_color_manual(values = grcolors) + theme_light() + ggtitle(xnames) + labs(color="") +
    theme(plot.title = element_text(hjust = 0.5, size = 12), axis.title=element_blank(),
          axis.text=element_text(size=10), axis.text.x = element_text(angle=90,hjust=1,vjust=0.2))
}

#'Combined statistics plot
#'@description Provide plot of outputs from 2 statistical analyses.
#'@usage combine_statplot(x_data, y_data, x_cutoff=0, y_cutoff=0, ptsize=3, plot_title="")
#'@param x_data a numerical vector or a data frame of statistical values from x analysis, see details.
#'@param y_data a numerical vector or a data frame of statistical values from y analysis, see details.
#'@param x_cutoff text indicating plot title.
#'@param y_cutoff text indicating plot title.
#'@param ptsize a number of geom_point size.
#'@param plot_title text indicating plot title.
#'@details A plot of outputs from two statistical analyses will generate.
#'A numerical vector or a one-column data frame of statistical values (e.g. p-values or VIP)
#'from each statistical analysis (e.g. t-test or PLS-DA) is needed for x_data and y_data.
#'Significance of a variable is denoted by its statistical value less than x_cutoff or more than y_cutoff.
#'@return ggplot object.
#'@examples
#'#sugar_dt = set_input_obj(sugar, 1,2,5)
#'#uout=univ_analyze(sugar_dt)
#'#mout=multiv_analyze(sugar_dt, method="pls", scale="pareto")
#'#combine_statplot(uout$p_adj,mout$vip_val,x_cutoff = 0.05,y_cutoff = 1)
#'@import ggplot2
#'@export
combine_statplot <- function(x_data, y_data, x_cutoff=0, y_cutoff=0, ptsize=3, plot_title=""){
  if(!is.data.frame(x_data)){
    x_data = data.frame(x_data)
  }
  if(!is.data.frame(y_data)){
    y_data = data.frame(y_data)
  }
  plotdata = merge(x_data, y_data, by = "row.names")
  colnames(plotdata) = c("Variable","x_stat","y_stat")
  plotdata$x_sum = unlist(lapply(plotdata$x_stat, FUN = function(x){x < x_cutoff}))
  plotdata$y_sum = unlist(lapply(plotdata$y_stat, FUN = function(x){x > y_cutoff}))
  plotdata$stat_sum = factor(paste0(plotdata$x_sum,"_",plotdata$y_sum),
                             levels = c("TRUE_TRUE", "TRUE_FALSE", "FALSE_TRUE","FALSE_FALSE"))
  ggplot(plotdata, aes(x_stat, y_stat, label=Variable, color=stat_sum, group = stat_sum)) + geom_point(alpha=0.5, size=ptsize) +
    geom_hline(yintercept=y_cutoff, linetype=5, colour="#3F3A65", alpha=0.7) +
    geom_vline(xintercept = x_cutoff, linetype=5, colour="#3F3A65", alpha=0.7) +
    scale_color_manual(name = "",
                       values = c("TRUE_TRUE"="#FF003A", "TRUE_FALSE"="#004DDF", "FALSE_TRUE"="#00BF15","FALSE_FALSE"="#464A44"),
                       label = c("Both","Only < x_cutoff","Only > y_cutoff","None")) +
    theme_bw() + ggtitle(plot_title) + labs(x = "x statistics", y = "y statistics") +
    theme(plot.title = element_text(size=12), axis.title=element_text(size=10), axis.text=element_text(size=10))
}

#'Correlation heatmap
#'@description Provide correlation heatmap.
#'@usage corrplot_heatmap(corr_data, plot_title="")
#'@param corr_data data frame result from \code{correlation_analyze()}.
#'The data frame must contain Var1,Var2,coeff,p_value,p_adj columns.
#'@param plot_title text indicating plot title.
#'@return ggplot object.
#'@examples
#'#sugar_dt = set_input_obj(sugar, 1,2,5)
#'#out=correlation_analyze(sugar_dt)
#'#corrplot_heatmap(data.frame(out$corr_data))
#'@import ggplot2
#'@import RColorBrewer
#'@export
corrplot_heatmap <- function(corr_data, plot_title=""){
  #check corr_data colnames
  tmparg = paste0(c("Var1","Var2","coeff","p_value","p_adj"),collapse = "") == paste0(colnames(corr_data),collapse ="")
  if(!tmparg){
    cat("\ncorr_data must contain column names: Var1,Var2,coeff,p_value,p_adj.\nReturn no plot.\n")
    return(FALSE)
  }
  corr_data$coeff=round(corr_data$coeff,2); corr_data$p_value=round(corr_data$p_value,2); corr_data$p_adj=round(corr_data$p_adj,2) #round to 2 points
  ggplot(corr_data, aes(x=Var1, y=Var2, fill=coeff, label=p_adj)) + geom_tile() + theme_light() + ggtitle(plot_title) +
    scale_fill_gradient2(low="blue",mid="white",high="red", limits=c(-1,1)) + labs(x = "", y = "", fill = "Coefficient") +
    theme(plot.title = element_text(size = 12), axis.title=element_text(size=10),
          axis.text=element_text(size=10), axis.text.x = element_text(angle=90,hjust=1,vjust=0.2)) #+ geom_text()
}

#'Density plot
#'@description Provide density overview.
#'@usage densityplot_overview(x_input, plotvar=TRUE, plot_title="")
#'@param x_input data frame of variables/features. Variables are in columns and samples are in rows.
#'@param plotvar a logical variable indicating density plot of variables/features (TRUE) or samples (FALSE).
#'@param plot_title text indicating plot title.
#'@return ggplot object.
#'@examples
#'#densityplot_overview(METBObj$X)
#'@export
densityplot_overview <- function(x_input, plotvar=TRUE, plot_title=""){
  dat = x_input #working data
  if(!is.data.frame(dat)){
    cat("\nError! x_input must be a data frame.\nReturn no plot.\n")
    return(FALSE)
  }
  if(plotvar){
    if (ncol(dat)<2){#need at least 2 points
      cat("\nError! in density.default, need at least 2 data points.\nReturn no plot.\n")
      return(FALSE)
    }
    dens = data.frame(variable=apply(dat, 2, mean, na.rm=TRUE))#variable density
  }else{
    if (nrow(dat)<2){#need at least 2 points
      cat("\nError! in density.default, need at least 2 data points.\nReturn no plot.\n")
      return(FALSE)
    }
    dens = data.frame(variable=apply(dat, 1, mean, na.rm=TRUE))#sample density
  }
  ggplot(dens, aes(x=variable)) + geom_density() + theme_bw() + ggtitle(plot_title) + labs(x = "", y = "Density")
}

#'Intensity plot
#'@description Provide intensity plot of a variable/feature
#'@usage intensityplot_overview(METBObj, xCol=1, classCol=NULL,
#'legend_title="Group", ptsize=3, plot_title="")
#'@param METBObj METBObj object contains list of data.
#'@param xCol a column number or column name indicating which variable/feature in the \code{METBObj$X} data frame to be plotted.
#'@param classCol a column number/index of category/factor column to show on the plot.
#'If not specified, \code{METBObj$classCol} is used by default.
#'@param legend_title text indicating the legend title.
#'@param ptsize a number of geom_point size.
#'@param plot_title text indicating plot title.
#'@return ggplot object.
#'@examples
#'#intensityplot_overview(qc_batch_example, xCol = 200, classCol = 3)
#'@import ggplot2
#'@import RColorBrewer
#'@export
intensityplot_overview <- function(METBObj, xCol=1, classCol=NULL, legend_title="Group", ptsize=3, plot_title=""){
  dat = data.frame(value=METBObj$X[,xCol]) #working data
  if(is.null(classCol)){#METBObj$classCol
    plotdata = cbind(Sample_ID=(METBObj$ID), data.frame(ClassCol=METBObj$inputdata[,METBObj$classCol]), dat) #for ggplot
  }else{#other classCol
    plotdata = cbind(Sample_ID=(METBObj$ID), data.frame(ClassCol=METBObj$inputdata[,classCol]), dat) #for ggplot
  }
  #Check type of category/factor column
  if (is.numeric(plotdata$ClassCol)) {#plot regression
    Y = plotdata$ClassCol #working data
    ggplot(plotdata, aes(Sample_ID, value, color = ClassCol)) + geom_point(size=ptsize, alpha=0.5) +
      scale_colour_gradientn(colours = rainbow(length(Y), start=0.17, end=1)) + theme_minimal() +
      labs(x = "Sample", y = "Intensity", color="") + ggtitle(colnames(METBObj$X)[xCol]) +
      theme(plot.title = element_text(size = 12), axis.title=element_text(size=10),
            axis.text=element_text(size=8), axis.text.x = element_text(angle=90,hjust=1,vjust=0.2))
  }else{#plot category/factor
    Y = as.factor(plotdata$ClassCol) #working data
    numlevels = nlevels(Y) #no. of levels
    if(numlevels <= 8){
      grcolors = mbcolors
    }else{
      grcolors = colorRampPalette(mbcolors)(numlevels)
    }
    ggplot(plotdata, aes(Sample_ID, value, color = ClassCol)) + geom_point(size=ptsize, alpha=0.5) +
      scale_color_manual(values = grcolors) + theme_minimal() +
      labs(x = "Sample", y = "Intensity", color="") + ggtitle(colnames(METBObj$X)[xCol]) +
      theme(plot.title = element_text(size = 12), axis.title=element_text(size=10),
            axis.text=element_text(size=10), axis.text.x = element_text(angle=90,hjust=1,vjust=0.2))
  }
}

#'Interaction plot between category/factor F1 and F2
#'@description Provide interaction plot between category/factor F1 and F2.
#'@usage interactionplot_byF1F2(x_input, xCol=1, factorLv1, factorLv2, ptsize=3)
#'@param x_input data frame of variables/features. Variables are in columns and samples are in rows.
#'@param xCol a column number or column name indicating which variable/feature in the \code{x_input} data frame to be plotted.
#'@param factorLv1 a vector of F1 category/factor levels.
#'@param factorLv2 a vector of F2 category/factor levels.
#'@param ptsize a number of geom_point size.
#'@return ggplot object.
#'@examples
#'#sugar_dt = set_input_obj(sugar, 1,2,5)
#'#interactionplot_byF1F2(sugar_dt$X, 1, sugar_dt$Y, sugar_dt$inputdata$time)
#'@import ggplot2
#'@import RColorBrewer
#'@export
interactionplot_byF1F2 <- function(x_input, xCol=1, factorLv1, factorLv2, ptsize=3){
  if(!is.data.frame(x_input)){
    cat("\nError! x_input must be a data frame.\nReturn no plot.\n")
    return(FALSE)
  }
  #Check category/factor
  if (is.factor(factorLv1)) {
    F1 = factorLv1 #working data
    cat('\nFind a factor F1 of ',nlevels(factorLv1),' levels.',sep='')
  }else{
    F1 = as.factor(factorLv1) #working data
    cat('\nConverting to a factor F1 of ',nlevels(factorLv1),' levels.',sep='')
  }
  if (is.factor(factorLv2)) {
    F2 = factorLv2 #working data
    cat('\nFind a factor F2 of ',nlevels(factorLv2),' levels.',sep='')
  }else{
    F2 = as.factor(factorLv2) #working data
    cat('\nConverting to a factor F2 of ',nlevels(factorLv2),' levels.',sep='')
  }
  dat = cbind(x_input, F1=F1, Groups=F2) #working data; add the F1 and F2 columns at the end of the dataframe
  #Set plot title
  if(is.numeric(xCol)){#xCol is index
    xnames = colnames(dat)[xCol]
  }else{#xCol is name
    xnames = as.character(xCol)
  }
  numlevels = max(nlevels(F1), nlevels(F2)) #no. of levels for color
  if(numlevels <= 8){
    grcolors = mbcolors
  }else{
    grcolors = colorRampPalette(mbcolors)(numlevels)
  }
  statsum = FSA::Summarize(dat[,xCol] ~ F1 + Groups, dat, digits=4)
  ggplot(statsum, aes(x = F1, y = median, colour = Groups, group = Groups)) + geom_point(size=ptsize, alpha=0.5) + geom_line() +
    #geom_errorbar(aes(ymin = mean - sd,ymax = mean + sd),width=0.1) +
    scale_color_manual(values = grcolors) + theme_light() + labs(title = xnames, x = "", y = "Median", color="") +
    theme(plot.title = element_text(hjust = 0.5, size = 12), axis.title=element_text(size=10),
          axis.text=element_text(size=10), axis.text.x = element_text(angle=90,hjust=1,vjust=0.2))
}

#'Loading plot
#'@description Provide loading plot for multivariate analysis.
#'@usage multiv_loadingplot(loading_data, oloading_data=NULL, ptsize=3, plot_title="")
#'@param loading_data a numerical matrix or a data frame of x loadings, see details.
#'@param oloading_data a numerical matrix or a data frame of orthogonal loadings, for OPLS only. Default is \code{NULL}.
#'@param ptsize a number of geom_point size.
#'@param plot_title text indicating plot title.
#'@details A two-column numerical matrix or a two-column data frame of x loadings is required for PCA and PLSDA.
#'@return ggplot object.
#'@examples
#'#adipose_dt = set_input_obj(adipose,2,3,4)
#'#out=multiv_analyze(adipose_dt, method="pls", scale="standard")
#'#multiv_loadingplot(out$loading_val)
#'@import ggplot2
#'@import RColorBrewer
#'@export
multiv_loadingplot <- function(loading_data, oloading_data=NULL, ptsize=3, plot_title=""){
  if(!is.data.frame(loading_data)){
    loading_data = data.frame(loading_data)
  }
  if(is.null(oloading_data)){
    plotdata = loading_data #for ggplot
  }else{
    if(!is.data.frame(oloading_data)){
      oloading_data = data.frame(oloading_data)
    }
    plotdata = cbind(loading_data, oloading_data) #orthogonal loadings for ggplot
  }
  plotdata = cbind(rownames(loading_data),plotdata)
  xlabel=colnames(plotdata)[2]; ylabel=colnames(plotdata)[3] #x,y label
  colnames(plotdata)[1] = 'Sample_ID'; colnames(plotdata)[2] = 'PCX'; colnames(plotdata)[3] = 'PCY' #change column name
  ggplot(plotdata, aes(PCX, PCY, label=Sample_ID)) + geom_point(color="red", alpha=0.5, size=ptsize) +
    geom_hline(yintercept=0, linetype=1) + geom_vline(xintercept = 0, linetype=1) + theme_minimal() + ggtitle(plot_title) +
    theme(plot.title = element_text(size=12), axis.title=element_text(size=10), axis.text=element_text(size=10)) +
    labs(x = xlabel, y = ylabel)
}

#'Loading plot by PC
#'@description Provide loading plot by a principal component.
#'@usage multiv_loadingplot_bypc(loading_data, pc=1, plot_title="")
#'@param loading_data a numerical matrix or a data frame of x loadings or orthogonal loadings.
#'@param pc a number indicating a principal component to plot. 1st PC is plotted by default.
#'@param plot_title text indicating plot title.
#'@return ggplot object.
#'@examples
#'#adipose_dt = set_input_obj(adipose,2,3,4)
#'#out=multiv_analyze(adipose_dt, method="pls", scale="standard")
#'#multiv_loadingplot_bypc(out$loading_val)
#'@import ggplot2
#'@import RColorBrewer
#'@export
multiv_loadingplot_bypc <- function(loading_data, pc=1, plot_title=""){
  plotdata = data.frame(vname=row.names(loading_data),loading=loading_data[,pc]) #for ggplot
  plotdata = plotdata[order(plotdata$loading),] #sort
  plotdata$vname = factor(plotdata$vname, levels = plotdata$vname) #set factor
  ggplot(plotdata, aes(x=vname, y=loading, fill=loading)) + geom_bar(stat = "identity", alpha=0.7) + theme_minimal() +
    scale_fill_gradient(low = "blue", high = "red") + labs(x = "", y = paste("Loading",pc), fill = "Loadings") + ggtitle(plot_title) +
    theme(plot.title = element_text(size=12), axis.title=element_text(size=10),
          axis.text=element_text(size=10), axis.text.x = element_text(angle=90,hjust=1,vjust=0.2))
}

#'Score plot
#'@description Provide score plot for multivariate analysis.
#'@usage multiv_scoreplot(METBObj,score_data,pcx=0,pcy=0,oscore_data=NULL,
#'shapeCol=NULL, ptsize=3, plot_title="", legend_title="Color", shape_title="Shape")
#'@param METBObj METBObj object contains list of data.
#'@param score_data a numerical matrix or a data frame of x scores, see details.
#'@param pcx a number indicating variance explained by a principal component (x-axis).
#'@param pcy a number indicating variance explained by a principal component (y-axis).
#'@param oscore_data a numerical matrix or a data frame of orthogonal scores, for OPLS only. Default is \code{NULL}.
#'@param shapeCol a column number/index of category/factor column to show as different point shapes on the plot.
#'If not specified, a filled circle is used by default.
#'@param ptsize a number of geom_point size.
#'@param plot_title text indicating plot title.
#'@param legend_title text indicating the legend title.
#'@param shape_title text indicating the shape title.
#'@details A two-column numerical matrix or a two-column data frame of x scores is required for PCA and PLSDA.
#'@return ggplot object.
#'@examples
#'#adipose_dt = set_input_obj(adipose,2,3,4)
#'#out=multiv_analyze(adipose_dt, method="opls", scale="standard")
#'#multiv_scoreplot(adipose_dt,out$score_val,pcx=out$model_summary$R2X[1],
#'#pcy=out$model_summary$R2X[2],out$oscore_val)
#'@import ggplot2
#'@import RColorBrewer
#'@export
multiv_scoreplot <- function(METBObj, score_data, pcx=0, pcy=0, oscore_data=NULL, shapeCol=NULL, ptsize=3, plot_title="", legend_title="Color", shape_title="Shape"){
  dat = METBObj$X #working data
  if(is.null(oscore_data)){
    plotdata = cbind(Sample_ID=METBObj$ID, data.frame(ClassCol=METBObj$inputdata[,METBObj$classCol], shapeCol=METBObj$inputdata[,shapeCol]), score_data) #for ggplot
  }else{
    plotdata = cbind(Sample_ID=METBObj$ID, data.frame(ClassCol=METBObj$inputdata[,METBObj$classCol], shapeCol=METBObj$inputdata[,shapeCol]), score_data, oscore_data) #orthogonal scores for ggplot
  }
  if(is.null(shapeCol)){
    xlabel=colnames(plotdata)[3]; ylabel=colnames(plotdata)[4] #x,y label
    colnames(plotdata)[3] = 'PCX'; colnames(plotdata)[4] = 'PCY' #change column name
  }else{
    xlabel=colnames(plotdata)[4]; ylabel=colnames(plotdata)[5] #x,y label
    colnames(plotdata)[4] = 'PCX'; colnames(plotdata)[5] = 'PCY' #change column name
  }
  #Check type of category/factor column
  if (is.numeric(plotdata$ClassCol)) {#plot regression
    Y = plotdata$ClassCol #working data
    if(!is.null(shapeCol) && shapeCol != which(colnames(plotdata) == "ClassCol")){# check shapeCol
      pl = ggplot(plotdata, aes(PCX, PCY, group=ClassCol, color = ClassCol, shape=shapeCol, label=Sample_ID))
    }else if(!is.null(shapeCol) && shapeCol == which(colnames(plotdata) == "ClassCol")){# check shapeCol
      pl = ggplot(plotdata, aes(PCX, PCY, group=ClassCol, color = ClassCol, shape=ClassCol, label=Sample_ID))
    }else{
      pl = ggplot(plotdata, aes(PCX, PCY, group=ClassCol, color = ClassCol, label=Sample_ID))
    }
    pl + geom_point(size=ptsize, alpha=0.75) +
      scale_colour_gradientn(colours = rainbow(length(Y), start=0.17, end=1)) +
      theme_minimal() + ggtitle(plot_title) +
      theme(plot.title = element_text(size=12), axis.title=element_text(size=10), axis.text=element_text(size=10)) +
      labs(color=legend_title, shape=shape_title, x = paste(xlabel,"[", round(pcx*100, 2), "%]", sep=""), y = paste(ylabel,"[", round(pcy*100, 2), "%]", sep=""))
  }else{#plot category/factor
    Y = as.factor(plotdata$ClassCol) #working data
    numlevels = nlevels(Y) #no. of levels
    if(numlevels <= 8){
      grcolors = mbcolors
    }else{
      grcolors = colorRampPalette(mbcolors)(numlevels)
    }
    if(!is.null(shapeCol) && shapeCol != which(colnames(plotdata) == "ClassCol")){# check shapeCol
      pl = ggplot(plotdata, aes(PCX, PCY, group=ClassCol, color = ClassCol, shape=shapeCol, label=Sample_ID))
    }else if(!is.null(shapeCol) && shapeCol == which(colnames(plotdata) == "ClassCol")){# check shapeCol
      pl = ggplot(plotdata, aes(PCX, PCY, group=ClassCol, color = ClassCol, shape=ClassCol, label=Sample_ID))
    }else{
      pl = ggplot(plotdata, aes(PCX, PCY, group=ClassCol, color = ClassCol, label=Sample_ID))
    }
    pl + geom_point(size=ptsize, alpha=0.75) +
      stat_ellipse(aes(color=ClassCol), type = "norm", size=0.3) + scale_color_manual(values = grcolors) +
      theme_minimal() + ggtitle(plot_title) +
      theme(plot.title = element_text(size=12), axis.title=element_text(size=10), axis.text=element_text(size=10)) +
      labs(color=legend_title, shape=shape_title, x = paste(xlabel,"[", round(pcx*100, 2), "%]", sep=""), y = paste(ylabel,"[", round(pcy*100, 2), "%]", sep=""))
  }
}

#'VIP plot
#'@description Provide VIP plot for multivariate analysis.
#'@usage multiv_vipplot(vip_data, plot_title="")
#'@param vip_data a numerical vector of VIP.
#'@param plot_title text indicating plot title.
#'@return ggplot object.
#'@examples
#'#adipose_dt = set_input_obj(adipose,2,3,4)
#'#out=multiv_analyze(adipose_dt, method="pls", scale="standard")
#'#multiv_vipplot(out$vip_val)
#'@import ggplot2
#'@import RColorBrewer
#'@export
multiv_vipplot <- function(vip_data, plot_title=""){
  if(!is.vector(vip_data)){
    cat("\nError! vip_data must be a a numerical vector.\nReturn no plot.\n")
    return(FALSE)
  }
  plotdata = data.frame(vname=names(vip_data),value=vip_data) #for ggplot
  plotdata = plotdata[order(plotdata$value, decreasing = TRUE),] #sort
  plotdata$vname = factor(plotdata$vname, levels = plotdata$vname) #set factor
  ggplot(plotdata, aes(x=vname, y=value, fill=value)) + geom_bar(stat = "identity", alpha=0.7) + theme_minimal() +
    scale_fill_gradient(low = "yellow", high = "red") + labs(x = "", y = "VIP", fill = "VIP") + ggtitle(plot_title) +
    theme(plot.title = element_text(size=12), axis.title=element_text(size=10),
          axis.text=element_text(size=10), axis.text.x = element_text(angle=90,hjust=1,vjust=0.2))
}

#'VIP and loading plot
#'@description Provide combined VIP and loading plot for multivariate analysis.
#'@usage multiv_viploadingplot(vip_data, loading_data, oloading_data=NULL, plot_title="")
#'@param vip_data a numerical vector of VIP.
#'@param loading_data a numerical matrix or a data frame of x loadings, see details.
#'@param oloading_data a numerical matrix or a data frame of orthogonal loadings, for OPLS only. Default is \code{NULL}.
#'@param plot_title text indicating plot title.
#'@return ggplot object.
#'@examples
#'#adipose_dt = set_input_obj(adipose,2,3,4)
#'#out=multiv_analyze(adipose_dt, method="pls", scale="standard")
#'#multiv_viploadingplot(out$vip_val,out$loading_data)
#'@import ggplot2
#'@export
multiv_viploadingplot <- function(vip_data, loading_data, oloading_data=NULL, plot_title=""){
  if(!is.vector(vip_data)){
    cat("\nError! vip_data must be a a numerical vector.\nReturn no plot.\n")
    return(FALSE)
  }
  if(!is.data.frame(loading_data)){
    loading_data = data.frame(loading_data)
  }
  if(is.null(oloading_data)){
    plotdata = loading_data #for ggplot
  }else{
    if(!is.data.frame(oloading_data)){
      oloading_data = data.frame(oloading_data)
    }
    plotdata = cbind(loading_data, oloading_data) #orthogonal loadings for ggplot
  }
  plotdata = cbind(rownames(loading_data),plotdata)
  xlabel=colnames(plotdata)[2]; ylabel=colnames(plotdata)[3] #x,y label
  colnames(plotdata)[1] = 'Sample_ID'; colnames(plotdata)[2] = 'PCX'; colnames(plotdata)[3] = 'PCY' #change column name
  vipdata = data.frame(vname=names(vip_data),VIP=vip_data) #for ggplot
  plotdata = merge(plotdata,vipdata,by.x = "Sample_ID",by.y="vname")
  ggplot(plotdata, aes(PCX, PCY, label=Sample_ID)) + geom_point(aes(fill=VIP,size=VIP), alpha=0.5, shape=21) +
    geom_hline(yintercept=0, linetype=1) + geom_vline(xintercept = 0, linetype=1) + theme_minimal() + ggtitle(plot_title) +
    theme(plot.title = element_text(size=12), axis.title=element_text(size=10), axis.text=element_text(size=10)) +
    labs(x = xlabel, y = ylabel)+scale_fill_gradientn(colours=rev(heat.colors(16)))
}

#'PCA plot
#'@description Provide PCA plot overview by category/factor.
#'@usage pcaplot_overview(METBObj, classCol=NULL, shapeCol=NULL, px=1, py=2,
#'scale=FALSE, ptsize=3, plot_title="", legend_title="Color", shape_title="Shape")
#'@param METBObj METBObj object contains list of data.
#'@param classCol a column number/index of category/factor column to show on the plot.
#'If not specified, \code{METBObj$classCol} is used by default.
#'@param shapeCol a column number/index of category/factor column to show as different point shapes on the plot.
#'If not specified, a filled circle is used by default.
#'@param px a number indicating a principal component to be plotted on x-axis (default=1).
#'@param py a number indicating a principal component to be plotted on y-axis (default=2).
#'@param scale a logical value indicating whether the variables should be scaled to have unit variance before the analysis.
#'Default is FALSE.
#'@param ptsize a number of geom_point size.
#'@param plot_title text indicating plot title.
#'@param legend_title text indicating the legend title.
#'@param shape_title text indicating the shape title.
#'@return ggplot object.
#'@examples
#'#pcaplot_overview(METBObj)
#'@import ggplot2
#'@import RColorBrewer
#'@export
pcaplot_overview <- function(METBObj, classCol=NULL, shapeCol=NULL, px=1, py=2, scale=FALSE, ptsize=3, plot_title="", legend_title="Color", shape_title="Shape"){
  dat = METBObj$X #working data
  metadat = METBObj$inputdata[,1:METBObj$xCol-1]
  nmetada = ncol(metadat); xpc = nmetada+1; ypc = nmetada+2; #plot 2 PCs
  out_data = tryCatch({
    pca = prcomp(dat, center = TRUE, scale. = scale)
    pcadata = cbind(metadat, pca$x[,c(px,py)]) #for ggplot
    if(is.null(classCol)){#METBObj$classCol
      colnames(pcadata)[METBObj$classCol] = 'ClassCol' #class column
    }else{#other classCol
      colnames(pcadata)[classCol] = 'ClassCol' #class column
    }
    xlabel=colnames(pcadata)[xpc]; ylabel=colnames(pcadata)[ypc] #x,y label
    colnames(pcadata)[METBObj$idCol] = 'Sample_ID'; colnames(pcadata)[xpc] = 'PCX'; colnames(pcadata)[ypc] = 'PCY' #change column name
    prop.pca = summary(pca)
    #Check type of category/factor column
    if (is.numeric(pcadata$ClassCol)) {#plot regression
      Y = pcadata$ClassCol #working data
      if(!is.null(shapeCol) && shapeCol != which(colnames(pcadata) == "ClassCol")){# check shapeCol
        colnames(pcadata)[shapeCol] = 'shapeCol' #shape column
        pl = ggplot(pcadata, aes(PCX, PCY, group=ClassCol, color = ClassCol, shape=shapeCol, key=Sample_ID))
      }else if(!is.null(shapeCol) && shapeCol == which(colnames(pcadata) == "ClassCol")){# check shapeCol
        pl = ggplot(pcadata, aes(PCX, PCY, group=ClassCol, color = ClassCol, shape=ClassCol, key=Sample_ID))
      }else{
        pl = ggplot(pcadata, aes(PCX, PCY, group=ClassCol, color = ClassCol, key=Sample_ID))
      }
      #geom_text(aes(label=pcadata[,1]),hjust=0.4, vjust=1.3) + #show label
      pl + geom_point(size=ptsize, alpha=0.75) +
        scale_colour_gradientn(colours = rainbow(length(Y), start=0.17, end=1)) + theme_minimal() + ggtitle(plot_title) +
        theme(plot.title = element_text(size=12), axis.title=element_text(size=10), axis.text=element_text(size=10)) +
        labs(col = legend_title, shape=shape_title, x = paste(xlabel, "[", round(prop.pca$importance[2,px]*100, 2), "%]", sep=""),
             y = paste(ylabel, "[", round(prop.pca$importance[2,py]*100, 2), "%]", sep=""))
    }else{#plot category/factor
      Y = as.factor(pcadata$ClassCol) #working data
      numlevels = nlevels(Y) #no. of levels
      if(numlevels <= 8){
        grcolors = mbcolors
      }else{
        grcolors = colorRampPalette(mbcolors)(numlevels)
      }
      if(!is.null(shapeCol) && shapeCol != which(colnames(pcadata) == "ClassCol")){# check shapeCol
        colnames(pcadata)[shapeCol] = 'shapeCol' #shape column
        pl = ggplot(pcadata, aes(PCX, PCY, group=ClassCol, color = ClassCol, shape=shapeCol, key=Sample_ID))
      }else if(!is.null(shapeCol) && shapeCol == which(colnames(pcadata) == "ClassCol")){# check shapeCol
        pl = ggplot(pcadata, aes(PCX, PCY, group=ClassCol, color = ClassCol, shape=ClassCol, key=Sample_ID))
      }else{
        pl = ggplot(pcadata, aes(PCX, PCY, group=ClassCol, color = ClassCol, key=Sample_ID))
      }
      #geom_text(aes(label=pcadata[,1]),hjust=0.4, vjust=1.3) + #show label
      pl + geom_point(size=ptsize, alpha=0.75) + stat_ellipse(aes(color=ClassCol), type = "norm", size=0.3) +
        scale_color_manual(values = grcolors) + theme_minimal() + ggtitle(plot_title) +
        theme(plot.title = element_text(size=12), axis.title=element_text(size=10), axis.text=element_text(size=10)) +
        labs(color=legend_title, shape=shape_title, x = paste(xlabel, "[", round(prop.pca$importance[2,px]*100, 2), "%]", sep=""),
             y = paste(ylabel, "[", round(prop.pca$importance[2,py]*100, 2), "%]", sep=""))
    }
  },
  error=function(e){
    cat(e$message)
    #message(e)
    ggplot()+ggtitle(label=paste("Could not generate a plot:",e$message))+theme(plot.title = element_text(color = "red"))
  })
  return(out_data)
}

#'p-value plot
#'@description Provide p-value plot overview.
#'@usage pvalplot_overview(x_input, cutoff=0.05,
#'plot_title="", ptsize=3)
#'@param x_input data frame of p-values/adjusted p-values. A data frame contains only one column using variables as \code{row.names}.
#'@param cutoff a number indicating statistical significance.
#'@param plot_title text indicating the plot title.
#'@param ptsize a number of geom_point size.
#'@return ggplot object.
#'@examples
#'#sugar_dt = set_input_obj(sugar, 1,2,5)
#'#uni_out=univ_analyze(sugar_dt)
#'#pvalplot_overview(data.frame(uni_out$p_adj),plot_title = "Adjusted p-value")
#'@import ggplot2
#'@import RColorBrewer
#'@export
pvalplot_overview <- function(x_input, cutoff=0.05, plot_title="", ptsize=3){
  if(!is.data.frame(x_input)){
    cat("\nError! x_input must be a data frame.\nReturn no plot.\n")
    return(FALSE)
  }
  dat = data.frame(variables=row.names(x_input),value=x_input) #working data
  colnames(dat)[2] = "value"
  dat$value = -log10(dat$value)
  dat$sig = factor(dat$value > -log10(cutoff))
  ggplot(dat, aes(x=variables, y=value, color=sig, key=variables)) + geom_point(alpha=0.5,size=ptsize) + geom_hline(yintercept=(-log10(cutoff)), linetype=2) +
    scale_color_manual(labels = c("Not significant","Significant"), values = c("#377EB8","#E41A1C")) + theme_light() +
    labs(title = plot_title, x = "", y = "-log10 Scale", color = "") +
    theme(plot.title = element_text(hjust = 0.5, size = 12), axis.title=element_text(size=10),
          axis.text=element_text(size=10), axis.text.x = element_text(angle=90,hjust=1,vjust=0.2))
}

#'RLA plot
#'@description Provide RLA plot overview to visualize sample variation.
#'@usage rlaplot_overview(METBObj, classCol=NULL, type="wg", dolog=TRUE,
#'doavg=FALSE, limitx=FALSE, plot_title="", legend_title="Color")
#'@param METBObj METBObj object contains list of data.
#'@param classCol a column number/index of category/factor column to show on the plot.
#'If not specified, \code{METBObj$classCol} is used by default.
#'@param type text indicating whether within group ("wg") or across group ("ag") to be plotted.
#'@param dolog a logical value indicating whether the variables should be generalized log2-transformed (glog2) prior to plotting.
#'Default is TRUE.
#'@param doavg a logical value indicating whether the standardized variables should be obtained by
#'removing the mean from each variable. Default is FALSE.
#'@param limitx a logical value to limit no. of samples on the plot.
#'@param plot_title text indicating plot title.
#'@param legend_title text indicating legend title.
#'@return ggplot object.
#'@references Rlaplots \url{https://rdrr.io/cran/metabolomics/man/RlaPlots.html}.
#'@examples
#'#sugar_dt = set_input_obj(sugar, 1,2,5)
#'#rlaplot_overview(sugar_dt)
#'@import ggplot2
#'@import RColorBrewer
#'@import dplyr
#'@export
rlaplot_overview <- function(METBObj, classCol=NULL, type="wg", dolog=TRUE, doavg=FALSE, limitx=FALSE, plot_title="", legend_title="Color"){
  if(dolog){
    dat = transform_input_data(METBObj,method="glog2")$X #working data
  }else{
    dat = METBObj$X #working data
  }
  metadat = METBObj$inputdata[,1:METBObj$xCol-1]
  if(is.null(classCol)){#METBObj$classCol
    inputdata = cbind(ClassCol=metadat[,METBObj$classCol], dat) #for ggplot
  }else{#other classCol
    inputdata = cbind(ClassCol=metadat[,classCol], dat) #for ggplot
  }
  if(is.numeric(inputdata$ClassCol)){type = "ag"}
  #within groups (type == "wg")
  if(type == "wg") {
    if(doavg){
      out_data = inputdata %>%
        group_by(ClassCol) %>%
        mutate_at(vars(-c("ClassCol")), list(~ .-mean(.)))
    }else{
      out_data = inputdata %>%
        group_by(ClassCol) %>%
        mutate_at(vars(-c("ClassCol")), list(~ .-median(.)))
    }
  #across groups (type == "ag")
  }else{
    if(doavg){
      out_data = inputdata %>% mutate_at(vars(-c("ClassCol")), list(~ .-mean(.)))
    }else{
      out_data = inputdata %>% mutate_at(vars(-c("ClassCol")), list(~ .-median(.)))
    }
  }
  out_data = out_data[order(out_data$ClassCol),] #order by class
  out_data$ind = factor(row.names(out_data), levels = row.names(out_data))
  dim_dt = dim(out_data)
  plot_data = reshape2::melt(out_data, id.vars=c("ind","ClassCol"))
  if (is.numeric(inputdata$ClassCol)) {#plot regression
    pl = ggplot(plot_data, aes(x = ind, y = value))
    if(limitx && dim_dt[1] > 100){
      pl = pl + scale_x_discrete(limits = factor(c(1:100))) #xlimit 1-100
    }
    pl + geom_boxplot(outlier.alpha = 0.5) + theme_light() + ggtitle(plot_title) + labs(x = "Sample", y="") +
      theme(plot.title = element_text(size = 12), axis.text=element_text(size=10), axis.text.y = element_blank()) + coord_flip() #by samples
  }else{#plot category/factor
    Y = as.factor(inputdata$ClassCol) #working data
    numlevels = nlevels(Y) #no. of levels
    if(numlevels <= 8){
      grcolors = mbcolors
    }else{
      grcolors = colorRampPalette(mbcolors)(numlevels)
    }
    pl = ggplot(plot_data, aes(x = ind, y = value, color=ClassCol))
    if(limitx && dim_dt[1] > 100){
      pl = pl + scale_x_discrete(limits = factor(c(1:100))) #xlimit 1-100
    }
    pl + geom_boxplot(outlier.alpha = 0.5)+
      scale_color_manual(values = grcolors) + theme_light() + ggtitle(plot_title) + labs(x = "Sample", y="", color=legend_title) +
      theme(plot.title = element_text(size = 12), axis.text=element_text(size=10), axis.text.y = element_blank()) + coord_flip() #by samples
  }
}
