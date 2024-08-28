#'Perform integrative analysis
#'@description perform integrative analysis of several data sets using a multi-block partial least squares discriminant analysis (MBPLSDA). The function wraps around the main functions of \pkg{\link{packMBPLSDA}} and \pkg{\link{ade4}}.
#'@usage mbplsda_analyze(class_data, input_data, scale=TRUE, option="none", nf=10, optdim=2,
#'nrepet=30, npermut=100, nboot=100, threshold=0.5, testmodel=FALSE, cpus=1)
#'@param class_data one-column data frame of class/category/factor data.
#'@param input_data a list of data frames containing 2 or more data sets of explanatory/independent variables. The data frames must have the same samples and unduplicated variables.
#'@param scale a logical value to standardize the explanatory variables
#'@param option a string specifying the option for the block weighting. It can be one of none (default), uniform. If none, the block weight is equal to the block inertia.
#'If uniform, the weight of each explanatory data set is equal to 1/number of explanatory data sets, and the weight of the Y-block is equal to 1.
#'@param nf a number of principle components to be calculated. Default is 10.
#'@param optdim an optimal number of principle components included in the model. Default is 2.
#'@param nrepet a number of cross-validation repetitions. See \code{\link{testdim_mbplsda}} and \code{\link{permut_mbplsda}}.
#'@param npermut a number of permutations. See \code{\link{permut_mbplsda}}.
#'@param nboot a number of bootstrap replications. See \code{\link{boot_mbplsda}}.
#'@param threshold a number indicating the prediction threshold, between 0 and 1 (Default is 0.5). See \code{\link{testdim_mbplsda}}.
#'@param testmodel a logical variable to perform model components testing and permutation. Default is FALSE.
#'@param cpus a number of CPU for parallel computing.
#'@details The function automatically performs all steps of MBPLSDA including generating base MBPLSDA model, determination of optimal components based on Error Rates (ER),
#'permutation test for model validity and parameter estimation.
#'@return a list of the following components:
#'
#'result = MBPLSDA results.
#'
#'details = a list of analysis details: testMethod, nclass, inputsize, cpus, testmodel, ncomp, noptcomp, call.
#'
#'@author Kwanjeera W \email{kwanjeera.wan@@mahidol.ac.th}
#'@references Brandolini-Bunlon M., et al. (2019) Multi-block PLS discriminant analysis. Metabolomics, 15(10):134.
#'@references Bougeard S. and Dray S. (2018) Supervised Multiblock Analysis in R with the ade4 Package.Journal of Statistical Software,86(1), 1-17.
#'@seealso \code{\link[packMBPLSDA:mbplsda]{packMBPLSDA::mbplsda()}}, \pkg{\link{packMBPLSDA}}, \pkg{\link{ade4}}
#'@examples
#'#sugar_dt = set_input_obj(sugar, 1,2,5)
#'#out = mbplsda_analyze(sugar_dt$inputdata[,2:3], list(omics = sugar_dt$X),
#'#nrepet=10, npermut=10, nboot=10)
#'#packMBPLSDA::plot_boot_mbplsda(out$result$res_boot)
#'#### Adipose example
#'#gc_dt = set_input_obj(read_input_data("data/adipose_GC.csv"),1,2,3)
#'#lc_pos = set_input_obj(read_input_data("data/adipose_LC_POS.csv"),1,2,3)
#'#lc_imp = impute_missing_data(lc_pos,method = "bpca",cutoff=30)
#'#lc_neg = set_input_obj(read_input_data("data/adipose_LC_NEG.csv"),1,2,3)
#'#out = mbplsda_analyze(data.frame(gc_dt$inputdata[,2]), list(gc = gc_dt$X,lcpos=lc_imp$X, lcneg=lc_neg$X),
#'#nrepet=10, npermut=3, nboot=3, testmodel = F,cpus = 2)
#' @export
mbplsda_analyze <- function(class_data, input_data, scale=TRUE, option="none", nf=10, optdim=2, nrepet=30, npermut=100, nboot=100, threshold=0.5, testmodel=FALSE, cpus=1){
  #Check argument
  if (class(class_data) != "data.frame"){
    stop("argument 'class_data' is not valid, a data frame is required.")
  }
  if (ncol(class_data) > 1){
    cat("\nargument 'class_data' contains more than one column, only the 1st column will be used.")
  }
  if (class(input_data) != "list"){
    stop("argument 'input_data' is not valid, a list is required.")
  }
  if (length(input_data) < 2){
    #stop("argument 'input_data' is not valid, two or more data sets are required.")
  }
  #Initialize parameters
  mbplsda_result = list(); methodls = list(); #Working data
  #Check category/factor column
  if (is.factor(class_data[,1])) {
    F1 = data.frame(class_data[,1]) #working data
    colnames(F1) = colnames(class_data)[1]
    cat('\nFind a factor F1 of ',nlevels(class_data[,1]),' levels.',sep='')
  }else if(is.numeric(class_data[,1])){
    F1 = data.frame(class_data[,1]) #working data
    colnames(F1) = colnames(class_data)[1]
    cat('\nFind a factor F1 of ',length(unique(class_data[,1])),' levels.',sep='')
  }else{
    F1 = data.frame(as.factor(class_data[,1])) #working data
    colnames(F1) = colnames(class_data)[1]
    cat('\nThe category/factor column is', class(class_data[,1]), '.\nConverting to a factor F1 of ',nlevels(F1[,1]),' levels.',sep='')
  }
  if(nf < 0 || nf > 10){maxCmp = 10}else{maxCmp = nf} #max no. of scanned components
  if(optdim > nf){maxOpt = maxCmp}else if(optdim < 0 || optdim > 10){maxOpt = maxCmp}else{maxOpt = optdim} #optimal no. of components

  cat("\nExecuting function ...\n")
  inst_pkg = NULL
  if(!requireNamespace("ade4", quietly = TRUE)){#check and install required package
    cat("\nMissing the required package 'ade4', trying to install the package ...\n")
    inst_pkg = install_pkgs('ade4')
  }
  if(!requireNamespace("packMBPLSDA", quietly = TRUE)){#check and install required package
    cat("\nMissing the required package 'packMBPLSDA', trying to install the package ...\n")
    inst_pkg = install_pkgs('packMBPLSDA')
  }
  if(length(unlist(inst_pkg))){
    mbplsda_res = list()
    cat("\nERROR! Could not install the required packages. Data was not analyzed.\n")
  }else{
    require(ade4, quietly=TRUE); require(packMBPLSDA, quietly=TRUE);
    base_mb = tryCatch({
      cat("Calculating base model ...\n")
      ls_X = ktab.list.df(input_data)
      dj_table = disjunctive(F1)
      dudi_pca = dudi.pca(dj_table , center = FALSE, scale = FALSE, scannf = FALSE)
      nlev = ncol(dj_table)
      mbplsda(dudi_pca, ls_X, scale = scale, option = option, scannf = FALSE, nf = maxCmp)
    },error=function(e){
      cat(e$message)
      #message(e)
      cat("\nERROR! Data was not analyzed.\n")
      list()
    })
    if(testmodel){
      test_dim = tryCatch({
        cat("Testing model components, this process might take long time ...\n")
        testdim_mbplsda(object = base_mb, nrepet = nrepet, threshold = threshold, bloY = nlev, cpus = cpus, algo = c("max"), outputs = c("ER"))
      },error=function(e){
        cat(e$message)
        #message(e)
        cat("\nERROR! Data was not analyzed.\n")
        list()
      })
      test_perm = tryCatch({
        cat("Performing permutation testing, this process might take long time ...\n")
        permut_mbplsda(base_mb, nrepet = nrepet, npermut = npermut, optdim = maxOpt, bloY = nlev, nbObsPermut = 10, cpus = cpus, algo = c("max"), outputs = c("ER"))
      },error=function(e){
        cat(e$message)
        #message(e)
        cat("\nERROR! Data was not analyzed.\n")
        list()
      })
    }else{
      test_dim=list();test_perm=list();
    }
    test_boot = tryCatch({
      cat("Performing bootstrapping ...\n")
      boot_mbplsda(base_mb, optdim = maxOpt, nrepet = nboot, cpus=cpus)
    },error=function(e){
      cat(e$message)
      #message(e)
      cat("\nERROR! Data was not analyzed.\n")
      list()
    })
    mbplsda_res = list(base_model = base_mb, res_optimal = test_dim, res_permut = test_perm, res_boot = test_boot)
    methodls$testMethod = "Multi-block partial least squares discriminant analysis"; methodls$nclass=nlevels(F1[,1]); methodls$inputsize=length(input_data);
    methodls$cpus = cpus; methodls$testmodel = testmodel;
    methodls$ncomp=mbplsda_res$base_model$nf; methodls$noptcomp=maxOpt;  methodls$call=match.call();
  }
  mbplsda_result$result = mbplsda_res; mbplsda_result$details = methodls;
  unloadNamespace("packMBPLSDA"); unloadNamespace("ade4")
  return(mbplsda_result)
}
