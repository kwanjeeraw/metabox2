#'Run MUVR
#'@description Perform ML analysis using MUVR package.
#'@usage run_muvr(METBObj, scale=FALSE, nRep=3, varRatio=0.75, partitionSize=0.7,
#'method="PLS", DA=TRUE, fitness="AUROC")
#'@param METBObj METBObj object contains list of data.
#'@param scale a logical variable indicating whether the variables should be scaled to have unit variance before PLS modelling. Default is FALSE.
#'@param nRep number of repetitions of double CV.
#'@param varRatio ratio of variables to include in inner loop iterations. Default is 0.75.
#'@param partitionSize number of data holding for training set. Default is 0.7.
#'@param method ML method. Choose one from the list: PLS, RF. Default is PLS.
#'@param DA a logical variable for performing classification. Default is TRUE. For performing regression, DA=FALSE.
#'@param fitness fitness function for model tuning and evaluation. Choose AUROC (default) or MISS for classification; or choose RMSEP for regression.
#'@details Based on MUVR version 0.0.976. (2023-03-09). The \code{partitionSize} is for calculating \code{nOuter}.
#'@return MVObject object.
#'@author Kwanjeera W \email{kwanjeera.wan@@mahidol.ac.th}
#'@references Lin Shi, et al. MUVR (2019) \url{https://doi.org/10.1093/bioinformatics/bty710}.
#'@seealso \code{\link[MUVR:MUVR]{MUVR::MUVR()}}, \code{\link[MUVR:vectSamp]{MUVR::vectSamp()}}
#'@examples
#'#res_muvr = run_muvr(METBObj, method="PLS")
#'#MUVR::plotVAL(res_muvr) #need more distinct colors
#'#MUVR::plotMV(res_muvr, model = 'min')
#'#MUVR::plotStability(res_muvr, model = 'min')
#'#MUVR::plotVIP(res_muvr, model = 'min') #error if nrep=1
#'#MUVR::getVIP(res_muvr, model = 'min')
#'@export
run_muvr <- function(METBObj, scale=FALSE, nRep=3, varRatio=0.75, partitionSize=0.7, method="PLS", DA=TRUE, fitness="AUROC"){
  #Check argument
  tmparg_method <- try(method <- match.arg(method, c("PLS","RF"), several.ok = FALSE), silent = TRUE)
  if (class(tmparg_method) == "try-error") {
    cat("\nERROR! Argument 'method' is not valid, choose one from the list: PLS,RF.\n")
    return(list())
  }
  tmparg_fitness <- try(fitness <- match.arg(fitness, c("AUROC","MISS","RMSEP"), several.ok = FALSE), silent = TRUE)
  if (class(tmparg_fitness) == "try-error") {
    cat("\nERROR! Argument 'method' is not valid, choose one from the list: AUROC,MISS,RMSEP.\n")
    return(list())
  }
  #Check argument
  if(sum(is.na(METBObj$X)) > 0){#Data contains missing values
    cat("\nThe data contains missing values. Data was not analyzed.\n")
    return(FALSE)
  }

  ###Parallel processing
  # library(foreach)
  # if (parallel) "%doVersion%"=get("%dopar%") else "%doVersion%"=get("%do%")


  ###INITIALIZE PARAMETERS
  # scale=scale
  # modReturn=TRUE
  # methParam=list(compMax=5,robust=0.05,ntreeIn=150,ntreeOut=300,mtryMaxIn=150,oneHot=FALSE,NZV=FALSE,rfMethod="randomForest")
  # #methParams
  # if (method=='PLS') {
  #   methParam$NZV <- TRUE
  # }
  # ML=FALSE
  # parallel=FALSE
  # if (method == 'RF') library(randomForest)
  # library(pROC)
  # #files.sources = list.files("MUVR",pattern="*.R", full.names=TRUE)
  # #sapply(files.sources, source)
  #
  # ###INITIALIZE VARIABLES
  # start.time=proc.time()[3]
  # modelReturn=list(call=match.call())
  ID = METBObj$ID #ID
  orgID = METBObj$orgID #orgID
  unik = METBObj$unik #unik
  unikID = METBObj$unikID #unikID
  classCol = METBObj$classCol #classCol
  classLs = METBObj$Y #classLs-list of classes for DA only
  X = METBObj$X #X
  Y = METBObj$Y #Y

  ###START CALCULATION
  #Find nOuter
  if (DA & identical(unikID,ID)) {##Independent samples => Classification
    cat('\nClassification; Assuming samples are independent.')
    nTest = (table(classLs))-(round(table(classLs) * partitionSize)) #finding test size
    nOuter = min(round((table(classLs))/nTest)) #finding no. of folds
  } else{##Repeated samples => Regression|Classification; Independent samples => Regression
    cat('\nClassification or Regression; Assuming samples are dependent.')
    nTest = (length(unikID))-(round(length(unikID) * partitionSize)) #finding test size
    nOuter = round(length(unikID)/nTest) #finding no. of folds
  }
  if(nOuter <= 1){
    nOuter = 2 #default no. of folds
    nInner = nOuter #no. of inner folds
  }else{
    nInner = nOuter-1 #default no. of inner folds
  }
  returnMuvr = MUVR::MUVR(X=X, Y=Y, scale=scale, nRep=nRep, nOuter=nOuter, nInner=nInner,
                          varRatio=varRatio, DA=DA, fitness = fitness, method=method)
  return(returnMuvr)
}
