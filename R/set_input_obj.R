#'Set input data object
#'@description Set input data object
#'@usage set_input_obj(inputdata, idCol, classCol, xCol)
#'@param inputdata dataframe of input data. Samples are in rows and variables/features are in columns.
#'@param idCol a column number/index of sample ID column which indicates independent or repeated sample structure of the experimental design. Default is 1.
#'@param classCol a column number/index of category/factor column (for hypothesis testing and classification)
#'or response column (for regression). Default is 2.
#'@param xCol a column number/index of the first variable/feature column. Default is 3.
#'@return METBObj object; a list of the following components:
#'
#'inputdata = data frame of input data.
#'
#'ID = vector of all ids.
#'
#'orgID = vector of original ids.
#'
#'unik = vector of logical variables indicating unique ids, for MUVR.
#'
#'unikID = vector of unique ids, for MUVR.
#'
#'idCol = a column number/index indicates independent or repeated sample structure of the experimental design.
#'
#'xCol = a column number/index of the first variable/feature column in the inputdata data frame.
#'
#'classCol = a column number/index of category/factor column (for hypothesis testing and classification)
#'or response column (for regression) in the inputdata data frame.
#'
#'X = data frame of variables/features.
#'
#'Y = vector of factors (for classification) or response variables (for regression).
#'
#'isRepeated = a logical variable indicating a repeated-measures design.
#'
#'details = a list of object details.
#'
#'@author Kwanjeera W \email{kwanjeera.wan@@mahidol.ac.th}
#'@examples
#'#out = set_input_obj(read_input_data('filename'),1,2,3)
#'#out = set_input_obj(adipose, 2,3,4)
#'@export
set_input_obj <- function(inputdata, idCol=1, classCol=2, xCol=3){
  #Initialize parameters
  input_obj = list(); #Working data
  #Set ID; ID is numeric; for MUVR
  if (missing(idCol)) {
    cat("\nMissing sample ID column; Assuming samples are independent.\n")
    ID = 1:nrow(inputdata) #list of samples
    orgID = ID #keep original IDs, for statistical analyses
  }else{
    ID = inputdata[,idCol] #list of ids
    orgID = ID #keep original IDs, for statistical analyses
    if (is.character(ID)) {
      ID=factor(ID, levels = ID)
    }else if(is.factor(ID)) {
      ID=ID
    }else {
      ID=as.numeric(ID) #convert IDs to numbers, for MUVR
    }
  }
  # #Set classCol; classCol is numeric
  # if (missing(classCol)) {
  #   cat("\nMissing class/response column; Assuming it is next to sample ID column.\n")
  #   classCol = max(which(sapply(inputdata, is.factor)))
  # }
  # #Set xCol; xCol is numeric
  # if (missing(xCol)) {
  #   cat("\nMissing the first variable/feature column; Assuming it is next to class/response column.\n")
  #   xCol = 3
  # }

  ###INITIALIZE VARIABLES
  X = inputdata[,-1:-(xCol-1)] #X;
  if(!is.data.frame(X)){ #if 1 variable;
    X=data.frame(matrix(inputdata[,-1:-(xCol-1)],ncol=1)) #X;
    colnames(X) = colnames(inputdata)[xCol]
  }
  Y = inputdata[,classCol] #Y;
  unik = !duplicated(ID)  #boolean of unique IDs; for MUVR
  unikID = ID[unik] #list of ids; for MUVR
  isRepeated = (!identical(unikID,ID)) #boolean of repeated experiment

  cat("\nUploaded data contains:\n-",
      nrow(X), "samples and", ncol(X), "variables.\n-",
      "repeated samples =",isRepeated,".\n-",
      "class/factor with",nlevels(Y), "levels.\n-",
      sum(apply(X, 2, function(x) { sum(x<0) }),na.rm = TRUE)," negative variables.\n-",
      sum(is.na(X)), "(",round((sum(is.na(X))/(nrow(X)*ncol(X)))*100,2),"%) missing values.\n")

  input_obj$inputdata = inputdata; input_obj$ID = ID;  input_obj$orgID = orgID; input_obj$unik = unik; input_obj$unikID = unikID; input_obj$idCol = idCol;
  input_obj$xCol = xCol; input_obj$classCol = classCol; input_obj$X = X; input_obj$Y = Y; input_obj$isRepeated = isRepeated; input_obj$details = list();
  class(input_obj)=c('METBObj')
  return(input_obj)
}
