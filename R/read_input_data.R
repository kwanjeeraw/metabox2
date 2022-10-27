#'Read CSV file
#'@description Read CSV file
#'@usage read_input_data(filename)
#'@param filename the file path to the CSV file.
#'The input file contains data, which samples are in rows and variables/features are in columns separated by comma (,).
#'Information such as batch, replication, sample type and injection order, each must be in a separate column.
#'@return data frame; row: variables/features, col: samples
#'@author Kwanjeera W \email{kwanjeera.wan@@mahidol.ac.th}
#'@examples
#'#out = read_input_data('filename')
#'#fat = set_input_obj(read_input_data('data/repeated_samples_classification.csv'),2,3,4)
#'#lung = set_input_obj(read_input_data('data/independent_samples_classification.csv'),,5,6)
#'#freelive = set_input_obj(read_input_data('data/regression_example.csv'),1,2,3)
#'@export
read_input_data <- function(filename){
  dat <- tryCatch(
    {
      read.csv(file = filename, stringsAsFactors=TRUE, check.names=FALSE)
    }, error=function(e){
      print(e)
      data.frame() #return value in case of error
    })

  if(any(dim(dat) == 0) || ncol(dat) == 1){
    cat("\nERROR! Fail to read the input data.\n"
                 ,"Check the following guidelines:\n"
                 ,"- the input data is a table of comma separated values with .csv file format.\n"
                 ,"- samples are in rows and variables/features are in columns.\n"
                 ,"- missing values should be blank or NA.\n")
    data.frame() #return value in case of error
  }
  cat("\nUploaded data with", nrow(dat), "rows and", ncol(dat), "columns.\n")
  return(dat)
}
