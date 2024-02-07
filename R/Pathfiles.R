#' Rawpath Export Dataset inside object into ".csv" files
#'
#' @param Dataset a Dataset  data files
#' @return a path to 01RawData project file
#' @export
#' @import utils
#' @import R.utils
#' @import Matrix
#' @examples "non"
Rawpath <- function (Dataset){
  return(attributes(Dataset)$File.path$Project.RawData)
}


#' Processpath Export Dataset inside object into ".csv" files
#'
#' @param Dataset a Dataset  data files
#' @return a path to 02process project file
#' @export
#' @import utils
#' @import R.utils
#' @import Matrix
#' @examples "non"
Processpath <- function (Dataset){

  return(attributes(Dataset)$File.path$Project.Processes)

}

#' Refpath Export Dataset inside object into ".csv" files
#'
#' @param Dataset a Dataset  data files
#' @return a path to 03Referecence project file
#' @export
#' @import utils
#' @import R.utils
#' @import Matrix
#' @examples "non"
Refpath <- function (Dataset){
  return(attributes(Dataset)$File.path$References)
}

#' Verifiedpath Export Dataset inside object into ".csv" files
#'
#' @param Dataset a Dataset  data files
#' @return a path to 04VerifiedDataSet project file
#' @export
#' @import utils
#' @import R.utils
#' @import Matrix
#' @examples "non"
Verifiedpath <- function (Dataset){
  return(attributes(Dataset)$File.path$Project.VerifiedDataset)

}

