#' PatientLog
#'
#' @param Metadata A metadata object
#' @import readr
#' @return logs
#' @export
#'
#' @examples 'non'
PatientLog = function(Metadata){


  dt = as.data.frame(na.omit(suppressWarnings(readr::read_csv(file.path(Processpath(Metadata),"Patients.CleanedProcess.txt"),
           skip = 1,show_col_types = F))))



  return(dt)
}

#' SampleLog
#'
#' @param Metadata A metadata object
#' @import readr
#' @return logs
#' @export
#'
#' @examples 'non'
SampleLog = function(Metadata){


  dt = as.data.frame(na.omit(suppressWarnings(readr::read_csv(file.path(Processpath(Metadata),"Samples.CleanedProcess.txt"),
                                                              skip = 1,show_col_types = F))))


  return(dt)
}

#' CellLog
#'
#' @param Metadata A metadata object
#' @import readr
#' @return logs
#' @export
#'
#' @examples 'non'
CellLog = function(Metadata){


  dt = as.data.frame(na.omit(suppressWarnings(readr::read_csv(file.path(Processpath(Metadata),"CellAnnotation.CleanedProcess.txt"),
                                                              skip = 1,show_col_types = F))))



  return(dt)
}
