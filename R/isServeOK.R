#' Check GDC TCGA server statu
#'
#' @return true of false
#' @export
#' @import TCGAbiolinks
#' @examples isServeOK()
#'
isServeOK <- function(){
  tryCatch({
    status <- TCGAbiolinks::getGDCInfo()$status
    if(status != "OK") stop("GDC server down, try to use this package later")
  },error = function(e) stop("GDC server down, try to use this package later"))
  return(TRUE)
}
