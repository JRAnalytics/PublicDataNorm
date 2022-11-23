#' ColNameClinic function to prepare constant colnames to find in clinical data from lexiue .txt
#'
#' @param Lexical_colnames_path file path to find lexique of colnames
#' @importFrom data.table fread
#' @importFrom stats na.omit
#' @return a list
#' @export
#'
#' @examples "none"
ColNameClinic <- function(Lexical_colnames_path){

cnameClinic <- data.table::fread(paste(Lexical_colnames_path,"Lexic.txt",sep = "/"))


cnameClinic <- as.data.frame(cnameClinic)

l <- list()
for (i in 1:ncol(cnameClinic)) {
  nl <- names(l)
  l <- c(l, list(na.omit(cnameClinic[,i])))
  names(l) <-  c(nl, colnames(cnameClinic)[i])

}

return(l)
}






