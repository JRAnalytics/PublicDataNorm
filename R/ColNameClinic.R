
#' ColNameClinic function to prepare constant colnames to find in clinical data from lexiue .txt
#'
#' @param list.files.path file path to find lexique of colnames
#'
#' @return a list
#' @export
#'
#' @examples
ColNameClinic <- function(list.files.path){

cnameClinic <- data.table::fread(paste(list.files.path$Parent,"References","colnames clinical data.txt",sep = "/"))


cnameClinic <- as.data.frame(cnameClinic)

l <- list()
for (i in 1:ncol(cnameClinic)) {
  nl <- names(l)
  l <- c(l, list(na.omit(cnameClinic[,i])))
  names(l) <-  c(nl, colnames(cnameClinic)[i])

}

return(l)
}






