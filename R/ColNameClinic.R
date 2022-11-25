#' ColNameClinic function to prepare constant colnames to find in clinical data from lexiue .txt
#'
#' @param Lexical_colnames_path file path to find lexique of colnames
#' @importFrom data.table fread
#' @import dplyr
#' @importFrom stats na.omit
#' @return a list
#' @export
#'
#' @examples "none"
ColNameClinic <- function(Lexical_colnames_path){



x <- scan(paste(Lexical_colnames_path,"Lexic.txt",sep = "/"), what="", sep="\n")%>%strsplit("[[:space:]]+")

names(x) <- sapply(x, `[[`, 1)
x<- lapply(x, `[`, -1)

return(x)
}






