#' Extract attributes from a gtf files to extract for example gene names ou type
#'
#' @param gtf_attributes column attribute in the GTF file
#' @param att_of_interest a character specifying the targeted atribute to extract
#'
#' @return attribute of interest from the targeted line of data.frame with attributes
#' @import data.table magrittr
#'
#'
#' @examples
#' "non"
#' @export
 extract_attributes <- function(gtf_attributes, att_of_interest){
  if(is.null(gtf_attributes)) stop(paste("GFT$attributes is missing"))
  if(is.null(att_of_interest)) stop(paste("No targeted attribute is mentionned"))

  att <- unlist(strsplit(gtf_attributes, " "))
  if(att_of_interest %in% att){
    return(gsub("\"|;","", att[which(att %in% att_of_interest)+1]))
  }else{
    return(NA)}
}


