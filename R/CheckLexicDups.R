#' CheckLexicDups
#'
#' @param Lexic Lexic to check if duplicated values (Targeted column) are in the lexic items.
#'
#' @return nothing
#' @export
#'
#' @examples "non"
CheckLexicDups = function(Lexic) {

dup = which(duplicated(unlist(lapply(Lexic, function(x){x = x[-1]}))))
columnT = unlist(lapply(Lexic, function(x){x = x[-1]}))[dup]

ItemsN = names(unlist(Lexic))[which(unlist(Lexic) %in% columnT)]

Items = unique(names(Lexic)[unlist(lapply(ItemsN, function(x) which(str_detect(x, names(Lexic)))))])



if(length(dup)>0){message(c(paste0("Duplicated column target, '",columnT,"' in Lexic items: "),paste0(Items,sep = ", ")))}


}



