#' AddKeyLexic
#'
#' @param lexic PatientLexic or SamplesLexic
#' @param key a character string: name to add into lexic
#' @param value a character string: target column in "raw" clinical data to rename as 'key' in cleaned clinical data.
#' @return a Lexic
#' @export
#'
#' @examples "None"
AddKeyLexic <- function(lexic=NULL, key = NULL,value= NULL){



  if(value%in%unlist(lexic)){


    ItemsN = names(unlist(lexic))[which(unlist(lexic) %in% value)]
    Items = unique(names(lexic)[unlist(lapply(ItemsN, function(x) which(str_detect(x, names(lexic)))))])


    message(c(paste0("Value '", value, "' already present in: "), paste0(Items, sep= ", "), "will not be added.
              Suppress Value from Key items.
              Try:  Lexic[[key]] = Lexic[[key]][Lexic[[key]]!= c(value)] "))} else{


    if(key%in%names(lexic)){

    if(!value%in%lexic[[key]]){
      message(paste("Adding",value, "in", key, "."))
    lexic[[key]] <-  c(lexic[[key]],value)
    } else { }



  } else {

    lexic[[key]] <- c(key,value)


    }}

  return(lexic)


}



#' addSeveralKeysToLexic
#'
#' @param vector a list of c(key=value) for
#' @param lexic lexic to add
#'
#' @return a lexic
#' @export
#'
#' @examples "not"

addSeveralKeysToLexic <- function(vector=NULL, lexic=NULL){

  if(is.null(vector)){stop("vector is null")}
  if(is.null(lexic)){stop("lexic is null")}
  if(is.null(names(vector))){stop("vector has no names. vector must be c('key' = 'value') to add in lexic. see ?AddKeyLexic")}


  for (i in 1:length(vector)) {
    lexic = AddKeyLexic(lexic = lexic,
                        key = names(vector)[i],
                        value =vector[[i]])
  }
  return(lexic)
}

