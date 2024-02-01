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


  if(key%in%names(lexic)){

    if(!value%in%lexic[[key]]){
      message(paste("Adding",value, "in", key, "."))
    lexic[[key]] <-  c(lexic[[key]],value)
    } else { message(paste(value, "in", key, "already present."))}



  } else {

    lexic[[key]] <- c(key,value)


    }

  return(lexic)






}
