#' AddLexic
#'
#' @param lexic PatientLexic or SamplesLexic
#' @param Param string of character c("A","B", etc) : "A" names of parameter, will be a colnames in cleaned clinic. "B"& etc targeted colname(s) from raw clinicil data
#'
#' @return a Lexic
#' @export
#'
#' @examples None
AddKeyLexic <- function(lexic, Param = c("IDname", "param target")){


  if(Param[1]%in%names(lexic)){

    if(!Param[-1]%in%lexic[[Param[1]]]){
      message(paste("Adding",Param[-1], "in", Param[1], "."))
    lexic[[Param[1]]] <-  c(lexic[[Param[1]]],Param[-1])
    } else { message(paste(Param[-1], "in", Param[1], "already present."))}



  } else {

    lexic[Param[1]] <- Param


    }

  return(lexic)






}
