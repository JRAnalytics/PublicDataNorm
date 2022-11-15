#' SubsetMetaData subsettinf Metadata with samples of interests.
#'
#' @param Metadata meta object
#' @param Samples.selection string of characters samples
#'
#' @return a meta object
#' @export
#'
#' @examples "none"
SubsetMetaData <- function(Metadata, Samples.selection= NULL){

  lm <- length(Metadata)

  if (!is.null(Samples.selection)){


      for (i in 1:lm){

        if(str_detect(names(Metadata[i]), "clinic|pheno")==T) {Metadata[[i]] <- Metadata[[i]][Samples.selection,] } else if(str_detect(names(Metadata[i]), "gene")==T) {}

          else {Metadata[[i]] <- Metadata[[i]][,Samples.selection] }
             }
  }

  else {
    for (i in 2:lm){

      if(str_detect(names(Metadata[i]), "clinic|pheno")==T) { Metadata[[i]] <- Metadata[[i]][colnames(Metadata[[1]]),]}  else if(str_detect(names(Metadata[i]), "gene")==T) {}
      else  { Metadata[[i]] <- Metadata[[i]][,colnames(Metadata[[1]])]}


    }
    }

  return(Metadata)





}
