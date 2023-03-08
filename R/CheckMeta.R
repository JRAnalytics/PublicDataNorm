#' CheckMeta : Checking samples accors MetaData files
#'
#' @param MetaData a mMetaData list
#'
#' @return none
#' @export
#'
#' @examples "none"
CheckMeta <- function(MetaData) {

  if(is.null(MetaData)){stop("Need a MetaData List file")}
  if(!is.list(MetaData)){stop("Need a MetaData List file")}


  l <-length(names(MetaData))
  MetaDataN <- names(MetaData)

  ColN <- colnames(MetaData[[1]])
  m <- which(attributes(Meta)$Data.Type=="Expression.Matrix")
  c <- which(attributes(Meta)$Data.Type=="Clinical.data")

  message("-------------------------")
  message("Checking colnames of matrices in MetaData from Raw gene expression matrix")
  message("-------------------------")
  for (i in m){



    if(all(ColN %in% colnames(MetaData[[i]]))==T) {   message(paste(MetaDataN[i]), " colnames : PASS") } else {  message(paste(MetaDataN[i]), " colnames : FAIL") }

  }
  message("-------------------------")
  message("Checking Matrices probes")
  message("-------------------------")


    for (i in m){



      if(all(str_detect(rownames(MetaData[[i]]), "ENSG")==T)) { message(paste(names(MetaData[i]), " gene probes as ENSEMBL"))

        suma <- summary(rownames(MetaData[[i]])%in%MetaData$geneAnnotation$EnsemblID)
        names(suma) <- c("Mode", "Gene not found", "Found")
        print(suma)
      }

      if(all(rownames(MetaData[[i]])%in%MetaData$geneAnnotation$Entrez.id)==T)  { message(paste(names(MetaData[i]), " gene probes as ENTREZ gene id"))

        suma <- summary(rownames(MetaData[[i]])%in%MetaData$geneAnnotation$Entrez.id)
        names(suma) <- c("Mode", "Gene not found", "Found")
        print(suma)
      }


       if(all(str_detect(rownames(MetaData[[i]]), "ILMN_")==T)) { message(paste(names(MetaData[i]), " gene probes as Illumina Bead Array Probes"))

        suma <- summary(rownames(MetaData[[i]])%in%rownames(MetaData$geneAnnotation))
        print(suma)
      }


      if(all(str_detect(rownames(MetaData[[i]]), "_at")==T)) { message(paste(names(MetaData[i]), " gene probes as Illumina Microarray Probes"))

        suma <- summary(rownames(MetaData[[i]])%in%rownames(MetaData$geneAnnotation))
        print(suma)
      }



      if(length(which(rownames(MetaData[[i]])%in%MetaData$geneAnnotation$GeneSymbol))>1)   {

        message(paste(names(MetaData[i]), " gene probes as genes Symbols"))

        suma <- summary(rownames(MetaData[[i]])%in%MetaData$geneAnnotation$GeneSymbol)
        names(suma) <- c("Mode", "Gene not found", "Found")
        print(suma)
      }
    }




  message("-------------------------")
  message("Checking Samples of data clinic in MetaDataData from Raw gene expression matrix")
  message("-------------------------")
  for (i in c){

    if(all(ColN %in% rownames(MetaData[[i]]))==T) {   message(paste(MetaDataN[i]), " colnames : PASS") } else {  message(paste(MetaDataN[i]), " rownames : FAIL") }

  }











}
