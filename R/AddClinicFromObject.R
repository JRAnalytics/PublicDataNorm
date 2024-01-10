#' AddClinicFromObject
#'
#' @param Metadata  a Metadata Object to fill
#' @param object an object to add
#' @param name name to apply in Metadata objet list
#' @param SamplesFilter Default F, select only samples found in both the clinical data and the Expression Matrix.
#' @param Data.type attribute  c("Clinic", "SamplesAnnot","CellAnnot" )
#' @param Export data to export after cleaning c("Yes", "No")
#' @import matrixStats
#'
#' @return a Metadata Object
#' @export
#'
#' @examples "none"
#'
AddClinicFromObject  <- function(Metadata,
                           object,
                           name,
                           SamplesFilter = F,
                           Data.type = c("Clinic", "SamplesAnnot","CellAnnot"),
                           Export = c("Yes", "No")){



  if(!is.list(Metadata)) { stop("Metadata object should be a list")}
  if(is.null(object)) { stop("Object is null")}
  if(is.null(name)) { stop("name should be specified, as character")}
  if(!is.character(name)) { stop("name should be specified as character")}
  if(is.null(Data.type)) { stop("Data.type should be specified")}
  if(is.null(Export)) { stop("Export should be specified")}
  if(!Data.type%in%c("Count","Clinic", "SamplesAnnot","geneAnnot","CellAnnot" )) { stop("Data.type should be specified from this values c('Count','Clinic', 'SamplesAnnot','geneAnnot','CellAnnot' ) ")}
  if(!Export%in%c("Yes", "No")) { stop("Export should be specified 'Yes' or 'No'")}



    l <- length(Metadata)


    zz <- which(attributes(Metadata)$Data.Type=="Count")[1]
    samples <- colnames(Metadata[[zz]])
    val <- samples[which(samples%in%as.matrix(object))[1]]
    colT <- which(matrixStats::colAnys(as.matrix(object),value = val))

    if(length(colT)>1){ colT <- colT[1]}

    message("Found Samples :")
    print(summary(samples%in%as.matrix(object)))
    message("Unfound Samples :")
    print(samples[!samples%in%as.matrix(object)])

    if(SamplesFilter==T){
      message("Selecting only Samples present in both Count and clinical.data.")

      # Diviser chaque élément de la colonne en un vecteur de sous-chaînes
      substrings <- strsplit(as.data.frame(object)[,colT], ";")
      # Vérifier si chaque valeur du vecteur est présente dans chaque vecteur de sous-chaînes
      est_present <- sapply(substrings, function(x) any(colnames(Metadata[[zz]]) %in% x))

      # Récupérer les index de position correspondants
      indices <- which(est_present)
      object <- object[indices,]


    }





  Metadata[[l+1]] <- object
  names( Metadata )[l+1] <- name

  attributes(Metadata)$Data.Type[l+1] <- Data.type

  attributes(Metadata)$Export[l+1] <- Export



  return(Metadata)

}
