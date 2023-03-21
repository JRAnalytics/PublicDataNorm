#' AddObjetToMeta
#'
#' @param Meta  a Meta Object to fill
#' @param object an object to add
#' @param name name to apply in Meta objet list
#' @param Data.type attribute c("Expression.Matrix","Patient.Clinical.data", "Samples.Clinical.data","geneAnnotation.file" )
#' @param Raw attribute c("Yes", "No")
#'
#' @return a Meta Object
#' @export
#'
#' @examples "none
AddObjetToMeta <- function(Meta, object, name ,Data.type = c("Expression.Matrix","Patient.Clinical.data", "Samples.Clinical.data","geneAnnotation.file" ), Raw = c("Yes", "No")){


  if(!is.list(Meta)) { stop("Meta object should be a list")}
  if(is.null(object)) { stop("Object is null")}
  if(is.null(name)) { stop("name should be specified, as character")}
  if(!is.character(name)) { stop("name should be specified as character")}
  if(is.null(Data.type)) { stop("Data.type should be specified")}
  if(is.null(Raw)) { stop("Raw should be specified")}
  if(!Data.type%in%c("Expression.Matrix","Patient.Clinical.data", "Samples.Clinical.data","geneAnnotation.file" )) { stop("Data.type should be specified from this values c('Expression.Matrix','Patient.Clinical.data', 'Samples.Clinical.data','geneAnnotation.file' ")}
  if(!Raw%in%c("Yes", "No")) { stop("Raw should be specified 'Yes' or 'No'")}


  l <- length(Meta)

  Meta[[l+1]] <- object
  names( Meta )[l+1] <- name

  attributes(Meta)$Data.Type[l+1] <- Data.type

  attributes(Meta)$Raw.data[l+1] <- Raw

return(Meta)

}
