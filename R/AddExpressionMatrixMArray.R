#' AddExpressionMatrix MicroArray to Meta object
#'
#' @param Metadata Meta object
#' @param name if loca=True, names to apply in Metadata object slot
#' @param ExpressionMatrix name file of interest in path directory
#' @param force.replace set as F. T : replace an already object with the same name
#' @importFrom utils menu
#' @importFrom Matrix readMM
#' @import data.table
#' @return a data.frame in the Meta Object
#' @export
#'
#' @examples "none"
AddExpressionMatrixMArray<- function(Metadata=NULL,
                                     name = NULL,
                                     ExpressionMatrix = NULL,
                                     force.replace=F ) {

  Omics.type = "Microarray"

  path = attributes(Metadata)$File.path$Project.RawData

  if(is.null(Metadata)){stop("A Metadata object must be created with CreateDataset() function. See ?CreateDataset.")}

  if(is.null(ExpressionMatrix)){stop("ExpressionMatrix must be a character string or an environement object.")}


  if(is.null(name)){ stop("A Metadata object must be created with CreateDataset() function. See ?CreateDataset.") }

  if(is.null(attributes(Metadata)$Omics.type)){attributes(Metadata)$Omics.type=Omics.type}

  if(!attributes(Metadata)$Omics.type==Omics.type){
    warning(paste("Omics.type is", Omics.type,"different than attributes(Metadata)$Omics.type",
                  attributes(Metadata)$Omics.type, "\nWill be replace"))

    attributes(Metadata)$Omics.type = Omics.type

  }

  if(!is.list(Metadata)){stop("Metadata should be a list.")}



  ## if character
  if(inherits(ExpressionMatrix, "character")){

    message("Local import.")
    l <-length(names(Metadata))
    lf <- list.files(path)

    if(length(lf)>1){print(c(message("There is more than one files in Dir :"),lf))}

    if(all(str_detect(lf, ".rds|.txt|.csv|.tsv|.mtx", negate = FALSE)==F)){stop("No '*.rds' or '*.txt' or '*.csv' or '*.mtx' files in set directory. \n change path or add file")}

    if(!ExpressionMatrix%in%lf){stop(paste(ExpressionMatrix, "is not found in",path ))}

    filepath <- paste(path,ExpressionMatrix,sep="/")

    message(paste("Loading", ExpressionMatrix, "file"))
    if(str_detect(ExpressionMatrix, ".rds", negate = FALSE)){
      dt <- readRDS(filepath)}
    else {
      if(str_detect(ExpressionMatrix, ".txt", negate = FALSE)){

        dt <- suppressWarnings(as.data.frame(data.table::fread(filepath)))

        if(length(dt[,1])!=length(unique( dt[,1]))) { } else {  rownames(dt) <- dt[,1]}


      }  else {
        if(str_detect(ExpressionMatrix, ".csv", negate = FALSE)){
          dt <- suppressWarnings(as.data.frame(data.table::fread(filepath)))
          if(length(dt[,1])!=length(unique( dt[,1]))) { } else {  rownames(dt) <- dt[,1]}
        }else {

          if(str_detect(ExpressionMatrix, ".tsv", negate = FALSE)){
            dt <- suppressWarnings(as.data.frame(data.table::fread(filepath)))
            if(length(dt[,1])!=length(unique( dt[,1]))) { } else {  rownames(dt) <- dt[,1]}
          } else {


          }# tsv


        }#csv
      }#txt
    }#rds
  } else {

    if(inherits_any(ExpressionMatrix, c("data.frame", "matrix", "dgCMatrix" ,"dgTMatrix"))){
      dt = ExpressionMatrix}else { stop("Object set in ExpressionMatrix is not of class 'data.frame', 'matrix', 'dgCMatrix' ,'dgTMatrix'")}}

  if(length(Metadata)>=1) {
    l = length(Metadata)
    if(!all(str_detect(names(Metadata),name))==F){
      message("An Object with the same name already exist in MetaObject")
      if(force.replace==F){stop("set force.replace==T to subset object.")}
      message("Subsetting object.")
      Metadata[[name]] <- dt    } else { Metadata[[l+1]] <- dt
      names(Metadata)[l+1] <- name}

    if(length(attributes(Metadata)$Data.Type)<length(Metadata)){
      attributes(Metadata)$Data.Type <-  c(attributes(Metadata)$Data.Type, "Count")

     attributes(Metadata)$Export <- c(attributes(Metadata)$Export,"Yes")
    }

    return(Metadata)} # Metadat >1
  else {

    Metadata$mat <-  dt
    names(Metadata)[1] <- name

    if(length(attributes(Metadata)$Data.Type)<length(Metadata)){
      attributes(Metadata)$Data.Type <-  c("Count")
attributes(Metadata)$Export <- "Yes"
    return(Metadata)}
  } # Metadata = 0




} # function

