#' AddExpressionMatrixSC from .mtx to Meta object
#'
#' @param Metadata Metaobject
#' @param name if loca=True, names to apply in Metadata object slot
#' @param ExpressionMatrix an object or a character string name file of interest in path directory from .mtx exclusive
#' @param Cell.file an object or a character string name of single cell cell annotation file specific to expression matrix. Not mandatory.
#' @param Genes.file an object or a character string name of single cell gene annotation file specific to expression matrix. Not mandatory.
#' @param force.replace set as F. T : replace an already object with the same name
#' @param setID.cellAnnotColumn  a character string or numeric : column in CellAnnot to fetch colnames from Count matrix.
#' @importFrom utils menu
#' @importFrom Matrix readMM
#' @import data.table
#' @import rlang
#' @return a data.frame in the Meta Object
#' @export
#'
#' @examples "none"
AddExpressionMatrixSC <- function(Metadata=NULL,
                                Cell.file=NULL,
                                Genes.file=NULL,
                                setID.cellAnnotColumn = NULL,
                                name,
                                ExpressionMatrix = NULL,
                                force.replace=F ) {
  Omics.type = "Single.Cell"
  path = Rawpath(Metadata)

  if(is.null(Metadata)){stop("A Metadata object must be created with CreateDataset() function. See ?CreateDataset.")}

  if(is.null(ExpressionMatrix)){stop("ExpressionMatrix must be a character string or an environement object.")}

  if(!is.list(Metadata)){stop("Metadata should be a list.")}

  if(is.null(attributes(Metadata)$Omics.type)){attributes(Metadata)$Omics.type=Omics.type}

  if(!attributes(Metadata)$Omics.type==Omics.type){warning(paste("Omics.type is", Omics.type,"different than attributes(Metadata)$Omics.type",
      attributes(Metadata)$Omics.type, "\nWill be replace"))

      attributes(Metadata)$Omics.type = Omics.type

    }



    l <-length(names(Metadata))

    if(inherits(ExpressionMatrix, "character")){
    lf <- list.files(path)

    if(length(lf)>1){print(c(message("There is more than one files in Dir :"),lf))}

    if(all(str_detect(lf, ".rds|.txt|.csv|.tsv|.mtx", negate = FALSE)==F)){stop("No '*.rds' or '*.txt' or '*.csv' or '*.mtx' files in set directory. \n change path or add file")}



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

              if(str_detect(ExpressionMatrix, ".mtx", negate = FALSE)){

                dt <- readMM(filepath)
                gc()
                #dt = as.matrix(dt)
                #gc()

            } # mtx
          }# tsv


        }#csv
      }#txt
    }#rds
      }#ExpressionMatrix caract
    else {

      if(rlang::inherits_any(ExpressionMatrix, c("data.frame", "matrix", "dgCMatrix" ,"dgTMatrix"))){
        dt = ExpressionMatrix

             }else { stop("Object set in ExpressionMatrix is not of class 'data.frame', 'matrix', 'dgCMatrix' ,'dgTMatrix'")}}



    if(is.null(colnames(dt))){
      message(paste(ExpressionMatrix,"has no colnames. A Cell.csv file may be associated in raw data directory."))} else {  colnames(dt) = gsub("_","-", colnames(dt))}


    if(!"CellsAnnot" %in% attributes(Metadata)){
      if(!is.null(Cell.file)){

        if(rlang::inherits_any(Cell.file, c("data.frame", "matrix"))){Cells = Cell.file} else {

          if(inherits(Cell.file, "character")){
            message(paste("Loading",Cell.file ))
            Cells <- as.data.frame(data.table::fread(file.path(path,Cell.file)))
            if("cell_name"%in%colnames(Cells)){
              rownames(Cells)  = Cells$cell_name} else {
                message("Cell.file has no colnames specified 'cell_name', the first collumn will be used.\n Please check file before adding cell file.")

                Cells <- as.data.frame(data.table::fread(file.path(path,Cell.file),header = F))
                rownames(Cells) = Cells[,1]}}else {stop("Cell.file is not a character string or an environment object as data.frame or matrix.")}}



        if(length(rownames(Cells))==length(colnames(dt))) {  colnames(dt) = rownames(Cells)} else
        { stop(paste(Cell.file, "has not the same number of cells than column of expression matrix."))
        }


        if(is.null(setID.cellAnnotColumn)){stop("setID.cellAnnotColumn mus be specify")}
        if(inherits(setID.cellAnnotColumn,"character")){
        if(!setID.cellAnnotColumn %in%colnames(Cells) ){stop(paste(setID.cellAnnotColumn, "is not found in colnames of Cell.File"))}
        Cells$CellsBarcode = Cells[,setID.cellAnnotColumn]
        Cells$CellsBarcode = gsub("[[:punct:]]","-", Cells$CellsBarcode)
        colnames(dt) = Cells$CellsBarcode}
        if(inherits(setID.cellAnnotColumn,"numeric")){Cells$CellsBarcode = Cells[,setID.cellAnnotColumn]
        Cells$CellsBarcode = gsub("[[:punct:]]","-", Cells$CellsBarcode)
        colnames(dt) = Cells$CellsBarcode}



      }
      }



    if(is.null(rownames(dt))){
      message(paste(ExpressionMatrix,"has no rownames A Genes.csv file may be associated in raw data directory."))}

        if(!is.null(Genes.file)){
          if(rlang::inherits_any(Genes.file, c("data.frame", "matrix"))){Genes =Genes.file }else {

            if(inherits(Genes.file, "character")){
              message(paste("Loading",Genes.file ))
              Genes <- as.data.frame(data.table::fread(file.path(path,Genes.file), header = F))}
            else {stop("Genes.file is not a character string or an environment object as data.frame or matrix.")}}

          if(length(Genes[,1])==nrow(dt)){

            rownames(dt) =  Genes[,1]
          } else { stop(paste(Genes.file, "has not the same length as rows of expression matrix"))}
        }





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
     attributes(Metadata)$Cleaned = c(attributes(Metadata)$Cleaned, "No")}


      if(!is.null(Cell.file)){
        Metadata$CellsAnnot = Cells
        attributes(Metadata)$Data.Type = c(attributes(Metadata)$Data.Type, "CellsAnnot")
        attributes(Metadata)$Export = c(attributes(Metadata)$Export, "No")
        attributes(Metadata)$Cleaned = c(attributes(Metadata)$Cleaned,"No")

      }

      if(!is.null(Genes.file)){
        Metadata$geneAnnotation = Genes
        attributes(Metadata)$Data.Type = c(attributes(Metadata)$Data.Type, "geneAnnot")
        attributes(Metadata)$Export = c(attributes(Metadata)$Export, "Yes")
        attributes(Metadata)$Cleaned = c(attributes(Metadata)$Cleaned,"No")

      }

      return(Metadata)} # Metadat >1


    else {

      Metadata$mat <-  dt
      names(Metadata)[1] <- name

      if(length(attributes(Metadata)$Data.Type)<length(Metadata)){
        attributes(Metadata)$Data.Type <-  c("Count")
       attributes(Metadata)$Export <- "Yes"
       attributes(Metadata)$Cleaned = c("No")}

      if(!is.null(Cell.file)){
        Metadata$CellsAnnot = Cells
        attributes(Metadata)$Data.Type = c(attributes(Metadata)$Data.Type, "CellsAnnot")
        attributes(Metadata)$Export = c(attributes(Metadata)$Export, "No")
        attributes(Metadata)$Cleaned = c(attributes(Metadata)$Cleaned,"No")

      }

      if(!is.null(Genes.file)){
        Metadata$geneAnnotation = Genes
        attributes(Metadata)$Data.Type = c(attributes(Metadata)$Data.Type, "geneAnnot")
        attributes(Metadata)$Export = c(attributes(Metadata)$Export, "Yes")
        attributes(Metadata)$Cleaned = c(attributes(Metadata)$Cleaned,"No")

      }

      return(Metadata)
    } # Metadata = 0


} # function
