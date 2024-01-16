#' AddExpressionMatrix to Meta object
#'
#' @param Metadata Metaobject
#' @param Raw.file.path dir path in which the raw single cell exrpession is saved
#' @param name if loca=True, names to apply in Metadata object slot
#' @param name.local.file name file of interest in path directory
#' @param Cell.file name of single cell cell annotation file specific to expression matrix. Not mandatory.
#' @param Genes.file name of single cell gene annotation file specific to expression matrix. Not mandatory.
#' @param force.replace set as F. T : replace an already object with the same name
#' @param Export TRUE or FALSE. If data to be Exported, set T.
#' @importFrom utils menu
#' @importFrom Matrix readMM
#' @import data.table
#' @return a data.frame in the Meta Object
#' @export
#'
#' @examples "none"
AddExpressionMatrixSC <- function(Metadata=NULL,
                                Cell.file=NULL,
                                Genes.file=NULL,
                                Raw.file.path,
                                name,
                                Export=T,
                                name.local.file = NULL,
                                force.replace=F ) {
  Omics.type = "Single.Cell"
  path =Raw.file.path

  if(is.null(Metadata)){
    Metadata = list()
    attributes(Metadata)$Omics.type = Omics.type
  }

  if(!is.null(Metadata)){

    if(!is.list(Metadata)){
      stop("Metadata should be a list.")}

    if(!attributes(Metadata)$Omics.type==Omics.type){
      warning(paste("Omics.type is", Omics.type,"different than attributes(Metadata)$Omics.type",
                    attributes(Metadata)$Omics.type, "\nWill be replace"))

      attributes(Metadata)$Omics.type = Omics.type

    }

    setwd(path)

    l <-length(names(Metadata))
    lf <- list.files(path)

    if(length(lf)>1){print(c(message("There is more than one files in Dir :"),lf))}

    if(all(str_detect(lf, ".rds|.txt|.csv|.tsv|.mtx", negate = FALSE)==F)){stop("No '*.rds' or '*.txt' or '*.csv' or '*.mtx' files in set directory. \n change path or add file")}

    if(!is.null(name.local.file)) {

      filepath <- paste(path,name.local.file,sep="/")
      message(paste("Loading", name.local.file, "file"))
      if(str_detect(name.local.file, ".rds", negate = FALSE)){
        dt <- readRDS(filepath)}
      else {
        if(str_detect(name.local.file, ".txt", negate = FALSE)){

          dt <- suppressWarnings(as.data.frame(data.table::fread(filepath)))

          if(length(dt[,1])!=length(unique( dt[,1]))) { } else {  rownames(dt) <- dt[,1]}


        }  else {
          if(str_detect(name.local.file, ".csv", negate = FALSE)){
            dt <- suppressWarnings(as.data.frame(data.table::fread(filepath)))
            if(length(dt[,1])!=length(unique( dt[,1]))) { } else {  rownames(dt) <- dt[,1]}
          }else {

            if(str_detect(name.local.file, ".tsv", negate = FALSE)){
              dt <- suppressWarnings(as.data.frame(data.table::fread(filepath)))
              if(length(dt[,1])!=length(unique( dt[,1]))) { } else {  rownames(dt) <- dt[,1]}
            } else {

              if(str_detect(name.local.file, ".mtx", negate = FALSE)){

                dt <- readMM(filepath)
                gc()
                #dt = as.matrix(dt)
                #gc()

                if(is.null(colnames(dt))){
                  message(paste(name.local.file,"has no colnames. A Cell.csv file may be associated in raw data directory."))

                  if(!is.null(Cell.file)){
                    message(paste("Loading",Cell.file ))
                    Cells <- fread(file.path(path,Cell.file))
                    if("cell_name"%in%colnames(Cells)){
                      rownames(Cells)  = Cells$cell_name} else {
                        message("Cell.file has no colnames specified 'cell_name', the first collumn will be used.\n Please check file before adding cell file.")

                        Cells <- as.data.frame(fread(file.path(path,Cell.file),header = F))
                        rownames(Cells) = Cells[,1]}


                    if(length(rownames(Cells)==length(colnames(dt)))) {  colnames(dt) = rownames(Cells)} else
                    { stop(paste(Cell.file, "has not the same number of cells than column of expression matrix."))
                    }



                  }
                }


                if(is.null(rownames(dt))){
                  message(paste(name.local.file,"has no rownames A Genes.csv file may be associated in raw data directory."))}

                if(!is.null(Genes.file)){
                  message(paste("Loading",Genes.file ))
                  Genes <- as.data.frame(fread(file.path(path,Genes.file), header = F))

                  if(length(Genes[,1])==nrow(dt)){

                    rownames(dt) =  Genes[,1]
                  } else { stop(paste(Genes.file, "has not the same length as rows of expression matrix"))}
                }
              }

            } # mtx
          }# tsv


        }#csv
      }#txt
    }#rds


    if(length(Metadata)>=1) {

      if(!all(str_detect(names(Metadata),paste0(name,".matrix")))==F){
        message("An Object with the same name already exist in MetaObject")
        if(force.replace==F){stop("set force.replace==T to subset object.")}
        message("Subsetting object.")
        Metadata[[paste0(name,".matrix")]] <- dt    } else { Metadata[[l+1]] <- dt
        names(Metadata)[l+1] <- paste0(name,".matrix")}

      if(length(attributes(Metadata)$Data.Type)<length(Metadata)){
        attributes(Metadata)$Data.Type <-  c(attributes(Metadata)$Data.Type, "Count")

        if(Export==T){attributes(Metadata)$Export <- c(attributes(Metadata)$Export,"Yes") } else {attributes(Metadata)$Export <- c(attributes(Metadata)$Export,"No") }
      }


      if(!is.null(Cell.file)){

        Metadata$CellsAnnot = Cells
        attributes(Metadata)$Data.Type = c(attributes(Metadata)$Data.Type, "CellsAnnot")
        attributes(Metadata)$Export = c(attributes(Metadata)$Export, "Yes")

      }

      if(!is.null(Genes.file)){
        Metadata$geneAnnot = Genes
        attributes(Metadata)$Data.Type = c(attributes(Metadata)$Data.Type, "geneAnnot")
        attributes(Metadata)$Export = c(attributes(Metadata)$Export, "Yes")

      }

      return(Metadata)} # Metadat >1


    else {

      Metadata$mat <-  dt
      names(Metadata)[1] <- paste0(name,".matrix")

      if(length(attributes(Metadata)$Data.Type)<length(Metadata)){
        attributes(Metadata)$Data.Type <-  c("Count")
        if(Export==T){attributes(Metadata)$Export <- "Yes" } else {attributes(Metadata)$Export <- "No" }}

      if(!is.null(Cell.file)){
        Metadata$CellsAnnot = Cells
        attributes(Metadata)$Data.Type = c(attributes(Metadata)$Data.Type, "CellsAnnot")
        attributes(Metadata)$Export = c(attributes(Metadata)$Export, "Yes")

      }

      if(!is.null(Genes.file)){
        Metadata$geneAnnot = Genes
        attributes(Metadata)$Data.Type = c(attributes(Metadata)$Data.Type, "geneAnnot")
        attributes(Metadata)$Export = c(attributes(Metadata)$Export, "Yes")

      }

      return(Metadata)
    } # Metadata = 0




  } else { stop("No local file name entered.")}

















} # function
