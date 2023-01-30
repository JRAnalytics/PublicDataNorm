#' AddClinic to Meta object
#'
#' @param Metadata Meta object
#' @param path dir path in which the GDC project is saved, or local files are saved
#' @param name.local.file if loca=True, names to apply in Metadata object slot
#' @param merge merge loaded data clinic with existing clinial data : cbind by rownames
#' @importFrom utils menu
#' @return a data.frame in the Meta Object
#' @export
#'
#' @examples "none"
AddClinic <- function(Metadata, path, merge = c(F,T), name.local.file) {
  if(!is.null(Metadata)){
    if(!is.list(Metadata)){
      stop("Metadata should be a list.")}

    setwd(path)

      l <-length(names(Metadata))
      lf <- list.files(path)
      lf <- lf[str_detect(lf, ".clinic")]

      if(all(str_detect(lf, ".rds", negate = FALSE)==F) & all(str_detect(lf, ".txt", negate = FALSE)==F) & all(str_detect(lf, "csv", negate = FALSE)==F) ){stop("No '*.rds' or '*.txt' '.csv' files in set directory. \n change path or add file")}


      if(str_detect(lf, ".rds", negate = FALSE)){Metadata[[l+1]] <- readRDS(lf) }

      if(str_detect(lf, ".txt|.csv|.tsv", negate = FALSE)){

        dt <- suppressWarnings(as.data.frame(data.table::fread(lf, na.strings = "")))
        rownames(dt) <- dt[,1]
        }

        if (merge == F){


            Metadata[[l+1]] <- dt

            names(Metadata)[l+1] <- name.local.file
        } else {

          clinic <- Metadata[which(str_detect(names(Metadata),"clinic"))][[1]]

          clinic <- cbind(clinic, dt[rownames(clinic),])
          Metadata[which(str_detect(names(Metadata),"clinic"))][[1]] <- clinic
          }

      return(Metadata)


   }else{ stop("No meta data object found")




  } #else meta data object not found



} # function
