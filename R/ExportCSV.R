#' ExportCSV Export MetaData inside object into ".csv" files
#'
#' @param MetaData
#'
#' @return ".csv" files into working directory
#' @export
#' @import utils
#' @examples "non"
ExportCSV <- function (MetaData, list.files.path){

  if(is.null(MetaData)){stop("Need a MetaData List file")}
  if(!is.list(MetaData)){stop("Need a MetaData List file")}

  if(is.null(list.files.path)){stop("Need a list file path for saving data")}
  if(!is.list(list.files.path)){stop(paste("list.files.path must be a list of file path whith Script, Raw genomic, Raw clinic, Processed and References directories in Parent Directory." ))}


  name <- names(MetaData)

  for (i in name) {

    if (!str_detect(toupper(i), "ANNOT")) {



    if(str_detect(toupper(i), "RAW.MATRIX")) {   z <- cbind("GeneSymbol" = rownames(MetaData[[i]]), MetaData[[i]])
                                                write.csv(z,row.names = F ,file = paste0(list.files.path$RawGenomic,"/",i,".csv"))}

      if(str_detect(toupper(i), "RAW.CLINIC")) {  z <- cbind( MetaData[[i]])
                                                   write.csv(z,row.names = F ,file = paste0(list.files.path$RawClinic,"/",i,".csv"))}

        if(str_detect(toupper(i), "SAMPLE.PHENO")) {  z <- cbind( MetaData[[i]])
                                                       write.csv(z,row.names = F ,file = paste0(list.files.path$Processed,"/",i,".csv"))}

          if(str_detect(toupper(i), "PATIENT.CLINIC")) {  z <- cbind( MetaData[[i]])
                                                          write.csv(z,row.names = F ,file = paste0(list.files.path$Processed,"/",i,".csv"))}


           if(str_detect(toupper(i), "NORMALIZED.MATRIX")) { z <- cbind("GeneSymbol" = rownames(MetaData[[i]]), MetaData[[i]])
                                                            write.csv(z,row.names = F ,file = paste0(list.files.path$Processed,"/",i,".csv"))}

    } else {

    z <- cbind("rownames" = rownames(MetaData[[i]]), MetaData[[i]][,colnames(MetaData[[i]])!="attributes"])
    if (str_detect(toupper(i), "GENEANNOTATION")){write.csv(z,row.names = F ,file = paste0(list.files.path$References,"/",i,".csv"))}
    }



    }


}
