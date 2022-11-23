#' ExportCSV Export MetaData inside object into ".csv" files
#'
#' @param MetaData a Meta  data files
#' @param list.files.path dirpath
#' @return ".csv" files into working directory
#' @export
#' @import utils
#' @import R.utils
#' @examples "non"
ExportCSV <- function (MetaData, list.files.path){

  if(is.null(MetaData)){stop("Need a MetaData List file")}
  if(!is.list(MetaData)){stop("Need a MetaData List file")}

  if(is.null(list.files.path)){stop("Need a list file path for saving data")}
  if(!is.list(list.files.path)){stop(paste("list.files.path must be a list of file path whith Script, Raw genomic, Raw clinic, Processed and References directories in Parent Directory." ))}


  for (i in list.files(list.files.path$RawDataDump)) {
        filename <- paste(c(list.files.path$RawDataDump,i),collapse = "/")
    gzip(filename, destname=sprintf("%s.gz", filename), overwrite=FALSE, remove=TRUE, BFR.SIZE=1e+07)

  }



  name <- names(MetaData)

  for (i in name) {

    if (!str_detect(toupper(i), "ANNOT")) {



         if(str_detect(toupper(i), "SAMPLE.PHENO")) {  filename <- paste0(list.files.path$PipelineDump,"/",i,".csv")
                                                        z <- cbind( MetaData[[i]])
                                                         write.csv(z,row.names = F ,file = filename)
                                                          gzip(filename, destname=sprintf("%s.gz", filename), overwrite=FALSE, remove=TRUE, BFR.SIZE=1e+07)}

          if(str_detect(toupper(i), "PATIENT.CLINIC")) {  filename <- paste0(list.files.path$PipelineDump,"/",i,".csv")
                                                            z <- cbind( MetaData[[i]])
                                                             write.csv(z,row.names = F ,file = filename)
                                                              gzip(filename, destname=sprintf("%s.gz", filename), overwrite=FALSE, remove=TRUE, BFR.SIZE=1e+07)}


           if(str_detect(toupper(i), "NORMALIZED.MATRIX")) { filename <- paste0(list.files.path$PipelineDump,"/",i,".csv")
                                                              z <- cbind("GeneSymbol" = rownames(MetaData[[i]]), MetaData[[i]])
                                                               write.csv(z,row.names = F ,file = filename)
                                                                gzip(filename, destname=sprintf("%s.gz", filename), overwrite=FALSE, remove=TRUE, BFR.SIZE=1e+07)}

    } else {

    z <- cbind("rownames" = rownames(MetaData[[i]]), MetaData[[i]][,colnames(MetaData[[i]])!="attributes"])
    if (str_detect(toupper(i), "GENEANNOTATION")){
      filename <- paste0(list.files.path$RawDataDump,"/",i,".csv")
      write.csv(z,row.names = F ,file = filename)
      gzip(filename, destname=sprintf("%s.gz", filename), overwrite=FALSE, remove=TRUE, BFR.SIZE=1e+07)}
    }



    }


}
