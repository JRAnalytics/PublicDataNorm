#' ExportTSV Export Metadata inside object into ".tsv" files
#'
#' @param Metadata a Metadata  data files
#' @param encoding see file encoding section "UTF-8" or "latin1"
#' @return ".tsv" files into working directory
#' @export
#' @import utils
#' @import R.utils
#' @import Matrix
#' @examples "non"
ExportTSV <- function (Metadata ,encoding = "UTF-8"){



  if(is.null(Metadata)){stop("Need a Metadata List file")}
  if(!is.list(Metadata)){stop("Need a Metadata List file")}


  list.files.path = attributes(Metadata)$File.path

  if(is.null(list.files.path)){stop("Need a list file path for saving data")}
  if(!is.list(list.files.path)){stop(paste("list.files.path must be a list of file path whith Script, Raw genomic, Raw clinic, Processed and References directories in Parent Directory." ))}


  project = attributes(Metadata)$Project

  count <- 0

  object <- summary(as.factor(attributes(Metadata)$Export))["Yes"]

  name <- names(Metadata)


  LF <- list.files(Verifiedpath(Metadata))
  if(length(LF)!=0){
    df <- file.info(list.files(Verifiedpath(Metadata), full.names = T))
    df$Filenames <- unlist(lapply(str_split(rownames(df),paste0(project,"/")),"[[",2))
    filepath <- rownames(df)
    filename <-  unlist(lapply(str_split(filepath,paste0(project,"/")),"[[",2))

    version <- na.omit(str_extract(filename,"V[0-9]"))

    Vnumber <- max(na.omit(as.numeric(str_extract(version,"([0-9]+).*$"))))+1

    if(Vnumber =="-Inf"){Vnumber=1}

  } else {
    Vnumber = 1
    attributes(Metadata)$Version <- "V1"}









  NB.raw.clinic <- which(c(attributes(Metadata)$Data.Type=="Clinic" | attributes(Metadata)$Data.Type=="SamplesAnnot" ) & attributes(Metadata)$Export=="No" )
  if(length(NB.raw.clinic)>0) {
    message("-------------------------------------------------")
    message(paste("Exporting", count, "/", object,"object: ",names(Metadata)[NB.raw.clinic],"data will not be exported"))
  }


  PatientLexic = lapply(ls(envir=.GlobalEnv), get)[lapply(lapply(ls(envir=.GlobalEnv), get), attr, "Lexic") == "Yes" & lapply(lapply(ls(envir=.GlobalEnv), get), attr, "Name")=="PatientsLexic" ][[1]]

  if(attributes(Metadata)$Omics.type=="Single.Cell"){

    CellsLexic = lapply(ls(envir=.GlobalEnv), get)[lapply(lapply(ls(envir=.GlobalEnv), get), attr, "Lexic") == "Yes" & lapply(lapply(ls(envir=.GlobalEnv), get), attr, "Name")=="CellsLexic" ][[1]]

  }

    SamplesLexic = lapply(ls(envir=.GlobalEnv), get)[lapply(lapply(ls(envir=.GlobalEnv), get), attr, "Lexic") == "Yes" & lapply(lapply(ls(envir=.GlobalEnv), get), attr, "Name")=="SamplesLexic" ][[1]]


    if(exists("PatientLexic", mode= "any" )) {object = object+1}
    if(exists("SamplesLexic", mode= "any" )) {object = object+1}
    if(exists("CellsLexic", mode= "any" )) {object = object+1}

    message(paste0("Exporting Version V", Vnumber))

    message(paste("Exporting", object, "objects"))

  if(exists("PatientLexic", mode= "any" )) {
    count <- count+1
    message("-------------------------------------------------")
    message(paste("Exporting", count, "/", object,"object: ","PatientLexic"))

    if (file.exists(paste0(list.files.path$Project.Processes,"/",project,".PatientLexic.txt"))) {
      #Delete file if it exists
      file.remove(paste0(list.files.path$Project.Processes,"/",project,".PatientLexic.txt"))
    }

    PatientLexic <- lapply(PatientLexic, function(x) {c(x[1],x)}) #Mandatory to duplicated listName in the listed values.
    lapply(PatientLexic, write, paste0(list.files.path$Project.Processes,"/",project,".PatientLexic.txt"), append=TRUE, ncolumns=1000 ) #write a ".txt" file without listNames
  }


  if(exists("SamplesLexic", mode= "any" )) {
    count <- count+1
    message("-------------------------------------------------")
    message(paste("Exporting", count, "/", object,"object: ","SamplesLexic"))

    if (file.exists(paste0(list.files.path$Project.Processes,"/",project,".SamplesLexic.txt"))) {
      #Delete file if it exists
      file.remove(paste0(list.files.path$Project.Processes,"/",project,".SamplesLexic.txt"))
    }

    SamplesLexic <- lapply(SamplesLexic, function(x) {c(x[1],x)}) #Mandatory to duplicated listName in the listed values.
    lapply(SamplesLexic, write, paste0(list.files.path$Project.Processes,"/",project,".SamplesLexic.txt"), append=TRUE, ncolumns=1000 ) #write a ".txt" file without listNames
  }

if(exists("CellsLexic", mode= "any" )) {
  count <- count+1
  message("-------------------------------------------------")
  message(paste("Exporting", count, "/", object,"object: ","CellsLexic"))

  if (file.exists(paste0(list.files.path$Project.Processes,"/",project,".CellsLexic.txt"))) {
    #Delete file if it exists
    file.remove(paste0(list.files.path$Project.Processes,"/",project,".CellsLexic.txt"))
  }

  CellsLexic <- lapply(CellsLexic, function(x) {c(x[1],x)}) #Mandatory to duplicated listName in the listed values.
  lapply(CellsLexic, write, paste0(list.files.path$Project.Processes,"/",project,".CellsLexic.txt"), append=TRUE, ncolumns=1000 ) #write a ".txt" file without listNames
}



  NB.Samples.Patients.pheno <-   which(c(attributes(Metadata)$Data.Type=="Clinic" |attributes(Metadata)$Data.Type=="SamplesAnnot" ) & attributes(Metadata)$Export=="Yes" )

  if(length(NB.Samples.Patients.pheno)!=0) {

    for (j in NB.Samples.Patients.pheno) {

      z <-  Metadata[[j]]

      count <- count+1
      message("-------------------------------------------------")
      message(paste("Exporting", count, "/", object,"object: ",names(Metadata)[j]))

      if(Vnumber==1){
        filename <- paste0(Verifiedpath(Metadata),"/",project,".",names(Metadata)[j],".V1",".tsv")
        write.table(z,row.names = F ,file = filename, sep = "\t", fileEncoding = encoding)

      } else {

        filename <- paste0(Verifiedpath(Metadata),"/",project,".",names(Metadata)[j],".V", Vnumber,".tsv")
        write.table(z,row.names = F ,file = filename, sep = "\t", fileEncoding = encoding)
      }



    } #for J in NBsamples

  } # if( length(NB.Samples.Patients.pheno)!=0)


  NB.Count <-  which(attributes(Metadata)$Data.Type=="Count")

  if(length(NB.Count)!=0) {

    for (j in NB.Count) {

      count <- count+1
      message("-------------------------------------------------")
      message(paste("Exporting", count, "/", object,"object: ",names(Metadata)[j], "file"))


      if(attributes(Metadata)$Omics.type!="Single.Cell"){ z <- try(cbind("GeneSymbol" = rownames(Metadata[[j]]), Metadata[[j]]),silent = T)}
      if(attributes(Metadata)$Omics.type=="Single.Cell"){ z = Metadata[[j]]}

      if(Vnumber==1){

        if(attributes(Metadata)$Omics.type!="Single.Cell"){
          filename <- paste0(Verifiedpath(Metadata),"/",project,".",names(Metadata)[j], ".V1.tsv")
          write.table(z,row.names = F ,file = filename, sep = "\t")
          message(paste("Compressing"))
          R.utils::gzip(filename, destname=sprintf("%s.gz", filename), overwrite=T, remove=TRUE, BFR.SIZE=1e+07)
        }

        if(attributes(Metadata)$Omics.type=="Single.Cell"){

          filename <- paste0(Verifiedpath(Metadata),"/","V1.matrix.mtx")
          if(class(Metadata[[j]])[1]== "dgCMatrix") {
            writeMM(Matrix((Metadata[[j]]), sparse = T),file = filename)} else{

            if(!class(Metadata[[j]])[1]=="dgTMatrix"){Metadata[[j]] = as.matrix(Metadata[[j]])}

            writeMM(Matrix(Metadata[[j]], sparse = T),file = filename)
          }




          message(paste("Compressing"))
          R.utils::gzip(filename, destname=sprintf("%s.gz", filename), overwrite=T, remove=TRUE, BFR.SIZE=1e+07)
          gc()

          if(!"geneAnnot"%in%attributes(Metadata)$Data.Type){
            message("No geneAnnot file found. Exporting geneAnnot from count matrix.")
            filename.genes <- paste0(Verifiedpath(Metadata),"/","V1.features.tsv")
            write.table(rownames(Metadata[[j]]),row.names = F ,col.names = F ,file = filename.genes, sep = "\t")
            R.utils::gzip(filename.genes, destname=sprintf("%s.gz", filename.genes), overwrite=T, remove=TRUE, BFR.SIZE=1e+07)

            count = count+1
            message("-------------------------------------------------")
            message(paste("Exporting", count, "/", object,"object: ","geneAnnot", "file"))}


          if(!"CellsAnnot"%in%attributes(Metadata)$Data.Type){
            message("No CellsAnnot file found. Exporting CellsAnnot from count matrix.")
            filename.cells <- paste0(Verifiedpath(Metadata),"/", "V1.barcodes.tsv")
            write.table(data.frame("Cells"= colnames(Metadata[[j]])),row.names = F ,file = filename.cells, sep = "\t")
            R.utils::gzip(filename.cells, destname=sprintf("%s.gz", filename.cells), overwrite=T, remove=TRUE, BFR.SIZE=1e+07)
            count = count+1
            message("-------------------------------------------------")
            message(paste("Exporting", count, "/", object,"object: ","CellsAnnot", "file"))} else {

              kk = which(attributes(Metadata)$Data.Type%in%"CellsAnnot"& attributes(Metadata)$Export=="Yes")
              if(length(kk)>0){
                for (i in kk){
                filename.cells <- paste0(Verifiedpath(Metadata),"/", "V1.barcodes.tsv")
                write.table(Metadata[[i]],row.names = F ,file = filename.cells, sep = "\t")
                R.utils::gzip(filename.cells, destname=sprintf("%s.gz", filename.cells), overwrite=T, remove=TRUE, BFR.SIZE=1e+07)

                count = count+1
                message("-------------------------------------------------")
                message(paste("Exporting", count, "/", object,"object: ","CellsAnnot", "file"))}}

            }}} else { #Vnumber ==1

              if(attributes(Metadata)$Omics.type!="Single.Cell") {

                filename <- paste0(Verifiedpath(Metadata),"/",project,".",names(Metadata)[j],".V",Vnumber, ".tsv")
                write.table(z,row.names = F ,file = filename, sep = "\t")
                message(paste("Compressing"))
                R.utils::gzip(filename, destname=sprintf("%s.gz", filename), overwrite=T, remove=TRUE, BFR.SIZE=1e+07)}


              if(attributes(Metadata)$Omics.type=="Single.Cell") {

                filename <- paste0(Verifiedpath(Metadata),"/","V",Vnumber,".matrix.mtx")


                if(class(Metadata[[j]])[1]== "dgCMatrix") {
                  writeMM(Matrix((Metadata[[j]]), sparse = T),file = filename)} else{

                    if(!class(Metadata[[j]])[1]=="dgTMatrix"){Metadata[[j]] = as.matrix(Metadata[[j]])}

                    writeMM(Matrix(Metadata[[j]], sparse = T),file = filename)
                  }


                writeMM(Metadata[[j]],file = filename)
                message(paste("Compressing"))
                R.utils::gzip(filename, destname=sprintf("%s.gz", filename), overwrite=T, remove=TRUE, BFR.SIZE=1e+07)
                gc()

                if(!"geneAnnot"%in%attributes(Metadata)$Data.Type){
                  message("No geneAnnot file found. Exporting geneAnnot from count matrix.")
                  filename.genes <- paste0(Verifiedpath(Metadata),"/","V",Vnumber,".features.tsv")
                  count = count+1
                  message("-------------------------------------------------")
                  message(paste("Exporting", count, "/", object,"object: ","geneAnnot", "file"))

                  write.table(rownames(Metadata[[j]]),row.names = F ,col.names = F ,file = filename.genes, sep = "\t")
                  R.utils::gzip(filename.genes, destname=sprintf("%s.gz", filename.genes), overwrite=T, remove=TRUE, BFR.SIZE=1e+07)
                }

                if(!"CellsAnnot"%in%attributes(Metadata)$Data.Type){
                  message("No CellsAnnot file found. Exporting CellsAnnot from count matrix.")
                  filename.cells <- paste0(Verifiedpath(Metadata),"/",project,".CellsAnnot","V",Vnumber,".barcodes.tsv")
                  count = count+1
                  message("-------------------------------------------------")
                  message(paste("Exporting", count, "/", object,"object: ","CellsAnnot", "file"))

                  write.table(data.frame("Cells"= colnames(Metadata[[j]])),row.names = F ,file = filename.cells, sep = "\t")
                  R.utils::gzip(filename.cells, destname=sprintf("%s.gz", filename.cells), overwrite=T, remove=TRUE, BFR.SIZE=1e+07)

                  } else {

                    filename.cells <- paste0(Verifiedpath(Metadata),"/","V",Vnumber,".barcodes.tsv")
                    kk = which(attributes(Metadata)$Data.Type%in%"CellsAnnot"& attributes(Metadata)$Export=="Yes")
                    if(length(kk)>0){
                      for (i in kk){
                      count = count+1
                      message("-------------------------------------------------")
                      message(paste("Exporting", count, "/", object,"object: ",names(Metadata)[[i]], "file"))

                      write.table(Metadata[[i]],row.names = F ,file = filename.cells, sep = "\t")
                      R.utils::gzip(filename.cells, destname=sprintf("%s.gz", filename.cells), overwrite=T, remove=TRUE, BFR.SIZE=1e+07)}}

                  }





              }



            } #Vnumber !=1

    }#j in NB.count
  } # if(length(NB.Count)!=0)


  NB.geneAnnot<-  which(attributes(Metadata)$Data.Type=="geneAnnot")
  if(length(NB.geneAnnot)!=0) {

    for (j in NB.geneAnnot) {

      count <- count+1
      message("-------------------------------------------------")
      message(paste("Exporting", count, "/", object,"object: ",names(Metadata)[j], "file"))

      z <- Metadata[[j]]

      if(Vnumber==1){
        if(attributes(Metadata)$Omics.type=="Single.Cell") {
        filename <- paste0(Verifiedpath(Metadata),"/","V1.features.tsv")
    }else{
          filename <- paste0(list.files.path$Project.VerifiedDataset,"/",project,".",names(Metadata)[j],".V1.tsv")
        }
        write.table(z,row.names = F ,file = filename, sep = ",")
        message(paste("Compressing"))
        R.utils::gzip(filename, destname=sprintf("%s.gz", filename), overwrite=T, remove=TRUE, BFR.SIZE=1e+07)
      } else { #Vnumber==1


        if(attributes(Metadata)$Omics.type=="Single.Cell") {
          filename <- paste0(list.files.path$Project.VerifiedDataset,"/","V",Vnumber,".features.tsv")
        }else{
        filename <- paste0(list.files.path$Project.VerifiedDataset,"/",project,".",names(Metadata)[j],".V",Vnumber,".tsv")}
        write.table(z,row.names = F ,file = filename, sep = ",")
        message(paste("Compressing"))
        R.utils::gzip(filename, destname=sprintf("%s.gz", filename), overwrite=T, remove=TRUE, BFR.SIZE=1e+07)

      } #Vnumber !=1
    } #j in NGgeneannot
  }# NBgeneAnnot


  attributes(Metadata)$Version = paste0("V", Vnumber)

  return(Metadata)


}#function


