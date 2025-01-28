#' ExportCSV Export Metadata inside object into ".csv" files
#'
#' @param Metadata a Metadata  data files
#' @return ".csv" files into working directory
#' @export
#' @import utils
#' @import R.utils
#' @import Matrix
#' @examples "non"
ExportCSV <- function (Metadata){



  if(is.null(Metadata)){stop("Need a Metadata List file")}
  if(!is.list(Metadata)){stop("Need a Metadata List file")}


  list.files.path = attributes(Metadata)$File.path
  if(is.null(list.files.path)){stop("Need a list file path for saving data")}
  if(!is.list(list.files.path)){stop(paste("list.files.path must be a list of file path whith Script, Raw genomic, Raw clinic, Processed and References directories in Parent Directory." ))}

  project = attributes(Metadata)$Project

  count <- 0
  count <- 0

  object <- summary(as.factor(attributes(Metadata)$Export))["Yes"]

  name <- names(Metadata)


  LF <- list.files(list.files.path$Project.VerifiedDataset)
  if(length(LF)!=0){
    df <- file.info(list.files(list.files.path$Project.VerifiedDataset, full.names = T))
    df$Filenames <- unlist(lapply(str_split(rownames(df),paste0(project,"/")),"[[",2))
    filepath <- rownames(df)
    filename <-  unlist(lapply(str_split(filepath,paste0(project,"/")),"[[",2))

    version <- na.omit(str_extract(filename,"V[0-9]"))

    Vnumber <- max(na.omit(as.numeric(str_extract(version,"([0-9]+).*$"))))+1


  } else {
      Vnumber = 1
      attributes(Metadata)$Version <- "V1"}







  NB.raw.clinic <- which(c(attributes(Metadata)$Data.Type=="Clinic" | attributes(Metadata)$Data.Type=="SamplesAnnot" ) & attributes(Metadata)$Export=="No" )
        if(length(NB.raw.clinic)>0) {
    count <- count+1
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


  NB.Samples.Patients.pheno <-   which(c(attributes(Metadata)$Data.Type=="Clinic" |attributes(Metadata)$Data.Type=="SamplesAnnot" ) & attributes(Metadata)$Export=="Yes" )

        if(length(NB.Samples.Patients.pheno)!=0) {

       for (j in NB.Samples.Patients.pheno) {

         z <-  Metadata[[j]]

         count <- count+1
         message("-------------------------------------------------")
         message(paste("Exporting", count, "/", object,"object: ",names(Metadata)[j]))

         if(Vnumber==1){
             filename <- paste0(list.files.path$Project.VerifiedDataset,"/",project,".",names(Metadata)[j],".V1",".csv")
             write.table(z,row.names = F ,file = filename, sep = ",")

             } else {

               filename <- paste0(list.files.path$Project.VerifiedDataset,"/",project,".",names(Metadata)[j],".V", Vnumber,".csv")
               write.table(z,row.names = F ,file = filename, sep = ",")
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
            filename <- paste0(list.files.path$Project.VerifiedDataset,"/",project,".",names(Metadata)[j], ".V1.csv")
            write.table(z,row.names = F ,file = filename, sep = ",")
            message(paste("Compressing"))
            R.utils::gzip(filename, destname=sprintf("%s.gz", filename), overwrite=T, remove=TRUE, BFR.SIZE=1e+07)
          }

          if(attributes(Metadata)$Omics.type=="Single.Cell"){

            filename <- paste0(list.files.path$Project.VerifiedDataset,"/",project,".",names(Metadata)[j],".V1.mtx")
            if(!class(Metadata[[j]])[1]=="dgTMatrix"){Metadata[[j]] = as.matrix(Metadata[[j]]) }
            writeMM(Matrix(as.matrix(Metadata[[j]]), sparse = T),file = filename)
            message(paste("Compressing"))
            R.utils::gzip(filename, destname=sprintf("%s.gz", filename), overwrite=T, remove=TRUE, BFR.SIZE=1e+07)
            gc()

            if(!"geneAnnot"%in%attributes(Metadata)$Data.Type){
            message("No geneAnnot file found. Exporting geneAnnot from count matrix.")
            filename.genes <- paste0(list.files.path$Project.VerifiedDataset,"/",project,".",names(Metadata)[j],".GenesAnnot.V1.csv")
            write.table(rownames(Metadata[[j]]),row.names = F ,col.names = F ,file = filename.genes, sep = ",")
            count = count+1
            message("-------------------------------------------------")
            message(paste("Exporting", count, "/", object,"object: ","geneAnnot", "file"))}


          if(!"CellsAnnot"%in%attributes(Metadata)$Data.Type){
            message("No CellsAnnot file found. Exporting CellsAnnot from count matrix.")
            filename.cells <- paste0(list.files.path$Project.VerifiedDataset,"/",project,".CellsAnnot", "V1.csv")
            write.table(data.frame("Cells"= colnames(Metadata[[j]])),row.names = F ,file = filename.cells, sep = ",")
            count = count+1
            message("-------------------------------------------------")
            message(paste("Exporting", count, "/", object,"object: ","CellsAnnot", "file"))} else {

              kk = which(attributes(Metadata)$Data.Type%in%"CellsAnnot"& attributes(Metadata)$Export=="Yes")
              if(length(kk)>0){
              filename.cells <- paste0(list.files.path$Project.VerifiedDataset,"/",project,".",names(Metadata)[kk], ".V1.csv")
              write.table(Metadata[[kk]],row.names = F ,file = filename.cells, sep = ",")
              count = count+1
              message("-------------------------------------------------")
              message(paste("Exporting", count, "/", object,"object: ","CellsAnnot", "file"))}

            }}} else { #Vnumber ==1

            if(attributes(Metadata)$Omics.type!="Single.Cell") {

              filename <- paste0(list.files.path$Project.VerifiedDataset,"/",project,".",names(Metadata)[j],".V",Vnumber, ".csv")
              write.table(z,row.names = F ,file = filename, sep = ",")
              message(paste("Compressing"))
              R.utils::gzip(filename, destname=sprintf("%s.gz", filename), overwrite=T, remove=TRUE, BFR.SIZE=1e+07)}


            if(attributes(Metadata)$Omics.type=="Single.Cell") {

              filename <- paste0(list.files.path$Project.VerifiedDataset,"/",project,".",names(Metadata)[j],".V",Vnumber,".mtx")

              if(!class(Metadata[[j]])[1]=="dgTMatrix"){Metadata[[j]] = as.matrix(Metadata[[j]]) }
              writeMM(Matrix(as.matrix(Metadata[[j]]), sparse = T),file = filename)
              message(paste("Compressing"))
              R.utils::gzip(filename, destname=sprintf("%s.gz", filename), overwrite=T, remove=TRUE, BFR.SIZE=1e+07)
              gc()

              if(!"geneAnnot"%in%attributes(Metadata)$Data.Type){
              message("No geneAnnot file found. Exporting geneAnnot from count matrix.")
              filename.genes <- paste0(list.files.path$Project.VerifiedDataset,"/",project,".",names(Metadata)[j],".GenesAnnot",".V",Vnumber,".csv")
              count = count+1
              message("-------------------------------------------------")
              message(paste("Exporting", count, "/", object,"object: ","geneAnnot", "file"))

              write.table(rownames(Metadata[[j]]),row.names = F ,col.names = F ,file = filename.genes, sep = ",")
              }

              if(!"CellsAnnot"%in%attributes(Metadata)$Data.Type){
                message("No CellsAnnot file found. Exporting CellsAnnot from count matrix.")
                filename.cells <- paste0(list.files.path$Project.VerifiedDataset,"/",project,".CellsAnnot",".V",Vnumber,".csv")
                count = count+1
                message("-------------------------------------------------")
                message(paste("Exporting", count, "/", object,"object: ","CellsAnnot", "file"))

                write.table(data.frame("Cells"= colnames(Metadata[[j]])),row.names = F ,file = filename.cells, sep = ",")} else {

                  filename.cells <- paste0(list.files.path$Project.VerifiedDataset,"/",project,".CellsAnnot",".V",Vnumber,".csv")
                  kk = which(attributes(Metadata)$Data.Type%in%"CellsAnnot"& attributes(Metadata)$Export=="Yes")
                  if(length(kk)>0){
                  count = count+1
                  message("-------------------------------------------------")
                  message(paste("Exporting", count, "/", object,"object: ","CellsAnnot", "file"))

                  write.table(Metadata[[kk]],row.names = F ,file = filename.cells, sep = ",")}

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

        filename <- paste0(list.files.path$Project.VerifiedDataset,"/",project,".",names(Metadata)[j],".V1.csv")
        write.table(z,row.names = F ,file = filename, sep = ",")
        message(paste("Compressing"))
        R.utils::gzip(filename, destname=sprintf("%s.gz", filename), overwrite=T, remove=TRUE, BFR.SIZE=1e+07)
      } else { #Vnumber==1



        filename <- paste0(list.files.path$Project.VerifiedDataset,"/",project,".",names(Metadata)[j],".V",Vnumber,".csv")
        write.table(z,row.names = F ,file = filename, sep = ",")
        message(paste("Compressing"))
        R.utils::gzip(filename, destname=sprintf("%s.gz", filename), overwrite=T, remove=TRUE, BFR.SIZE=1e+07)

        } #Vnumber !=1
            } #j in NGgeneannot
         }# NBgeneAnnot


  attributes(Metadata)$Version = paste0("V", Vnumber)

return(Metadata)


  }#function


