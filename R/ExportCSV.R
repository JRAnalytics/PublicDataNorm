#' ExportCSV Export MetaData inside object into ".csv" files
#'
#' @param MetaData a MetaData  data files
#' @param list.files.path dirpath
#' @param project project
#' @return ".csv" files into working directory
#' @export
#' @import utils
#' @import R.utils
#' @import Matrix
#' @examples "non"
ExportCSV <- function (MetaData, list.files.path, project){

  if(is.null(MetaData)){stop("Need a MetaData List file")}
  if(!is.list(MetaData)){stop("Need a MetaData List file")}

  if(is.null(list.files.path)){stop("Need a list file path for saving data")}
  if(!is.list(list.files.path)){stop(paste("list.files.path must be a list of file path whith Script, Raw genomic, Raw clinic, Processed and References directories in Parent Directory." ))}

  count <- 0
  object <- length(MetaData)+2
  if(attributes(MetaData)$Omics.type=="Single.Cell") {object = object+1}
  name <- names(MetaData)

  Vnumber <- NA

  LF <- list.files(list.files.path$Propject.VerifiedDataset)
  if(length(LF)!=0){
    df <- file.info(list.files(list.files.path$Propject.VerifiedDataset, full.names = T))
    df$Filenames <- unlist(lapply(str_split(rownames(df),paste0(project,"/")),"[[",2))
    filepath <- rownames(df)
    filename <-  unlist(lapply(str_split(filepath,paste0(project,"/")),"[[",2))

    if(attributes(MetaData)$Omics.type!="Single.Cell") {
      filename = filename[str_detect(filename,".csv")]
      filename2 <- unlist(lapply(str_split(filename,".csv"),"[[",1))}
    if(attributes(MetaData)$Omics.type=="Single.Cell") {
      filename = filename[str_detect(filename,".mtx")]
      filename2 <- unlist(lapply(str_split(filename,".mtx"),"[[",1))

    }

    version <- str_extract(filename2,"V[0-9]*")
    Vnumber <- max(as.numeric(str_extract(version,"([0-9]+).*$")))+1
    Vnumber2 <- Vnumber
  }


  if(!c(is.null(Vnumber)|is.na(Vnumber))){ message(paste("Exporting Version V", Vnumber))} else {

    if(c(is.null(Vnumber)|is.na(Vnumber)) & length(LF)!=0) {   message(paste("Exporting Version V2")) } else {
    message(paste("Exporting Version V1")) }}


  message(paste("Exporting", object, "objects"))

  NB.raw.clinic <- which(c(attributes(MetaData)$Data.Type=="Patient.Clinical.data" | attributes(MetaData)$Data.Type=="Samples.Clinical.data" ) & attributes(MetaData)$Raw.data=="Yes" )
  if(length(NB.raw.clinic)>0) {
    count <- count+1
    message("-------------------------------------------------")
    message(paste("Exporting", count, "/", object,"object: ","Raw.clinical data will not be exported"))
  }


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




     NB.Samples.Patients.pheno <-   which(c(attributes(MetaData)$Data.Type=="Patient.Clinical.data" |attributes(MetaData)$Data.Type=="Samples.Clinical.data" ) & attributes(MetaData)$Raw.data=="No" )

        if(length(NB.Samples.Patients.pheno)!=0) {

       for (j in NB.Samples.Patients.pheno) {

         LF <- list.files(list.files.path$Propject.VerifiedDataset)
         lengthSTR <- length(LF[str_detect(LF,names(MetaData)[j])])

         if(lengthSTR==0){
           count <- count+1
           message("-------------------------------------------------")
           message(paste("Exporting", count, "/", object,"object: ",names(MetaData)[j]))
           filename <- paste0(list.files.path$Propject.VerifiedDataset,"/",names(MetaData)[j],".csv")
           z <-  MetaData[[j]]


           if(is.na(Vnumber)| is.null(Vnumber)){
             write.table(z,row.names = F ,file = filename, sep = "\t")
             } else {

               if(is.na(Vnumber)) {Vnumber2 = 1}
               if(Vnumber2==1){ Vnumber2 = 2}
               if(!is.null(Vnumber) & !is.na(Vnumber)){Vnumber2 = Vnumber }

               filename <- paste0(list.files.path$Propject.VerifiedDataset,"/",names(MetaData)[j],".V", Vnumber2,".csv")
               write.table(z,row.names = F ,file = filename, sep = "\t")


             }

           attributes(MetaData)$Version <- "V1"

         } else{
           count <- count+1
           message("-------------------------------------------------")
           message(paste("Exporting", count, "/", object,"object: ",names(MetaData)[j]))
           message(paste(names(MetaData)[j],"exported file exist. Versionning data"))



           df <- file.info(list.files(list.files.path$Propject.VerifiedDataset, full.names = T))
           df$Filenames <- unlist(lapply(str_split(rownames(df),paste0(project,"/")),"[[",2))
           df <- df[str_detect(df$Filenames, names(MetaData)[j]),]

           filepath <- rownames(df)[which.max(df$mtime)]

           filename <-  unlist(lapply(str_split(filepath,paste0(project,"/")),"[[",2))
           filename2 <- unlist(lapply(str_split(filename,".csv"),"[[",1))
            extension <- unlist(lapply(str_split(filename,".csv"),"[[",2))




             if(is.na(Vnumber)) {
               Vnumber2 = 1
               filepath2 <- paste0(list.files.path$Propject.VerifiedDataset,"/",filename2,".V",Vnumber2,".csv",extension)
               file.rename(from = filepath, to = filepath2)


             }

             #adding version file history

             if(!is.null(Vnumber) & !is.na(Vnumber)){Vnumber2 = Vnumber }

             if(Vnumber2==1){ Vnumber2 = 2}

             attributes(MetaData)$Version <- paste0("V", Vnumber2)

             message(paste0("Exporting ", count, " / ", object,"object: ",names(MetaData)[j]))
             filename3 <- unlist(lapply(str_split(filename2,".V"),"[[",1))
             filepath2 <- paste0(list.files.path$Propject.VerifiedDataset,"/",filename3,".V",Vnumber2,".csv")
             z <- cbind( MetaData[[j]])
             write.table(z,row.names = F ,file = filepath2, sep = "\t")

         }


       } #for i in NBsamples

     } # if( length(NB.Samples.Patients.pheno)!=0)


     NB.Expression.Matrix <-  which(attributes(MetaData)$Data.Type=="Expression.Matrix")

        if(length(NB.Expression.Matrix)!=0) {

          for (j in NB.Expression.Matrix) {

        LF <- list.files(list.files.path$Propject.VerifiedDataset)
        lengthSTR <- length(LF[str_detect(LF,names(MetaData)[j])])

        if(lengthSTR==0){
          count <- count+1
          message("-------------------------------------------------")
          message(paste("Exporting", count, "/", object,"object: ",names(MetaData)[j]))
          if(attributes(MetaData)$Omics.type!="Single.Cell"){
            filename <- paste0(list.files.path$Propject.VerifiedDataset,"/",project,".",names(MetaData)[j], ".csv")
          }

          z <- cbind("GeneSymbol" = rownames(MetaData[[j]]), MetaData[[j]])


          if(attributes(MetaData)$Omics.type=="Single.Cell"){
            filename <- paste0(list.files.path$Propject.VerifiedDataset,"/",project,".",names(MetaData)[j],".mtx")
            filename.genes <- paste0(list.files.path$Propject.VerifiedDataset,"/",project,".GenesInfo.csv")
          }



          if(is.na(Vnumber)| is.null(Vnumber)){

            if(attributes(MetaData)$Omics.type!="Single.Cell") {write.table(z,row.names = F ,file = filename, sep = "\t")
              message(paste("Compressing"))
              R.utils::gzip(filename, destname=sprintf("%s.gz", filename), overwrite=T, remove=TRUE, BFR.SIZE=1e+07)}


            if(attributes(MetaData)$Omics.type=="Single.Cell") {

              count <- count+1
              message("-------------------------------------------------")
              message(paste("Exporting", count, "/", object,"object: ","GeneInfo", "file"))
              write.table(rownames(MetaData[[j]]),row.names = F ,file = filename.genes, sep = "\t")

              writeMM(Matrix(MetaData[[j]], sparse = T),file = filename)
              message(paste("Compressing"))
              R.utils::gzip(filename, destname=sprintf("%s.gz", filename), overwrite=T, remove=TRUE, BFR.SIZE=1e+07)

            }



              } else {

              if(is.na(Vnumber)) {Vnumber2 = 1}
              if(Vnumber2==1){ Vnumber2 = 2}
              if(!is.null(Vnumber) & !is.na(Vnumber)){Vnumber2 = Vnumber }

              filename <- paste0(list.files.path$Propject.VerifiedDataset,"/",project,".",names(MetaData)[j],".V", Vnumber2)

              if(attributes(MetaData)$Omics.type=="Single.Cell"){
                filename <- paste0(list.files.path$Propject.VerifiedDataset,"/",project,".",names(MetaData)[j],".V", Vnumber2, ".mtx")
                filename.genes <- paste0(list.files.path$Propject.VerifiedDataset,"/",project,".GenesInfo",".V", Vnumber2, ".csv")
              }




              if(attributes(MetaData)$Omics.type!="Single.Cell") {write.table(z,row.names = F ,file = paste0(filename,".csv"), sep = "\t")
                message(paste("Compressing"))
                R.utils::gzip(filename, destname=sprintf("%s.gz", filename), overwrite=T, remove=TRUE, BFR.SIZE=1e+07)}

              if(attributes(MetaData)$Omics.type=="Single.Cell") {
                write.table(rownames(MetaData[[j]]),row.names = F ,file = filename.genes, sep = "\t")
                writeMM(Matrix(MetaData[[j]], sparse = T),file = paste0(filename,".mtx"))
                message(paste("Compressing"))
                R.utils::gzip(filename, destname=sprintf("%s.gz", filename), overwrite=T, remove=TRUE, BFR.SIZE=1e+07)

                }





            }



        } else{
          count <- count+1
          message("-------------------------------------------------")
          message(paste("Exporting", count, "/", object,"object: ",names(MetaData)[j]))

          if( attributes(MetaData)$Raw.data[j]=="Yes"){ Raw = "Raw" } else { Raw = "Normalized"}
          message(paste(Raw, attributes(MetaData)$Data.Type[j],"exported file exist. Versionning data."))

          df <- file.info(list.files(list.files.path$Propject.VerifiedDataset, full.names = T))
          df$Filenames <- unlist(lapply(str_split(rownames(df),paste0(project,"/")),"[[",2))
          df <- df[str_detect(df$Filenames, names(MetaData)[j]),]
          filepath <- rownames(df)[which.max(df$mtime)]

          filename <-  unlist(lapply(str_split(filepath,paste0(project,"/")),"[[",2))

          if(attributes(MetaData)$Omics.type!="Single.Cell") {
            filename2 <- unlist(lapply(str_split(filename,".csv"),"[[",1))
            extension <- unlist(lapply(str_split(filename,".csv"),"[[",2))}

          if(attributes(MetaData)$Omics.type=="Single.Cell") {
            filename2 <- unlist(lapply(str_split(filename,".mtx"),"[[",1))
            extension <- unlist(lapply(str_split(filename,".mtx"),"[[",2))}




            if(is.na(Vnumber)) {
              Vnumber2 = 1

              if(attributes(MetaData)$Omics.type!="Single.Cell") {
                filepath2 <- paste0(list.files.path$Propject.VerifiedDataset,"/",filename2,".V",Vnumber2,".csv",extension)
                file.rename(from = filepath, to = filepath2)
              }

              if(attributes(MetaData)$Omics.type=="Single.Cell") {
                filepath2 <- paste0(list.files.path$Propject.VerifiedDataset,"/",filename2,".V",Vnumber2,".mtx",extension)
                file.rename(from = filepath, to = filepath2)


                message(paste(project, "Single cell GenesInfo exported file exist. Versionning data."))

                file.rename(from =  paste0(list.files.path$Propject.VerifiedDataset,"/",project,".GenesInfo.csv"),
                            to = paste0(list.files.path$Propject.VerifiedDataset,"/",project,".GenesInfo",".V", Vnumber2, ".csv") )


              }


            }

            #adding version file history
            if(!is.null(Vnumber) & !is.na(Vnumber)){Vnumber2 = Vnumber }
            if(Vnumber2==1){ Vnumber2 = 2}

            message(paste0("Exporting ", count, " / ", object," object: ",names(MetaData)[j]))
            filename3 <- unlist(lapply(str_split(filename2,".V"),"[[",1))
                   filepath2 <- paste0(list.files.path$Propject.VerifiedDataset,"/",filename3,".V",Vnumber2,".csv")
            z <- cbind("GeneSymbol" = rownames(MetaData[[j]]), MetaData[[j]])


            if(attributes(MetaData)$Omics.type=="Single.Cell"){
              filepath2 <- paste0(list.files.path$Propject.VerifiedDataset,"/",project,".",names(MetaData)[j],".V", Vnumber2, ".mtx")
              filename.genes <- paste0(list.files.path$Propject.VerifiedDataset,"/",project,".",names(MetaData)[j],".GenesInfo",".V", Vnumber2, ".csv")
            }



            if(attributes(MetaData)$Omics.type!="Single.Cell") { write.table(z,row.names = F ,file = filepath2,sep = "\t")

              message(paste("Compressing"))
              R.utils::gzip(filepath2, destname=sprintf("%s.gz", filepath2), overwrite=T, remove=TRUE, BFR.SIZE=1e+07)}

            if(attributes(MetaData)$Omics.type=="Single.Cell") {
            count <- count+1
            message("-------------------------------------------------")
            message(paste("Exporting", count, "/", object,"object: ","GeneInfo", "file"))
            write.table(rownames(MetaData[[j]]),row.names = F ,file = filename.genes, sep = "\t")
            message("-------------------------------------------------")

            writeMM(Matrix(MetaData[[j]], sparse = T),file = filepath2)

              message(paste("Compressing"))
              R.utils::gzip(filepath2, destname=sprintf("%s.gz", filepath2), overwrite=T, remove=TRUE, BFR.SIZE=1e+07)

            }




        }

          }
      } # if(length(NB.Expression.Matrix)!=0)


     NB.geneAnnot<-  which(attributes(MetaData)$Data.Type=="geneAnnotation.file")

         if(length(NB.Expression.Matrix)!=0) {

            for (j in NB.geneAnnot) {

      LF <- list.files(list.files.path$Propject.VerifiedDataset)
      lengthSTR <- length(LF[str_detect(LF,names(MetaData)[j])])

      if(lengthSTR==0){
        count <- count+1
        message("-------------------------------------------------")
        message(paste("Exporting", count, "/", object,"object: ",names(MetaData)[j], "file"))
        filename <- paste0(list.files.path$Propject.VerifiedDataset,"/",project,".",names(MetaData)[j],".csv")
        z <- cbind(MetaData[[j]])

        if(is.na(Vnumber)| is.null(Vnumber)){
        write.table(z,row.names = F ,file = filename, sep = "\t")
        message(paste("Compressing"))
        R.utils::gzip(filename, destname=sprintf("%s.gz", filename), overwrite=T, remove=TRUE, BFR.SIZE=1e+07)} else {

          if(is.na(Vnumber)) {Vnumber2 = 1}
          if(Vnumber2==1){ Vnumber2 = 2}
          if(!is.null(Vnumber) & !is.na(Vnumber)){Vnumber2 = Vnumber }

          filename <- paste0(list.files.path$Propject.VerifiedDataset,"/",project,".",names(MetaData)[j],".V", Vnumber2,".csv")
          write.table(z,row.names = F ,file = filename, sep = "\t")
          message(paste("Compressing"))
          R.utils::gzip(filename, destname=sprintf("%s.gz", filename), overwrite=T, remove=TRUE, BFR.SIZE=1e+07)

          }




      } else{
        count <- count+1
        message("-------------------------------------------------")
        message(paste("Exporting", count, "/", object,"object: ",names(MetaData)[j], "file"))

        message(paste(attributes(MetaData)$Data.Type[j],"exported file exist. Versionning data."))

        df <- file.info(list.files(list.files.path$Propject.VerifiedDataset, full.names = T))
        df$Filenames <- unlist(lapply(str_split(rownames(df),paste0(project,"/")),"[[",2))
        df <- df[str_detect(df$Filenames, names(MetaData)[j]),]
        filepath <- rownames(df)[which.max(df$mtime)]


        filename <-  unlist(lapply(str_split(filepath,paste0(project,"/")),"[[",2))


          filename2 <- unlist(lapply(str_split(filename,".csv"),"[[",1))
          extension <- unlist(lapply(str_split(filename,".csv"),"[[",2))


          if(is.na(Vnumber)) {
            Vnumber2 = 1

            filepath2 <- paste0(list.files.path$Propject.VerifiedDataset,"/",filename2,".V",Vnumber2,".csv",extension)
            file.rename(from = filepath, to = filepath2)
          }

          #adding version file history
          if(!is.null(Vnumber) & !is.na(Vnumber)){Vnumber2 = Vnumber }
          if(Vnumber2==1){ Vnumber2 = 2}

          message(paste0("Exporting ", count, " / ", object,"object: ",names(MetaData)[j]))
          filename3 <- unlist(lapply(str_split(filename2,".V"),"[[",1))
          filepath2 <- paste0(list.files.path$Propject.VerifiedDataset,"/",filename3,".V",Vnumber2,".csv")
          z <- cbind( MetaData[[j]])
          write.table(z,row.names = F ,file = filepath2, sep = "\t")
          message(paste("Compressing"))
          R.utils::gzip(filepath2, destname=sprintf("%s.gz", filepath2), overwrite=T, remove=TRUE, BFR.SIZE=1e+07)



      }}
           }# if length(NB.geneAnnot)

     objs =  mget(ls(envir=.GlobalEnv), envir=.GlobalEnv)
     NO <- names(Filter(function(i) inherits(i, "list"), objs))[str_detect(toupper(names(Filter(function(i) inherits(i, "list"), objs))),"META")]

     pos <- 1
     envir = as.environment(pos)

     assign(NO, MetaData, envir = envir)




  }#function


