#' ExportCSV Export MetaData inside object into ".csv" files
#'
#' @param MetaData a MetaData  data files
#' @param list.files.path dirpath
#' @return ".csv" files into working directory
#' @export
#' @import utils
#' @import R.utils
#' @examples "non"
ExportCSV <- function (MetaData, list.files.path, project){

  if(is.null(MetaData)){stop("Need a MetaData List file")}
  if(!is.list(MetaData)){stop("Need a MetaData List file")}

  if(is.null(list.files.path)){stop("Need a list file path for saving data")}
  if(!is.list(list.files.path)){stop(paste("list.files.path must be a list of file path whith Script, Raw genomic, Raw clinic, Processed and References directories in Parent Directory." ))}

  count <- 0
  object <- length(MetaData)+2
  name <- names(MetaData)

  message(paste("Exporting", object, "objects"))

  if(!all(str_detect(toupper(names(MetaData)), "RAW.CLINIC")==F)) {
    count <- count+1
    message("-------------------------------------------------")
    message(paste("Exporting", count, "/", object,"object: ","Raw.clinical data will not be exported"))
  }


  if(exists("LexicClinic", mode= "any" )) {
    count <- count+1
    message("-------------------------------------------------")
    message(paste("Exporting", count, "/", object,"object: ","LexicClinic"))

    if (file.exists(paste0(list.files.path$Project.Processes,"/",project,".Lexic.txt"))) {
      #Delete file if it exists
      file.remove(paste0(list.files.path$Project.Processes,"/",project,".Lexic.txt"))
    }

    LexicClinic <- lapply(LexicClinic, function(x) {c(x[1],x)}) #Mandatory to duplicated listName in the listed values.
    lapply(LexicClinic, write, paste0(list.files.path$Project.Processes,"/",project,".Lexic.txt"), append=TRUE, ncolumns=1000 ) #write a ".txt" file without listNames
  }


  if(exists("SamplesOrPatients", mode= "any" )) {

    count <- count+1
    message("-------------------------------------------------")
    message(paste("Exporting", count, "/", object,"object: ","SamplesOrPatients"))

    write.table(SamplesOrPatients, file= paste0(list.files.path$Project.Processes,"/",project,"SamplesOrPatients.txt"),quote = F, sep = "\t", row.names = F)}



     NB.Samples.Patients.pheno <-  which(attributes(MetaData)$Data.Type=="Clinical.data" & attributes(MetaData)$Raw.data=="No" )

        if(length(NB.Samples.Patients.pheno)!=0) {

       for (j in NB.Samples.Patients.pheno) {

         LF <- list.files(list.files.path$Propject.VerifiedDataset)
         lengthSTR <- length(LF[str_detect(LF,names(MetaData)[j])])

         if(lengthSTR==0){
           count <- count+1
           message("-------------------------------------------------")
           message(paste("Exporting", count, "/", object,"object: ",names(MetaData)[j]))
           filename <- paste0(list.files.path$Propject.VerifiedDataset,"/",project,".",names(MetaData)[j],".csv")
           z <- cbind( MetaData[[j]])
           write.csv(z,row.names = F ,file = filename)
           message(paste("Compressing"))
           R.utils::gzip(filename, destname=sprintf("%s.gz", filename), overwrite=T, remove=TRUE, BFR.SIZE=1e+07)

         } else{
           count <- count+1
           message("-------------------------------------------------")
           message(paste("Exporting", count, "/", object,"object: ",names(MetaData)[j]))
           message("Samples' pheno exported file exist. Loading to compare if different to Metaobject Version")
           df <- file.info(list.files(list.files.path$Propject.VerifiedDataset, full.names = T))
           df$Filenames <- unlist(lapply(str_split(rownames(df),paste0(project,"/")),"[[",2))
           df <- df[str_detect(df$Filenames, names(MetaData)[j]),]
           filepath <- rownames(df)[which.max(df$mtime)]

           loaded.file <- as.data.frame(data.table::fread(filepath, na.strings = "NA"))
           loaded.file <- as.data.frame(loaded.file)
           loaded.file <- apply(loaded.file, MARGIN = 2,function(x) as.character(x))
           loaded.file <- as.data.frame(loaded.file)
           rownames(loaded.file) <- loaded.file[,1]

           source.file <- apply(MetaData[[j]], MARGIN = 2,function(x) as.character(x))
           source.file <- as.data.frame(source.file)
           rownames(source.file) <- source.file[,1]


           filename <-  unlist(lapply(str_split(filepath,paste0(project,"/")),"[[",2))

           if(!dplyr::all_equal(loaded.file,source.file)==T) { #loaded.file ans MetaData[[i]] are not the same data.frame
             message(paste(names(MetaData)[j], "data are not the same original data. A new version will be created"))
             filename2 <- unlist(lapply(str_split(filename,".csv"),"[[",1))
             extension <- unlist(lapply(str_split(filename,".csv"),"[[",2))
             version <- str_extract(filename2,"V[0-9]*")
             Vnumber <- as.numeric(str_extract(version,"([0-9]+).*$"))+1

             if(is.na(Vnumber)) {
               Vnumber = 1
               message(paste(filename, " already exist. Adding Version to the latest exportedfile"))
               filepath2 <- paste0(list.files.path$Propject.VerifiedDataset,"/",filename2,".V",Vnumber,".csv",extension)
               file.rename(from = filepath, to = filepath2)
             }

             #adding version file history

             if(Vnumber==1){ Vnumber = 2}

             message(paste0("Exporting ", count, " / ", object,"object: ",names(MetaData)[j], "V", Vnumber))
             filename3 <- unlist(lapply(str_split(filename2,".V"),"[[",1))
             filepath2 <- paste0(list.files.path$Propject.VerifiedDataset,"/",filename3,".V",Vnumber,".csv")
             z <- cbind( MetaData[[j]])
             write.csv(z,row.names = F ,file = filepath2)
             message(paste("Compressing"))
             R.utils::gzip(filepath2, destname=sprintf("%s.gz", filepath2), overwrite=T, remove=TRUE, BFR.SIZE=1e+07)
           } else {


             message(names(MetaData)[j]," is the same data frame as ",filename)
             message("Will not be saved")}
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
          filename <- paste0(list.files.path$Propject.VerifiedDataset,"/",project,".",names(MetaData)[j],".csv")
          z <- cbind("GeneSymbol" = rownames(MetaData[[j]]), MetaData[[j]])
          write.csv(z,row.names = F ,file = filename)
          message(paste("Compressing"))
          R.utils::gzip(filename, destname=sprintf("%s.gz", filename), overwrite=T, remove=TRUE, BFR.SIZE=1e+07)

        } else{
          count <- count+1
          message("-------------------------------------------------")
          message(paste("Exporting", count, "/", object,"object: ",names(MetaData)[j]))

          if( attributes(MetaData)$Raw.data[j]=="Yes"){ Raw = "Raw" } else { Raw = "Normalized"}
          message(paste(Raw, attributes(MetaData)$Data.Type[j],"exported file exist. Loading to compare if different to Metaobject Version"))

          df <- file.info(list.files(list.files.path$Propject.VerifiedDataset, full.names = T))
          df$Filenames <- unlist(lapply(str_split(rownames(df),paste0(project,"/")),"[[",2))
          df <- df[str_detect(df$Filenames, names(MetaData)[j]),]
          filepath <- rownames(df)[which.max(df$mtime)]

          loaded.file <- as.data.frame(data.table::fread(filepath, na.strings = "NA"))
          loaded.file <- as.data.frame(loaded.file)
          rownames(loaded.file) <- loaded.file[,1]
          loaded.file <- loaded.file[,-1]

          source.file <- as.data.frame(MetaData[[j]])

          filename <-  unlist(lapply(str_split(filepath,paste0(project,"/")),"[[",2))

          if(!dplyr::all_equal(loaded.file,source.file)==T) { #loaded.file ans MetaData[[i]] are not the same data.frame
            message(paste(names(MetaData)[j], "data are not the same original data. A new version will be created"))
            filename2 <- unlist(lapply(str_split(filename,".csv"),"[[",1))
            extension <- unlist(lapply(str_split(filename,".csv"),"[[",2))
            version <- str_extract(filename2,"V[0-9]*")
            Vnumber <- as.numeric(str_extract(version,"([0-9]+).*$"))+1

            if(is.na(Vnumber)) {
              Vnumber = 1
              message(paste(filename, " already exist. Adding Version to the latest exportedfile"))
              filepath2 <- paste0(list.files.path$Propject.VerifiedDataset,"/",filename2,".V",Vnumber,".csv",extension)
              file.rename(from = filepath, to = filepath2)
            }

            #adding version file history

            if(Vnumber==1){ Vnumber = 2}

            message(paste0("Exporting ", count, " / ", object,"object: ",names(MetaData)[j], "V", Vnumber))
            filename3 <- unlist(lapply(str_split(filename2,".V"),"[[",1))
            filepath2 <- paste0(list.files.path$Propject.VerifiedDataset,"/",filename3,".V",Vnumber,".csv")
            z <- cbind("GeneSymbol" = rownames(MetaData[[j]]), MetaData[[j]])
            write.csv(z,row.names = F ,file = filepath2)
            message(paste("Compressing"))
            R.utils::gzip(filepath2, destname=sprintf("%s.gz", filepath2), overwrite=T, remove=TRUE, BFR.SIZE=1e+07)
          } else {


            message(names(MetaData)[j]," is the same data frame as ",filename)
            message("Will not be saved")}
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
        write.csv(z,row.names = F ,file = filename)
        message(paste("Compressing"))
        R.utils::gzip(filename, destname=sprintf("%s.gz", filename), overwrite=T, remove=TRUE, BFR.SIZE=1e+07)

      } else{
        count <- count+1
        message("-------------------------------------------------")
        message(paste("Exporting", count, "/", object,"object: ",names(MetaData)[j], "file"))

        message(paste(attributes(MetaData)$Data.Type[j],"exported file exist. Loading to compare if different to Metaobject Version"))

        df <- file.info(list.files(list.files.path$Propject.VerifiedDataset, full.names = T))
        df$Filenames <- unlist(lapply(str_split(rownames(df),paste0(project,"/")),"[[",2))
        df <- df[str_detect(df$Filenames, names(MetaData)[j]),]
        filepath <- rownames(df)[which.max(df$mtime)]

        loaded.file <- as.data.frame(data.table::fread(filepath, na.strings = "NA"))
        loaded.file <- as.data.frame(loaded.file)
        loaded.file <- apply(loaded.file, MARGIN  = 2, FUN = function(x) as.character(x))
        loaded.file <- as.data.frame(loaded.file)

        if(length(which(colnames(loaded.file)=="EnsemblID"))==0){     rownames(loaded.file) <- loaded.file[,1]

          } else {        rownames(loaded.file) <- loaded.file[,"EnsemblID"]   }




        source.file <- as.data.frame(MetaData[[j]])
        source.file <- apply(source.file, MARGIN = 2, FUN = function(x) as.character(x))
        source.file <- as.data.frame(source.file)

        if(length(which(colnames(source.file)=="EnsemblID"))==0){     rownames(source.file) <- source.file[,1]

        } else {        rownames(source.file) <- source.file[,"EnsemblID"]   }


        source.file$attributes <- NULL
        loaded.file$attributes <- NULL

        source.file$Entrez.id <- NULL
        loaded.file$Entrez.id <- NULL

        filename <-  unlist(lapply(str_split(filepath,paste0(project,"/")),"[[",2))

        if(!dplyr::all_equal(loaded.file,source.file, ignore_col_order = T, ignore_row_order = T)==T) { #loaded.file ans MetaData[[i]] are not the same data.frame
          message(paste(names(MetaData)[j], "data are not the same original data. A new version will be created"))
          filename2 <- unlist(lapply(str_split(filename,".csv"),"[[",1))
          extension <- unlist(lapply(str_split(filename,".csv"),"[[",2))
          version <- str_extract(filename2,"V[0-9]*")
          Vnumber <- as.numeric(str_extract(version,"([0-9]+).*$"))+1

          if(is.na(Vnumber)) {
            Vnumber = 1
            message(paste(filename, " already exist. Adding Version to the latest exportedfile"))
            filepath2 <- paste0(list.files.path$Propject.VerifiedDataset,"/",filename2,".V",Vnumber,".csv",extension)
            file.rename(from = filepath, to = filepath2)
          }

          #adding version file history

          if(Vnumber==1){ Vnumber = 2}

          message(paste0("Exporting ", count, " / ", object,"object: ",names(MetaData)[j], "V", Vnumber))
          filename3 <- unlist(lapply(str_split(filename2,".V"),"[[",1))
          filepath2 <- paste0(list.files.path$Propject.VerifiedDataset,"/",filename3,".V",Vnumber,".csv")
          z <- cbind( MetaData[[j]])
          write.csv(z,row.names = F ,file = filepath2)
          message(paste("Compressing"))
          R.utils::gzip(filepath2, destname=sprintf("%s.gz", filepath2), overwrite=T, remove=TRUE, BFR.SIZE=1e+07)
        } else {


          message(names(MetaData)[j]," is the same data frame as ",filename)
          message("Will not be saved")}
      }


    }} # if length(NB.geneAnnot)


  }#function


