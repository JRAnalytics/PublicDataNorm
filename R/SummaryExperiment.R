SummaryExperiment <- function(MetaData) {

  if(is.null(MetaData)){stop("Need a MetaData List file")}
  if(!is.list(MetaData)){stop("Need a MetaData List file")}


  name <- names(MetaData)

  if(!file.exists("Readme.txt")){

    tme <- Sys.Date()
    tme <- format(tme, format="%B %d %Y")

    dt <- data.frame("Type"="File created the : " ,"Description"=tme)
    sp <- data.frame("Type"="---------:" ,"Description"="---------")

    dt <- rbind(dt,data.frame("Type" = "Files included in folder : ","Description" = "-"),sp)



    name <- names(MetaData)

    for (i in name){


      if(str_detect(i, c("matrix"))==T){
        nr <- nrow(MetaData[[i]])
        nc <- ncol(MetaData[[i]])



        if(str_detect(i, "Raw")){ Assay = "Raw counts"}
        if(str_detect(i, "TPM")){ Assay="TPM normalization"}
        if(str_detect(i, "FKPM")){ Assay="FKPM normalization"}
        if(str_detect(i, "Normalized")){ Assay="Normalized gene expression"}

        ltest <- data.frame("Type"=c("File: " , "Class: " ,"Dimension: ", "Assay: ", "Rownames: ","Colnames: "),"Description" = c(
          paste0(i,".csv"),
          class(MetaData[[i]]),
          paste(nr,"x",nc),
          Assay,
          paste(rownames(MetaData[[i]])[1],"...",rownames(MetaData[[i]])[nrow(MetaData[[i]])]),
          paste(colnames(MetaData[[i]])[2],"...",colnames(MetaData[[i]])[ncol(MetaData[[i]])])
        ))


        dt <- rbind(dt, ltest,sp)

      }

      if(str_detect(i, c("gene"))==T){
        nr <- nrow(MetaData[[i]])
        nc <- ncol(MetaData[[i]])



        if(str_detect(i, "Annotation")){ Assay ="Gene annotation"}


        ltest <- data.frame("Type"=c("File: " , "Class: " ,"Dimension: ", "Assay: ", "Rownames: ","Colnames: "),"Description" = c(
          paste0(i,".csv"),
          class(MetaData[[i]]),
          paste(nr,"x",nc),
          Assay,
          paste(rownames(MetaData[[i]])[1],"...",rownames(MetaData[[i]])[nrow(MetaData[[i]])]),
          paste(colnames(MetaData[[i]])[2],"...",colnames(MetaData[[i]])[ncol(MetaData[[i]])])
        ))



        dt <- rbind(dt, ltest,sp)


      }

      if(str_detect(i, c("clinic"))==T){
        nr <- nrow(MetaData[[i]])
        nc <- ncol(MetaData[[i]])



        if(str_detect(i, "clinic_")){ Assay =  "Original clinical data"}
        if(str_detect(i, "Patient")){ Assay =  "Patient's clinical data"}
        if(str_detect(i, "Sample")){ Assay = "Samples pathological records"}

        ltest <- data.frame("Type"=c("File: " , "Class: " ,"Dimension: ", "Assay: ", "Rownames: ","Colnames: "),"Description" = c(
          paste0(i,".csv"),
          class(MetaData[[i]]),
          paste(nr,"x",nc),
          Assay,
          paste(rownames(MetaData[[i]])[1],"...",rownames(MetaData[[i]])[nrow(MetaData[[i]])]),
          paste(colnames(MetaData[[i]])[2],"...",colnames(MetaData[[i]])[ncol(MetaData[[i]])])
        ))



        dt <- rbind(dt, ltest,sp)


      }

      if(str_detect(i, c("pheno"))==T){
        nr <- nrow(MetaData[[i]])
        nc <- ncol(MetaData[[i]])


        if(str_detect(i, "clinic_")){ Assay =  "Original clinical data"}
        if(str_detect(i, "Patient")){ Assay =  "Patient's clinical data"}
        if(str_detect(i, "Sample")){ Assay = "Samples pathological records"}

        ltest <- data.frame("Type"=c("File: " , "Class: " ,"Dimension: ", "Assay: ", "Rownames: ","Colnames: "),"Description" = c(
          paste0(i,".csv"),
          class(MetaData[[i]]),
          paste(nr,"x",nc),
          Assay,
          paste(rownames(MetaData[[i]])[1],"...",rownames(MetaData[[i]])[nrow(MetaData[[i]])]),
          paste(colnames(MetaData[[i]])[2],"...",colnames(MetaData[[i]])[ncol(MetaData[[i]])])
        ))



        dt <- rbind(dt, ltest,sp)

      }

    }

    write.table(dt,"Readme.txt",row.names = F, quote = FALSE)

    } else {



      sp <- data.frame("Type"="---------:" ,"Description"="---------")

      tme <- Sys.Date()
      tme <- format(tme, format="%B %d %Y")

      mod <- data.frame("Type"="MetaData modified the: " ,"Description"=tme)
      dt <- rbind(sp,mod,sp)

      for (i in name){


        if(str_detect(i, c("matrix"))==T){
          nr <- nrow(MetaData[[i]])
          nc <- ncol(MetaData[[i]])



          if(str_detect(i, "Raw")){ Assay = "Raw counts"}
          if(str_detect(i, "TPM")){ Assay="TPM normalization"}
          if(str_detect(i, "FKPM")){ Assay="FKPM normalization"}
          if(str_detect(i, "Normalized")){ Assay="Normalized gene expression"}

          ltest <- data.frame("Type"=c("File: " , "Class: " ,"Dimension: ", "Assay: ", "Rownames: ","Colnames: "),"Description" = c(
            paste0(i,".csv"),
            class(MetaData[[i]]),
            paste(nr,"x",nc),
            Assay,
            paste(rownames(MetaData[[i]])[1],"...",rownames(MetaData[[i]])[nrow(MetaData[[i]])]),
            paste(colnames(MetaData[[i]])[2],"...",colnames(MetaData[[i]])[ncol(MetaData[[i]])])
          ))


          dt <- rbind(dt, ltest,sp)

        }

        if(str_detect(i, c("gene"))==T){
          nr <- nrow(MetaData[[i]])
          nc <- ncol(MetaData[[i]])



          if(str_detect(i, "Annotation")){ Assay ="Gene annotation"}


          ltest <- data.frame("Type"=c("File: " , "Class: " ,"Dimension: ", "Assay: ", "Rownames: ","Colnames: "),"Description" = c(
            paste0(i,".csv"),
            class(MetaData[[i]]),
            paste(nr,"x",nc),
            Assay,
            paste(rownames(MetaData[[i]])[1],"...",rownames(MetaData[[i]])[nrow(MetaData[[i]])]),
            paste(colnames(MetaData[[i]])[2],"...",colnames(MetaData[[i]])[ncol(MetaData[[i]])])
          ))



          dt <- rbind(dt, ltest,sp)


        }

        if(str_detect(i, c("clinic"))==T){
          nr <- nrow(MetaData[[i]])
          nc <- ncol(MetaData[[i]])



          if(str_detect(i, "clinic_")){ Assay =  "Original clinical data"}
          if(str_detect(i, "Patient")){ Assay =  "Patient's clinical data"}
          if(str_detect(i, "Sample")){ Assay = "Samples pathological records"}

          ltest <- data.frame("Type"=c("File: " , "Class: " ,"Dimension: ", "Assay: ", "Rownames: ","Colnames: "),"Description" = c(
            paste0(i,".csv"),
            class(MetaData[[i]]),
            paste(nr,"x",nc),
            Assay,
            paste(rownames(MetaData[[i]])[1],"...",rownames(MetaData[[i]])[nrow(MetaData[[i]])]),
            paste(colnames(MetaData[[i]])[2],"...",colnames(MetaData[[i]])[ncol(MetaData[[i]])])
          ))



          dt <- rbind(dt, ltest,sp)


        }

        if(str_detect(i, c("pheno"))==T){
          nr <- nrow(MetaData[[i]])
          nc <- ncol(MetaData[[i]])


          if(str_detect(i, "clinic_")){ Assay =  "Original clinical data"}
          if(str_detect(i, "Patient")){ Assay =  "Patient's clinical data"}
          if(str_detect(i, "Sample")){ Assay = "Samples pathological records"}

          ltest <- data.frame("Type"=c("File: " , "Class: " ,"Dimension: ", "Assay: ", "Rownames: ","Colnames: "),"Description" = c(
            paste0(i,".csv"),
            class(MetaData[[i]]),
            paste(nr,"x",nc),
            Assay,
            paste(rownames(MetaData[[i]])[1],"...",rownames(MetaData[[i]])[nrow(MetaData[[i]])]),
            paste(colnames(MetaData[[i]])[2],"...",colnames(MetaData[[i]])[ncol(MetaData[[i]])])
          ))



          dt <- rbind(dt, ltest,sp)

        }

      }

      write.table(dt,"Readme.txt",row.names = F,quote = FALSE, append = T)
    }




    file.show("Readme.txt")
    closeAllConnections()

}




