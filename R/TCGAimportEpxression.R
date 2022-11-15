#' TCGAimportEpxression allow to import and build the rnaseq data from the TCGA
#'
#' @param project a character of TCGA project name
#' @param data.category a character of data categorie
#' @param data.type a character of data type
#' @param experimental.strategy a character of experimental.strategy
#' @param sample.type a character of sample.type
#' @param data.norm a character "Raw"  = raw count, "TPM" = tpm normalisation gene expression, "FKPM" = fkpm normalisation gene expression if it exist
#' @param remove.files.prepared FALSE or TRUE to remove data downloaded
#' @param platform platform
#' @param file.type filetype
#' @param legacy legacy T or F
#' @param query query object exported
#' @param bcr_patient_barcode grouped by
#' @param Meta a Meta object. If NULL, will create the Meta object. If exist, will add objects from queries.
#' @import data.table
#' @import TCGAbiolinks
#' @import dplyr
#' @import stringr
#' @import DT
#' @return a data frame with build rna expression raw or normalized
#' @export
#'
#' @examples "none"
#'
TCGAimportEpxression <- function(Meta=NULL,
                                 project =c("TCGA-PAAD"),
                                   data.category = "Transcriptome Profiling",
                                   data.type = "Gene Expression Quantification",
                                   experimental.strategy = "RNA-Seq",
                                   sample.type = "Primary Tumor",
                                    platform = NULL,
                                    file.type = NULL,
                                    legacy = F,
                                   data.norm = c("Raw", "TPM", "FKPM"),remove.files.prepared = F){

  require(dplyr)
  require(data.table)
  require(TCGAbiolinks)
  require(DT)
  require(stringr)
if(isServeOK()==FALSE){stop("Connection to server GDC failled")}


message("Querying", paste(project, "Meta"))

query <<- GDCquery(
  project = project ,
  data.category = data.category,
  data.type = data.type,
  legacy = legacy,
  experimental.strategy = experimental.strategy,
  platform = platform,
  file.type = file.type,
  sample.type = sample.type
)







query.clinic <- GDCquery(project = project,
                         data.category = "Clinical",
                         data.type = "Clinical Supplement",
                         data.format = "BCR Biotab")


source <- ifelse(query$legacy,"legacy","harmonized")


files <- file.path(
  query$results[[1]]$project, source,
  gsub(" ","_",query$results[[1]]$data_category),
  gsub(" ","_",query$results[[1]]$data_type),
  gsub(" ","_",query$results[[1]]$file_id),
  gsub(" ","_",query$results[[1]]$file_name)
)

if(!dir.exists(file.path("GDCdata", project, source, data.norm))){
  message("No ", paste(project, "data found."))
  message("Starting", paste(project, "loading. May take a few times"))
GDCdownload(
  query = query,
  method = "api",
  files.per.chunk = 10)
}

if(dir.exists(file.path("GDCdata",project, source))){message(paste(project, "data found"))}

files <- file.path("GDCdata", files)


cases <- ifelse(grepl("TCGA|TARGET",query$results[[1]]$project %>% unlist()),query$results[[1]]$cases,query$results[[1]]$sample.submitter_id)


message("Building gene expression data")

data <- TCGA.build(
  ask = data.norm,legacy = legacy,
  files = files,
  cases = cases,
  genome = ifelse(query$legacy,"hg19","hg38"),
  experimental.strategy = unique(query$results[[1]]$experimental_strategy)

)

data <- data[-which(duplicated(unlist(lapply(str_split(rownames(data),"[|]" ), "[[", 1)))),]

rownames(data) <-  unlist(lapply(str_split(rownames(data),"[|]" ), "[[", 1))



if(remove.files.prepared){message("Removing downloaded data")
  # removes files and empty directories
  removefilesrecursively(files)
}




if(!dir.exists(file.path("GDCdata",project,source,"Clinical"))){

  message(paste(project, "clinical data  not found"))
  message(paste("Downloading", project, "clinical data."))
  GDCdownload(query.clinic)
  }




if(is.null(Meta)){
if(dir.exists(file.path("GDCdata",project,source,"Clinical"))){

  message(paste(project, "clinical data found"))
  message(paste("Building",project, "clinical data"))}



clinic2 <- GDCprepare(query.clinic)


for (j in 1:length(clinic2)){

  for (i in 1:ncol(clinic2[[j]])){

    zz <- which(clinic2[[j]][,i]=="[Discrepancy]"| clinic2[[j]][,i]=="[Unknown]" |clinic2[[j]][,i]=="[Not Available]" |clinic2[[j]][,i]=="[Not Applicable]" | clinic2[[j]][,i]=="[Not Evaluated]")
    clinic2[[j]][zz,i] <- NA
  }

}




clinic3 <- as.data.frame(clinic2[[1]])
for (i in 2:length(clinic2)) {
  dup <- c(colnames(clinic3),colnames(clinic2[[i]]))[which(duplicated(c(colnames(clinic3),colnames(clinic2[[i]]))))]



  alone <- dup[! dup%in% c("bcr_patient_barcode","bcr_patient_uuid")]




  if(length(alone)==0){ clinic3 <- full_join(clinic3,clinic2[[i]], by =c("bcr_patient_barcode","bcr_patient_uuid")) } else if(ncol(as.data.frame(clinic2[[i]])[,dup])==ncol(as.data.frame(clinic2[[i]]))){

    clinic3 <- bind_rows(clinic3[,],clinic2[[i]][-c(1:2),])

  } else if(length(alone)>0) { clinic3 <- full_join(clinic3,clinic2[[i]], by =c("bcr_patient_barcode","bcr_patient_uuid", alone))



  } else { clinic3 <- full_join(clinic3[,!colnames(clinic3)%in%alone],clinic2[[i]][,!colnames(clinic2[[i]])%in%alone], by =c("bcr_patient_barcode","bcr_patient_uuid")) }

}


clinic3 <- clinic3[order(clinic3$bcr_patient_barcode),]
clinic3 <- clinic3[-c(1:2),]

clinic3_rolled <- clinic3 %>%

  # create groups by name
  group_by(bcr_patient_barcode) %>%

  summarise(across(everything(), ~paste0((na.omit(.x)), collapse = ";")))
clinic3_rolled <- as.data.frame(clinic3_rolled)
rownames(clinic3_rolled) <- clinic3_rolled[,"bcr_patient_barcode"]




Meta <- list("DF" = data,
             "clinic"= clinic3_rolled)

names(Meta) <- c(paste0(data.norm,"-", project,"-matrix"), paste0("clinic_",project))

} else {

  l <- length(Meta)
  name <- names(Meta)
  Meta[[l+1]] <- data
  names(Meta) <- c(name,paste0(data.norm,"-", project,"-matrix"))

}




if(!file.exists("Readme.txt")){

  tme <- Sys.Date()
  tme <- format(tme, format="%B %d %Y")

  dt <- data.frame("Type"="File created the : " ,"Description"=tme)
  sp <- data.frame("Type"="---------:" ,"Description"="---------")

  dt <- rbind(dt,data.frame("Type" = "Files included in folder : ","Description" = "-"),sp)



  name <- names(Meta)

  for (i in name){


    if(str_detect(i, c("matrix"))==T){
      nr <- nrow(Meta[[i]])
      nc <- ncol(Meta[[i]])



      if(str_detect(i, "Raw")){ Assay = "Raw counts"}
      if(str_detect(i, "TPM")){ Assay="TPM normalization"}
      if(str_detect(i, "FKPM")){ Assay="FKPM normalization"}
      if(str_detect(i, "Normalized")){ Assay="Normalized gene expression"}

      ltest <- data.frame("Type"=c("File: " , "Class: " ,"Dimension: ", "Assay: ", "Rownames: ","Colnames: "),"Description" = c(
        paste0(i,".csv"),
        class(Meta[[i]]),
        paste(nr,"x",nc),
        Assay,
        paste(rownames(Meta[[i]])[1],"...",rownames(Meta[[i]])[nrow(Meta[[i]])]),
        paste(colnames(Meta[[i]])[2],"...",colnames(Meta[[i]])[ncol(Meta[[i]])])
      ))


      dt <- rbind(dt, ltest,sp)

    }

    if(str_detect(i, c("gene"))==T){
      nr <- nrow(Meta[[i]])
      nc <- ncol(Meta[[i]])



      if(str_detect(i, "Annotation")){ Assay ="Gene annotation"}


      ltest <- data.frame("Type"=c("File: " , "Class: " ,"Dimension: ", "Assay: ", "Rownames: ","Colnames: "),"Description" = c(
        paste0(i,".csv"),
        class(Meta[[i]]),
        paste(nr,"x",nc),
        Assay,
        paste(rownames(Meta[[i]])[1],"...",rownames(Meta[[i]])[nrow(Meta[[i]])]),
        paste(colnames(Meta[[i]])[2],"...",colnames(Meta[[i]])[ncol(Meta[[i]])])
      ))



      dt <- rbind(dt, ltest,sp)


    }

    if(str_detect(i, c("clinic"))==T){
      nr <- nrow(Meta[[i]])
      nc <- ncol(Meta[[i]])



      if(str_detect(i, "clinic_")){ Assay =  "Original clinical data"}
      if(str_detect(i, "Patient")){ Assay =  "Patient's clinical data"}
      if(str_detect(i, "Sample")){ Assay = "Samples pathological records"}

      ltest <- data.frame("Type"=c("File: " , "Class: " ,"Dimension: ", "Assay: ", "Rownames: ","Colnames: "),"Description" = c(
        paste0(i,".csv"),
        class(Meta[[i]]),
        paste(nr,"x",nc),
        Assay,
        paste(rownames(Meta[[i]])[1],"...",rownames(Meta[[i]])[nrow(Meta[[i]])]),
        paste(colnames(Meta[[i]])[2],"...",colnames(Meta[[i]])[ncol(Meta[[i]])])
      ))



      dt <- rbind(dt, ltest,sp)


    }

    if(str_detect(i, c("pheno"))==T){
      nr <- nrow(Meta[[i]])
      nc <- ncol(Meta[[i]])


      if(str_detect(i, "clinic_")){ Assay =  "Original clinical data"}
      if(str_detect(i, "Patient")){ Assay =  "Patient's clinical data"}
      if(str_detect(i, "Sample")){ Assay = "Samples pathological records"}

      ltest <- data.frame("Type"=c("File: " , "Class: " ,"Dimension: ", "Assay: ", "Rownames: ","Colnames: "),"Description" = c(
        paste0(i,".csv"),
        class(Meta[[i]]),
        paste(nr,"x",nc),
        Assay,
        paste(rownames(Meta[[i]])[1],"...",rownames(Meta[[i]])[nrow(Meta[[i]])]),
        paste(colnames(Meta[[i]])[2],"...",colnames(Meta[[i]])[ncol(Meta[[i]])])
      ))



      dt <- rbind(dt, ltest,sp)

    }

  }
  write.table(dt,"Readme.txt",row.names = F, quote = FALSE)
  file.show("Readme.txt")

  closeAllConnections()
} else {


  tme <- Sys.Date()
  tme <- format(tme, format="%B %d %Y")

  name <- paste0(data.norm,"-", project,"-matrix")

  sp <- data.frame("Type"="---------:" ,"Description"="---------")
  mod <- data.frame("Type"=paste(name, "added the: ") ,"Description"=tme)

  dt <- rbind(sp,mod,sp)


  if(str_detect(name, c("matrix"))==T){
    nr <- nrow(Meta[[name]])
    nc <- ncol(Meta[[name]])



    if(str_detect(name, "Raw")){ Assay = "Raw counts"}
    if(str_detect(name, "TPM")){ Assay="TPM normalization"}
    if(str_detect(name, "FKPM")){ Assay="FKPM normalization"}
    if(str_detect(name, "Normalized")){ Assay="Normalized gene expression"}

    ltest <- data.frame("Type"=c("File: " , "Class: " ,"Dimension: ", "Assay: ", "Rownames: ","Colnames: "),"Description" = c(
      paste0(name,".csv"),
      class(Meta[[name]]),
      paste(nr,"x",nc),
      Assay,
      paste(rownames(Meta[[name]])[1],"...",rownames(Meta[[name]])[nrow(Meta[[name]])]),
      paste(colnames(Meta[[name]])[2],"...",colnames(Meta[[name]])[ncol(Meta[[name]])])
    ))


    dt <- rbind(dt, ltest,sp)

  }



  if(str_detect(name, c("clinic"))==T){
    nr <- nrow(Meta[[name]])
    nc <- ncol(Meta[[name]])



    if(str_detect(name, "clinic_")){ Assay =  "Original clinical data"}
    if(str_detect(name, "Patient")){ Assay =  "Patient's clinical data"}
    if(str_detect(name, "Sample")){ Assay = "Samples pathological records"}

    ltest <- data.frame("Type"=c("File: " , "Class: " ,"Dimension: ", "Assay: ", "Rownames: ","Colnames: "),"Description" = c(
      paste0(name,".csv"),
      class(Meta[[name]]),
      paste(nr,"x",nc),
      Assay,
      paste(rownames(Meta[[name]])[1],"...",rownames(Meta[[name]])[nrow(Meta[[name]])]),
      paste(colnames(Meta[[name]])[2],"...",colnames(Meta[[name]])[ncol(Meta[[name]])])
    ))



    dt <- rbind(dt, ltest,sp)


  }


  write.table(dt,"Readme.txt",row.names = F,quote = FALSE, append = T, col.names=FALSE)
  file.show("Readme.txt")
  }


  closeAllConnections()



return(Meta)
}




