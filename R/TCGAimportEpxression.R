#' TCGAimportEpxression allow to import and build the rnaseq data from the TCGA
#'
#' @param project a character of TCGA project name
#' @param data.category a character of data categorie
#' @param list.files.path list.files.path
#' @param alt.dir aternaticve dir path for Windows user. Length of path can't exceed 250 character, Downaloading may be altered. The alt.dir path will be remove after moving GDCdata to list.file.path
#' @param data.type a character of data type
#' @param experimental.strategy a character of experimental.strategy
#' @param sample.type a character of sample.type
#' @param data.norm a character "Raw"  = raw count, "TPM" = tpm normalisation gene expression, "FKPM" = fkpm normalisation gene expression if it exist
#' @param platform platform
#' @param file.type filetype
#' @param legacy legacy T or F
#' @param Meta a Meta object. If NULL, will create the Meta object. If exist, will add objects from queries.
#' @param workflow.type workflow.type
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
                                 list.files.path= NULL,
                                 alt.dir = NULL,
                                 project =NULL,
                                   data.category = NULL,
                                   data.type = NULL,
                                   experimental.strategy = NULL,
                                   sample.type = NULL,
                                    platform = NULL,
                                    file.type = NULL,
                                    legacy = F,
                                    workflow.type = NULL,
                                   data.norm = c("Raw", "TPM", "FKPM", "FKPM_UQ")){


  if(is.null(list.files.path)) { stop("A list.files.path is mandatory.")}
if(isServeOK()==FALSE){stop("Connection to server GDC failled")}

work.dir= getwd()

message("Querying", paste(project, "Meta"))

if(!exists(x = "query")){

query <- GDCquery(
  project = project ,
  data.category = data.category,
  data.type = data.type,
  legacy = legacy,
  experimental.strategy = experimental.strategy,
  platform = platform,
  file.type = file.type,
  workflow.type =workflow.type,
  sample.type = sample.type
)

pos <- 1
envir = as.environment(pos)
assign("query", query, envir = envir)
}




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


  ses=sessionInfo()

  if(str_detect(ses$running, "Windows")) {

    if(str_detect(ses$running, "Windows")& nchar(files[1]>150) & is.null(alt.dir)){ stop("File path too long. Try  using an alternative directory to put downloaded data.\n(set 'alt.dir')")}

    if(!is.null(alt.dir)){
      setwd(file.path(paste0(alt.dir,"/")))
      message("setwd(file.path(alt.dir))")

      }else { setwd(list.files.path$Project.RawData)

      alt.dir = list.files.path$Project.RawData}

    pb <- txtProgressBar(min = 0,      # Minimum value of the progress bar
                         max = nrow( query$results[[1]]), # Maximum value of the progress bar
                         style = 3,    # Progress bar style (also available style = 1 and style = 2)
                         width = 50,   # Progress bar width. Defaults to getOption("width")
                         char = "=")   # Character used to create the bar


    for(i in 1:nrow( query$results[[1]])) {

      query2 = query
      query2$results[[1]]=query2$results[[1]][i,]


      suppressMessages(capture.output(
        GDCdownload(
          query = query2,
          method = "api",
          files.per.chunk = NULL)))




      # Sets the progress bar to the current state
      setTxtProgressBar(pb, i)
    }

    close(pb)


    query$results[[1]]$Rename = NA
    for (i in c(1:length(query$results[[1]]$file_name))){

      paste0(alt.dir,"/GDCdata")
      files = file.path(paste0(alt.dir,"/GDCdata"),
                        query$results[[1]]$project, source,
                        gsub(" ","_",query$results[[1]]$data_category),
                        gsub(" ","_",query$results[[1]]$data_type),
                        gsub(" ","_",query$results[[1]]$file_id),
                        gsub(" ","_",query$results[[1]]$file_name)
      )[i]

      file2 = paste0(file.path(paste0(alt.dir,"/GDCdata"),query$results[[1]]$project, source,
                               gsub(" ","_",query$results[[1]]$data_category),
                               gsub(" ","_",query$results[[1]]$data_type),
                               gsub(" ","_",query$results[[1]]$file_id))[i],"/",query[[1]][[1]]$sample.submitter_id[i],".tsv")

      file.rename(from=files,to= file2)
      query$results[[1]]$Rename[i] = paste0(query[[1]][[1]]$sample.submitter_id[i],".tsv")

    }


    file.copy(from = paste0("D:","/GDCdata"),to = list.files.path$Project.RawData,
              recursive = T,overwrite = T )






  } else {

    setwd(list.files.path$Project.RawData)
    GDCdownload(
      query = query,
      method = "api",
      files.per.chunk = 10)
    }





}

setwd(list.files.path$Project.RawData)

if(dir.exists(file.path("GDCdata",project, source))){message(paste(project, "data found"))}



if(str_detect(ses$running, "Windows")) {

  files = file.path("GDCdata",
                    query$results[[1]]$project, source,
                    gsub(" ","_",query$results[[1]]$data_category),
                    gsub(" ","_",query$results[[1]]$data_type),
                    gsub(" ","_",query$results[[1]]$file_id),
                    gsub(" ","_",query$results[[1]]$Rename)
  )



}else {

  files <- file.path(
  query$results[[1]]$project, source,
  gsub(" ","_",query$results[[1]]$data_category),
  gsub(" ","_",query$results[[1]]$data_type),
  gsub(" ","_",query$results[[1]]$file_id),
  gsub(" ","_",query$results[[1]]$file_name)
)

  files <- file.path("GDCdata", files) }



cases <- ifelse(grepl("TCGA|TARGET",query$results[[1]]$project %>% unlist()),query$results[[1]]$cases,query$results[[1]]$sample.submitter_id)


message("Building gene expression data")



data <- TCGA.build(
  ask = data.norm,
  legacy = legacy,
  files = files,
  cases = cases,
  genome = ifelse(query$legacy,"hg19","hg38"),
  experimental.strategy = unique(query$results[[1]]$experimental_strategy)

)



rownames(data) <-  unlist(lapply(str_split(rownames(data),"[|]" ), "[[", 1))








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


attributes(Meta)$Data.Type <- c("Expression.Matrix","Samples.Clinical.data")

attributes(Meta)$Raw.data <- c(ifelse(data.norm=="Raw", "Yes", "No"),"Yes")


names(Meta) <- c(paste0(data.norm,"-", project,"-matrix"), paste0("clinic_",project))

} else {

  l <- length(Meta)
  name <- names(Meta)
  Meta[[l+1]] <- data
  names(Meta) <- c(name,paste0(data.norm,"-", project,"-matrix"))

  attributes(Meta)$Data.Type[[l+1]] <- c("Expression.Matrix")
  attributes(Meta)$Raw.data[[l+1]] <- c(ifelse(data.norm=="Raw", "Yes", "No"))

}


pos <- 1
envir = as.environment(pos)
assign("query", query, envir = envir)

setwd(work.dir)

return(Meta)
}




