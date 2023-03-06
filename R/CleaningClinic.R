#' CleaningClinic function : Cleaning clinical data to map samples and patients in Sequençing data
#'
#' @param Metadata Metadata object
#' @param type "Sample Pheno" or "Patients' clinical data" for building clean clinical data from raw clinical data.
#' @param Lexical_colnames_path file path to find lexique of colnames
#' @importFrom utils menu
#' @import dplyr
#' @return a data frame of Samples pheno or patients clinical data. If Sample ID and Patients ID are the sames, so Samples.pheno and Patient_clinic are the same data frame
#' @export
#'
#' @examples "none"
CleaningClinic <- function(Metadata, type = c("Sample", "Patients"), Lexical_colnames_path){


if(!exists("LexicClinic", mode = "any")){
  pos <- 1
  envir = as.environment(pos)
  assign("LexicClinic", ColNameClinic(Lexical_colnames_path), envir = envir)
}

   NB <- which(attributes(Metadata)$Data.Type=="Clinical.data" & attributes(Metadata)$Raw.data=="Yes")

  if(length(NB)==0){stop("No clinical data in meta object")}



  clinic <- Metadata[[NB]]

  clcl <-  data.frame(matrix(nrow = nrow(clinic), ncol = length(LexicClinic)))
  colnames(clcl) = names(LexicClinic)
  LexicClinic <- lapply(LexicClinic, toupper)




for (i in 1:ncol(clinic)) {
  pat <- toupper(colnames(clinic)[i])
  col <- grep(paste("\\b",pat,"\\b", sep=""), LexicClinic)

  if(!length(col)==0){
    print(colnames(clinic[,i]))
    clcl[,col] <- clinic[,i]

  }
}

  if(!exists("SamplesOrPatients", mode = "any")){

  SamplesOrPatients <- data.table::fread(paste(Lexical_colnames_path,"SamplesOrPatients.txt",sep = "/"))
  pos <- 1
  envir = as.environment(pos)
  assign("SamplesOrPatients", SamplesOrPatients, envir = envir)


  }

  if(type=="Sample"){

    cc <-  SamplesOrPatients[which(SamplesOrPatients$Type=="Samples.Pheno" | SamplesOrPatients$Type=="Both"),]$Descirption
    cc <- cc[!cc%in%("SamplesID")]

    if(all(is.na(clcl$SamplesID))){
      message("No SamplesID found in raw clinical data. Using PatientID instead")
      clcl$SamplesID <- clcl$PatientsID
      }

    if(length(which(duplicated(clcl$SamplesID)))==0) {

       clinic2 <-  clcl[,c("SamplesID", cc)]

      clinic2[clinic2==""] <- NA
      clinic2[clinic2=="NA"] <- NA

      Metadata$Sample.pheno <- clinic2

      if(length(attributes(Metadata)$Data.Type)<length(Metadata)){
      attributes(Metadata)$Data.Type <-  c(attributes(Metadata)$Data.Type,"Clinical.data")
      attributes(Metadata)$Raw.data <- c(attributes(Metadata)$Raw.data,"No")}


    } else {

      cl_rolled <- clcl %>%

      # create groups by name
      group_by(SamplesID) %>%

        dplyr::summarise(across(everything(), ~paste0(unique(na.omit(.x)), collapse = ";")))
        cl_rolled <- as.data.frame(cl_rolled)


      cl_rolled <-  cl_rolled[,c("SamplesID", cc)]

      cl_rolled[cl_rolled==""] <- NA
      cl_rolled[cl_rolled=="NA"] <- NA
      Metadata$Sample.pheno <- cl_rolled

      if(length(attributes(Metadata)$Data.Type)<length(Metadata)){
      attributes(Metadata)$Data.Type <-  c(attributes(Metadata)$Data.Type,"Clinical.data")
      attributes(Metadata)$Raw.data <- c(attributes(Metadata)$Raw.data,"No")}
      }



    } else if(type=="Patients")
      {

      cc <-  SamplesOrPatients[which(SamplesOrPatients$Type=="Patients.Pheno" | SamplesOrPatients$Type=="Both"),]$Descirption
      cc <- cc[!cc%in%("PatientsID")]

      if(all(is.na(clcl$PatientsID))){
        message("No PatientsID found in raw clinical data. Using PatientID instead")
        clcl$PatientsID <- clcl$SamplesID
      }

      if(length(which(duplicated(clcl$PatientsID)))==0) {

        if(all(is.na(clcl$PatientsID))){
          message("No PatientsID found in raw clinical data. Using PatientID instead")
          clinic2 <-  clcl[,c("SamplesID", cc)]
          } else {
            clinic2 <-  clcl[,c("PatientsID", cc)]
           }

        clinic2 <-  clcl[,c("PatientsID", cc)]
        clinic2[clinic2==""] <- NA
        clinic2[clinic2=="NA"] <- NA
        Metadata$Patient.clinic <- clinic2

        if(length(attributes(Metadata)$Data.Type)<length(Metadata)){
        attributes(Metadata)$Data.Type <-  c(attributes(Metadata)$Data.Type,"Clinical.data")
        attributes(Metadata)$Raw.data <- c(attributes(Metadata)$Raw.data,"No")}

        } else {




          cl_rolled <- clcl %>%

            # create groups by name
            group_by(PatientsID) %>%

            dplyr::summarise(across(everything(), ~paste0(unique(na.omit(.x)), collapse = ";")))

          isNA <- which(is.na( cl_rolled$PatientsID))

          if(length(isNA)>0){ cl_rolled <- as.data.frame(cl_rolled[-isNA,])  }

          cl_rolled <- as.data.frame(cl_rolled)



          cl_rolled <-  cl_rolled[,c("PatientsID", cc)]
          cl_rolled[cl_rolled==""] <- NA
          cl_rolled[cl_rolled=="NA"] <- NA
          Metadata$Patient.clinic <- cl_rolled

          if(length(attributes(Metadata)$Data.Type)<length(Metadata)){
          attributes(Metadata)$Data.Type <-  c(attributes(Metadata)$Data.Type,"Clinical.data")
          attributes(Metadata)$Raw.data <- c(attributes(Metadata)$Raw.data,"No")}






          }


      } else {
   stop("Choose type = 'Sample Pheno' or 'Patients' clinical data' ")}



  return(Metadata)

}

