#' CleaningClinic function : Cleaning clinical data to map samples and patients in Sequençing data
#'
#' @param Metadata Metadata object
#' @param type "Sample Pheno" or "Patients' clinical data" for building clean clinical data from raw clinical data.
#' @param list.files.path file path to find lexique of colnames
#' @param project project
#' @param ForceCleaning If TRUE, Force a cleaning method from Samples.clinical data into Patient.clinical data and vice versa.
#' @importFrom utils menu
#' @import dplyr
#' @return a data frame of Samples pheno or patients clinical data. If Sample ID and Patients ID are the sames, so Samples.pheno and Patient_clinic are the same data frame
#' @export
#'
#' @examples "none"
CleaningClinic <- function(Metadata, type = c("Sample", "Patients"), list.files.path, project, ForceCleaning = F){




   NB <- which(str_detect(attributes(Metadata)$Data.Type ,"Clinical.data") & attributes(Metadata)$Raw.data=="Yes")

  if(length(NB)==0){stop("No clinical data in meta object")}





  if(type=="Sample"){

    NBS <- which(attributes(Metadata)$Data.Type=="Samples.Clinical.data" & attributes(Metadata)$Raw.data=="Yes")

    if(ForceCleaning==F){
    if(length(NBS)==0){stop("No Samples.Clinical.data found in Metadata object. Set ForceCleaning=T to force cleaning from Patient.Clinical.Data")}
    } else {  if(length(NBS)==0){NBS <- which(str_detect(attributes(Metadata)$Data.Type ,"Clinical.data") & attributes(Metadata)$Raw=="Yes")    }}

    for(k in NBS){

    clinic <- as.data.frame(Metadata[[k]])

    LexicClinic <- SamplesLexic

    clcl <-  data.frame(matrix(nrow = nrow(clinic), ncol = length(LexicClinic)))
    colnames(clcl) = names(LexicClinic)
    LexicClinic <- lapply(LexicClinic, toupper)

    colnames(clinic) <-gsub("[.]", "_",colnames(clinic))
    LexicClinic <-  lapply(LexicClinic,function(x) gsub("[.]", "_",x))

    if(file.exists(paste0(list.files.path$Project.Processes,"/Samples.CleanedProcess.txt"))){ file.remove(paste0(list.files.path$Project.Processes,"/Samples.CleanedProcess.txt"))}

    cat("Samples.CleanedProcess" , file=paste0(list.files.path$Project.Processes,"/Samples.CleanedProcess.txt"),sep="\n",append = T)
    cat("Raw.Clinic colnames Origin,Clean.Called" , file=paste0(list.files.path$Project.Processes,"/Samples.CleanedProcess.txt"),sep="\n",append = T)
    for (i in 1:ncol(clinic)) {
      pat <- toupper(colnames(clinic)[i])
      col <- grep(paste("\\b",pat, "\\b",sep=""), LexicClinic)



      cat(paste(pat,",", names(LexicClinic)[col]) , file=paste0(list.files.path$Project.Processes,"/Samples.CleanedProcess.txt"),sep="\n",append = T)


      if(!length(col)==0){

        clcl[,col] <- clinic[,i]

      }
    }

    cc <- names(LexicClinic)
    cc <- cc[!cc%in%("SamplesID")]

    if(all(is.na(clcl$SamplesID))){
      message("No SamplesID found in raw clinical data. Using PatientID instead")
      clcl$SamplesID <- clcl$PatientsID
      }

    if(length(which(duplicated(clcl$SamplesID)))==0) {

       clinic2 <-  clcl[,c("SamplesID", cc)]

       clinic2 <- as.data.frame(clinic2)
      clinic2[clinic2==""] <- NA
      clinic2[clinic2=="NA"] <- NA

      l <- length(Metadata)

      if(paste0("Cleaned.Samples.pheno.",names(Metadata)[k])%in%names(Metadata)){

        Metadata[[paste0("Cleaned.Samples.pheno.",names(Metadata)[k])]] <- clinic2

        } else {

      Metadata[[l+1]] <- clinic2
      names(Metadata)[l+1] <- paste0("Cleaned.Samples.pheno.",names(Metadata)[k])
      if(length(attributes(Metadata)$Data.Type)<length(Metadata)){
        attributes(Metadata)$Data.Type <-  c(attributes(Metadata)$Data.Type,"Samples.Clinical.data")
        attributes(Metadata)$Raw.data <- c(attributes(Metadata)$Raw.data,"No")}}





    } else {

      cl_rolled <- clcl %>%

      # create groups by name
      group_by(SamplesID) %>%

        dplyr::summarise(across(everything(), ~paste0(unique(na.omit(.x)), collapse = ";")))
        cl_rolled <- as.data.frame(cl_rolled)


      cl_rolled <-  cl_rolled[,c("SamplesID", cc)]

      cl_rolled[cl_rolled==""] <- NA
      cl_rolled[cl_rolled=="NA"] <- NA

      l <- length(Metadata)

      if("Samples.pheno"%in%names(Metadata)){

        Metadata[["Samples.pheno"]] <- cl_rolled

      } else {

        Metadata[[l+1]] <- cl_rolled
        names(Metadata)[l+1] <- paste0("Samples.pheno")
        if(length(attributes(Metadata)$Data.Type)<length(Metadata)){
          attributes(Metadata)$Data.Type <-  c(attributes(Metadata)$Data.Type,"Samples.Clinical.data")
          attributes(Metadata)$Raw.data <- c(attributes(Metadata)$Raw.data,"No")}}

      if(length(attributes(Metadata)$Data.Type)<length(Metadata)){
      attributes(Metadata)$Data.Type <-  c(attributes(Metadata)$Data.Type,"Samples.Clinical.data")
      attributes(Metadata)$Raw.data <- c(attributes(Metadata)$Raw.data,"No")}
      }

    file.show(paste0(list.files.path$Project.Processes,"/Samples.CleanedProcess.txt"))
}
    } else if(type=="Patients")
      {

      NBP <- which(attributes(Metadata)$Data.Type=="Patient.Clinical.data" & attributes(Metadata)$Raw.data=="Yes")



      if(ForceCleaning==F){
        if(length(NBP)==0){stop("No Patient.Clinical.data found in Metadata object. Set ForceCleaning=T to force cleaning from Samples.Clinical.Data")}
      } else { if(length(NBP)==0){ NBP <- which(str_detect(attributes(Metadata)$Data.Type ,"Clinical.data") & attributes(Metadata)$Raw=="Yes")    }}


      for(k in NBP){

        clinic <- as.data.frame(Metadata[[k]])

      LexicClinic <- PatientLexic

      clcl <-  data.frame(matrix(nrow = nrow(clinic), ncol = length(LexicClinic)))
      colnames(clcl) = names(LexicClinic)
      LexicClinic <- lapply(LexicClinic, toupper)

      colnames(clinic) <-gsub("[.]", "_",colnames(clinic))
      LexicClinic <-  lapply(LexicClinic,function(x) gsub("[.]", "_",x))

      if(file.exists(paste0(list.files.path$Project.Processes,"/Patients.CleanedProcess.txt"))){ file.remove(paste0(list.files.path$Project.Processes,"/Patients.CleanedProcess.txt"))}

      cat("Patients.CleanedProcess" , file=paste0(list.files.path$Project.Processes,"/Patients.CleanedProcess.txt"),sep="\n",append = T)
      cat("Raw.Clinic colnames Origin,Clean.Called" , file=paste0(list.files.path$Project.Processes,"/Patients.CleanedProcess.txt"),sep="\n",append = T)

      for (i in 1:ncol(clinic)) {
        pat <- toupper(colnames(clinic)[i])
        col <- grep(paste("\\b",pat, "\\b",sep=""), LexicClinic)

        cat(paste(pat,",", names(LexicClinic)[col]) , file=paste0(list.files.path$Project.Processes,"/Patients.CleanedProcess.txt"),sep="\n",append = T)

        if(!length(col)==0){

          clcl[,col] <- clinic[,i]

        }
      }

      cc <- names(LexicClinic)
      cc <- cc[!cc%in%("PatientsID")]

      if(all(is.na(clcl$PatientsID))){
        message("No PatientsID found in raw clinical data. Using PatientID instead")
        clcl$PatientsID <- clcl$SamplesID
      }

      if(length(which(duplicated(clcl$PatientsID)))==0) {


        clinic2 <-  clcl[,c("PatientsID", cc)]
        clinic2 <- as.data.frame(clinic2)
        clinic2[clinic2==""] <- NA
        clinic2[clinic2=="NA"] <- NA



        l <- length(Metadata)
        if("Patient.Clinical.data"%in%names(Metadata)){

          Metadata[[names(Metadata)[k]]] <- clinic2

        } else {

          Metadata[[l+1]] <- clinic2
          names(Metadata)[l+1] <- "Patient.Clinical.data"
          if(length(attributes(Metadata)$Data.Type)<length(Metadata)){
            attributes(Metadata)$Data.Type <-  c(attributes(Metadata)$Data.Type,"Patient.Clinical.data")
            attributes(Metadata)$Raw.data <- c(attributes(Metadata)$Raw.data,"No")}}



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
          l <- length(Metadata)

          if(paste0("Cleaned.Patient.Clinical.data.",names(Metadata)[k])%in%names(Metadata)){

            Metadata[[paste0("Cleaned.Patient.Clinical.data.",names(Metadata)[k])]] <- cl_rolled

          } else {

            Metadata[[l+1]] <- cl_rolled
            names(Metadata)[l+1] <- paste0("Cleaned.Patient.Clinical.data.",names(Metadata)[k])
            if(length(attributes(Metadata)$Data.Type)<length(Metadata)){
              attributes(Metadata)$Data.Type <-  c(attributes(Metadata)$Data.Type,"Patient.Clinical.data")
              attributes(Metadata)$Raw.data <- c(attributes(Metadata)$Raw.data,"No")}}

          if(length(attributes(Metadata)$Data.Type)<length(Metadata)){
          attributes(Metadata)$Data.Type <-  c(attributes(Metadata)$Data.Type,"Patient.Clinical.data")
          attributes(Metadata)$Raw.data <- c(attributes(Metadata)$Raw.data,"No")}






        }

      file.show(paste0(list.files.path$Project.Processes,"/Patients.CleanedProcess.txt"))
      }

      } else {
   stop("Choose type = 'Sample Pheno' or 'Patients' clinical data' ")}





  return(Metadata)

}

