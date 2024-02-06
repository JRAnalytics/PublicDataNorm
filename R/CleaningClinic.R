#' CleaningClinic function : Cleaning clinical data to map samples and patients in Sequençing data
#'
#' @param Metadata Metadata object
#' @param ClinicToClean a character string of name of clinical data to clean.
#' @param exportname name to apply in Metadata object list
#' @param type "c("Samples", "Patients") for building clean clinical data from raw clinical data.
#' @param Lexic Lexic to use for cleaning. Created from CreateLexic function .
#' @param FilterSamples default F, if T, keep only retrieved samples in Count matrix
#' @param FilterPatients default F, if T, keep only retrieved patients in Count matrix
#' @param CleanFromOtherType If TRUE, Force a cleaning method from Samples.clinical data into Patient.clinical data and vice versa.
#' @param force.replace set as F. T : replace an already object with the same name
#' @param keep.all.column default F, if T, copy all column from clinic.
#' @importFrom utils menu
#' @import dplyr
#' @return a data frame of Samples pheno or patients clinical data. If Sample ID and Patients ID are the sames, so Samples.pheno and Patient_clinic are the same data frame
#' @export
#'
#' @examples "none"
CleaningClinic <- function(Metadata = NULL,
                           ClinicToClean = NULL,
                           Lexic = NULL,
                           exportname = NULL,
                           type = c("Samples", "Patients"),
                           FilterSamples= F,
                           FilterPatients = F,
                           CleanFromOtherType = F,
                           force.replace = F,
                           keep.all.column = F){


  list.files.path = attributes(Metadata)$File.path



  if(is.null(ClinicToClean)){stop("ClinicToClean must be a character")}
  if(!inherits(ClinicToClean,what ="character" )){ stop("ClinicToClean must be a character")}
  if(is.null(exportname)){stop("exportname must be a character")}
  if(!inherits(exportname,what ="character" )){ stop("exportname must be a character")}
  if(!ClinicToClean%in%names(Metadata)) {stop(paste0(ClinicToClean,"is not in Metaobject")) }
  if(is.null(Lexic)){stop("A Lexic must be referred to. Create one with CreateLexic() function.")}



  if(!all(str_detect(names(Metadata),exportname)==F)){
    message("An Object with the same exportname already exist in MetaObject")
    if(force.replace==F){stop("set force.replace==T to subset object.")}}



    if(type=="Samples"){

    NBS <- which(attributes(Metadata)$Data.Type=="SamplesAnnot")


    if(CleanFromOtherType==F &length(NBS)==0 ){stop("No SamplesAnnot found in Metadata object. Set CleanFromOtherType=T to force cleaning from Clinic")
    }
    #  if(NBS==which(names(Metadata)%in%ClinicToClean)){

    clinic <- as.data.frame(Metadata[[ClinicToClean]])

    if(keep.all.column==T){

      for (i in colnames(clinic)){

        if (!i %in% names(Lexic)){Lexic <- AddKeyLexic(lexic = Lexic, key =i  ,value = i ) }

      }

    LexicClinic=Lexic

    }else {  LexicClinic <- Lexic }


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



      if(FilterPatients==T){stop("Data type is not 'Samples', you can't substet it by SamplesID. Try FilterPatients = T.")}
      if(FilterSamples==T){


        message("Selecting only Samples present in both Count and SamplesAnnot.")


        zz = which(attributes(Metadata)$Data.Type=="Count")[1]

        # Diviser chaque élément de la colonne en un vecteur de sous-chaînes
        substrings <- clinic2$SamplesID
        # Vérifier si chaque valeur du vecteur est présente dans chaque vecteur de sous-chaînes
        est_present <- sapply(substrings, function(x) any(colnames(Metadata[[zz]]) %in% x))

        if(all(est_present==F)){
          substrings <- clinic2$PatientsID
          est_present <- sapply(substrings, function(x) any(colnames(Metadata[[zz]]) %in% x))}

        # Récupérer les index de position correspondants
        indices <- which(est_present)
        clinic2 <- clinic2[indices,]
      }


      if(!all(str_detect(names(Metadata),exportname)==F)){
        if(force.replace==F){stop("set force.replace==T to subset object.")}
        message("Subsetting object.")

        Metadata[[exportname]] <- clinic2

        t= which(str_detect(names(Metadata),exportname))

        attributes(Metadata)$Data.Type[t] <-  c("SamplesAnnot")
        attributes(Metadata)$Export[t] <- c("Yes")




      }else {

        l <- length(Metadata)
        Metadata[[l+1]] <- clinic2
        names(Metadata)[l+1] <- exportname
        attributes(Metadata)$Data.Type <-  c(attributes(Metadata)$Data.Type,"SamplesAnnot")
        attributes(Metadata)$Export <- c(attributes(Metadata)$Export,"Yes")}





    } else {

      cl_rolled <- clcl %>%

        # create groups by name
        group_by(SamplesID) %>%

        dplyr::summarise(across(everything(), ~paste0(unique(na.omit(.x)), collapse = ";")))
      cl_rolled <- as.data.frame(cl_rolled)


      cl_rolled <-  cl_rolled[,c("SamplesID", cc)]

      cl_rolled[cl_rolled==""] <- NA
      cl_rolled[cl_rolled=="NA"] <- NA


      if(!all(str_detect(names(Metadata),exportname)==F)){
        if(force.replace==F){stop("set force.replace==T to subset object.")}
        message("Subsetting object.")


        if(FilterSamples==T){
          message("Selecting only Samples present in both Count and SamplesAnnot.")


          zz = which(attributes(Metadata)$Data.Type=="Count")[1]

          # Diviser chaque élément de la colonne en un vecteur de sous-chaînes
          substrings <- cl_rolled$SamplesID
          # Vérifier si chaque valeur du vecteur est présente dans chaque vecteur de sous-chaînes
          est_present <- sapply(substrings, function(x) any(colnames(Metadata[[zz]]) %in% x))

          if(all(est_present==F)){
            substrings <- cl_rolled$PatientsID
            est_present <- sapply(substrings, function(x) any(colnames(Metadata[[zz]]) %in% x))}

          # Récupérer les index de position correspondants
          indices <- which(est_present)
          cl_rolled <- cl_rolled[indices,]
        }


        Metadata[[exportname]] <- cl_rolled

        t= which(str_detect(names(Metadata),exportname))

        attributes(Metadata)$Data.Type[t] <-  c("SamplesAnnot")
        attributes(Metadata)$Export[t] <- c("Yes")


      }else {

        l <- length(Metadata)
        Metadata[[l+1]] <- cl_rolled
        names(Metadata)[l+1] <- exportname
        attributes(Metadata)$Data.Type <-  c(attributes(Metadata)$Data.Type,"SamplesAnnot")
        attributes(Metadata)$Export <- c(attributes(Metadata)$Export,"Yes")}
    }

  #  file.show(paste0(list.files.path$Project.Processes,"/Samples.CleanedProcess.txt"))

  } else if(type=="Patients")
  {

    NBP <- which(attributes(Metadata)$Data.Type=="Clinic")

    if(CleanFromOtherType==F & length(NBP)==0){stop("No Clinic found in Metadata object. Set CleanFromOtherType=T to force cleaning from SamplesAnnot")}


    clinic <- as.data.frame(Metadata[[ClinicToClean]])


    if(keep.all.column==T){

      for (i in colnames(clinic)){

        if (!i %in% names(Lexic)){Lexic <- AddKeyLexic(lexic = Lexic, key =i  ,value = i) }

      }
      LexicClinic <- Lexic

    }else {  LexicClinic <- Lexic}



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


      if(FilterSamples ==T){stop("Data type is not 'Patients',you can't substet it by PatientsID. Try FilterSamples = T.")}
      if(FilterPatients==T){


        message("Selecting only Patient present in both Count and clinical data.")

        zz = which(attributes(Metadata)$Data.Type=="Count")[1]

        # Diviser chaque élément de la colonne en un vecteur de sous-chaînes
        substrings <- clinic2$PatientsID
        # Vérifier si chaque valeur du vecteur est présente dans chaque vecteur de sous-chaînes
        est_present <- sapply(substrings, function(x) any(colnames(Metadata[[zz]]) %in% x))

        if(all(est_present==F)){
        substrings <- clinic2$SamplesID
        est_present <- sapply(substrings, function(x) any(colnames(Metadata[[zz]]) %in% x))}

        # Récupérer les index de position correspondants
        indices <- which(est_present)
        clinic2 <- clinic2[indices,]
      }





      if(!all(str_detect(names(Metadata),exportname)==F)){
        if(force.replace==F){stop("set force.replace==T to subset object.")}
        message("Subsetting object.")

        Metadata[[exportname]] <- clinic2

        t= which(str_detect(names(Metadata),exportname))

        attributes(Metadata)$Data.Type[t] <-  c("Clinic")
        attributes(Metadata)$Export[t] <- c("Yes")


      }else {

        l <- length(Metadata)
        Metadata[[l+1]] <- clinic2
        names(Metadata)[l+1] <- exportname
        attributes(Metadata)$Data.Type <-  c(attributes(Metadata)$Data.Type,"Clinic")
        attributes(Metadata)$Export <- c(attributes(Metadata)$Export,"Yes")}



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


      if(FilterPatients==T){
        message("Selecting only Patient present in both Count and clinical data.")

        zz = which(attributes(Metadata)$Data.Type=="Count")[1]

        # Diviser chaque élément de la colonne en un vecteur de sous-chaînes
        substrings <- cl_rolled$PatientsID
        # Vérifier si chaque valeur du vecteur est présente dans chaque vecteur de sous-chaînes
        est_present <- sapply(substrings, function(x) any(colnames(Metadata[[zz]]) %in% x))

        if(all(est_present==F)){
          substrings <- cl_rolled$SamplesID
          est_present <- sapply(substrings, function(x) any(colnames(Metadata[[zz]]) %in% x))}

        # Récupérer les index de position correspondants
        indices <- which(est_present)
        cl_rolled <- cl_rolled[indices,]
      }




      if(!all(str_detect(names(Metadata),exportname)==F)){
        if(force.replace==F){stop("set force.replace==T to subset object.")}
        message("Subsetting object.")

        Metadata[[exportname]] <- cl_rolled

        t= which(str_detect(names(Metadata),exportname))

        attributes(Metadata)$Data.Type[t] <-  c("Clinic")
        attributes(Metadata)$Export[t] <- c("Yes")


      }else {

        l <- length(Metadata)
        Metadata[[l+1]] <- cl_rolled
        names(Metadata)[l+1] <- exportname
        attributes(Metadata)$Data.Type <-  c(attributes(Metadata)$Data.Type,"Clinic")
        attributes(Metadata)$Export <- c(attributes(Metadata)$Export,"Yes")}





    }

    #file.show(paste0(list.files.path$Project.Processes,"/Patients.CleanedProcess.txt"))


  } else {
    stop("Choose type = c('Samples', 'Patients')")}





  return(Metadata)

}



