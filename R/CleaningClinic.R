#' CleaningClinic function : Cleaning clinical data to map samples and patients in Sequençing data
#'
#' @param Metadata Metadata object
#' @param ClinicToClean a character string of name of clinical data to clean.
#' @param exportname name to apply in Metadata object list
#' @param type "c("Samples", "Patients", "Cells) for building clean clinical data from raw clinical data.
#' @param Lexic Lexic to use for cleaning. Created from CreateLexic function .
#' @param FilterSamples default F, if T, keep only retrieved samples in SamplesAnnotation file
#' @param FilterPatients default F, if T, keep only retrieved patients in Count SamplesAnnotation file
#' @param CleanFromOtherType If TRUE, Force a cleaning method from Samples.clinical data into Patient.clinical data and vice versa.
#' @param force.replace set as F. T : replace an already object with the same name
#' @param keep.all.column default T, copy all column from clinic, else only column from Lexic are built.
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
                           type = c("Samples", "Patients", "Cells"),
                           FilterSamples= F,
                           FilterPatients = F,
                           CleanFromOtherType = F,
                           force.replace = F,
                           keep.all.column = T){






  if(is.null(ClinicToClean)){stop("ClinicToClean must be a character")}
  if(!inherits(ClinicToClean,what ="character" )){ stop("ClinicToClean must be a character")}
  if(is.null(exportname)){stop("exportname must be a character")}
  if(!inherits(exportname,what ="character" )){ stop("exportname must be a character")}
  if(!ClinicToClean%in%names(Metadata)) {stop(paste0(ClinicToClean,"is not in Metaobject")) }
  if(is.null(Lexic)){stop("A Lexic must be referred to. Create one with CreateLexic() function.")}



  if(exportname%in%names(Metadata)){
    message("An Object with the same exportname already exist in MetaObject")
    if(force.replace==F){stop("set force.replace==T to subset object.")}}



    if(type=="Samples"){
      message("------------------------")
      message("Samples annotation cleaning")

    NBS <- which(attributes(Metadata)$Data.Type=="SamplesAnnot" & attributes(Metadata)$Cleaned == "No")


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
    colnames(clinic) <-gsub(" ", "_",colnames(clinic))
    colnames(clinic) <-gsub("[:]", "_",colnames(clinic))
    colnames(clinic) <-gsub("[(]", "",colnames(clinic))
    colnames(clinic) <-gsub("[)]", "",colnames(clinic))
    colnames(clinic) <-gsub("[,]", "",colnames(clinic))
    colnames(clinic) = gsub("[[:punct:]]","-", colnames(clinic))
    LexicClinic <-  lapply(LexicClinic,function(x) gsub("[.]", "_",x))
    LexicClinic <-  lapply(LexicClinic,function(x) gsub("[:]", "_",x))
    LexicClinic <-  lapply(LexicClinic,function(x) gsub(" ", "_",x))
    LexicClinic <-  lapply(LexicClinic,function(x) gsub("[(]", "",x))
    LexicClinic <-  lapply(LexicClinic,function(x) gsub("[)]", "",x))
    LexicClinic <-  lapply(LexicClinic,function(x) gsub("[,]", "",x))
    LexicClinic = lapply(LexicClinic,function(x) gsub("[[:punct:]]","-",x))

    if(file.exists(paste0(Processpath(Metadata),"/Samples.CleanedProcess.txt"))){ file.remove(paste0(Processpath(Metadata),"/Samples.CleanedProcess.txt"))}

    cat("Samples.CleanedProcess" , file=paste0(Processpath(Metadata),"/Samples.CleanedProcess.txt"),sep="\n",append = T)
    cat("Raw.Clinic colnames Origin,Clean.Called" , file=paste0(Processpath(Metadata),"/Samples.CleanedProcess.txt"),sep="\n",append = T)





 for (i in 1:ncol(clinic)) {
      pat <- toupper(colnames(clinic)[i])
      col <- grep(paste("\\b",pat, "\\b",sep=""), LexicClinic)


    #      if(length(col)==0){ col <- grep(pat, LexicClinic)}

      cat(paste(pat,",", names(LexicClinic)[col]) , file=paste0(Processpath(Metadata),"/Samples.CleanedProcess.txt"),sep="\n",append = T)

      if(!length(col)==0){


        substrings <-   na.omit(clcl[,col])
        retrivedcols = na.omit(colnames(clcl)[apply(clcl,1, function(x)  x %in% substrings  )])

        if(all(is.na(clcl[,col] ))){
          clcl[,col] <- clinic[,i]} else {message(paste(names(LexicClinic)[col],
                                                        " already entered previously.\nCheck logs ('SampleLog()') and modify Lexic to match your needs.\n"))}

        if(length(retrivedcols)>1){
          message(paste(colnames(clcl)[col], "values are duplicated in a previous column:" ,
                        retrivedcols[retrivedcols!=colnames(clcl)[col]],"\n"))
        }

      }
    }


    cc <- names(LexicClinic)
    cc <- cc[!cc%in%("samplesID")]

    if(all(is.na(clcl$samplesID))){
      message("No SamplesID found in raw clinical data. Using PatientID instead")
      clcl$samplesID <- clcl$patientsID
    }

    if(length(which(duplicated(clcl$samplesID)))==0) {

      clinic2 <-  clcl[,c("samplesID", cc)]

      clinic2 <- as.data.frame(clinic2)
      clinic2[clinic2==""] <- NA
      clinic2[clinic2=="NA"] <- NA



      if(FilterPatients==T){stop("Data type is not 'Samples', you can't substet it by SamplesID. Try FilterPatients = T.")}
      if(FilterSamples==T){

        if(attributes(Metadata)$Omics.type != "Single.Cell"){
        message("Selecting only Samples present in both Count and SamplesAnnot.")


          zz = which(attributes(Metadata)$Data.Type=="Count")[1]

        # Diviser chaque élément de la colonne en un vecteur de sous-chaînes
        substrings <- clinic2$samplesID
        # Vérifier si chaque valeur du vecteur est présente dans chaque vecteur de sous-chaînes
        est_present <- sapply(substrings, function(x) any(colnames(Metadata[[zz]]) %in% x))

        if(all(est_present==F)){
          substrings <- clinic2$patientsID
          est_present <- sapply(substrings, function(x) any(colnames(Metadata[[zz]]) %in% x))}

        # Récupérer les index de position correspondants
        indices <- which(est_present)
        clinic2 <- clinic2[indices,]
        } else {

          zz = which(attributes(Metadata)$Data.Type=="Count")[1]
          substrings <- clinic2$samplesID
          names(substrings) = substrings
          for(z in substrings){
            if("TRUE" %in% str_detect(pattern = paste0(z,"-"), colnames(Metadata[[zz]]))){substrings[z] = T }else { substrings[z] = F}
          }
          clinic2 <- clinic2[which(substrings==T),]


        }}

      rownames(clinic2) = clinic2$samplesID

      if(exportname%in%names(Metadata)){
        if(force.replace==F){stop("set force.replace==T to subset object.")}
        message("Subsetting object.")

        Metadata[[exportname]] <- clinic2

        t= which(str_detect(names(Metadata),exportname))

        attributes(Metadata)$Data.Type[t] <-  c("SamplesAnnot")
        attributes(Metadata)$Export[t] <- c("Yes")
        attributes(Metadata)$Cleaned[t] = c("Yes")



      }else {

        l <- length(Metadata)
        Metadata[[l+1]] <- clinic2
        names(Metadata)[l+1] <- exportname
        attributes(Metadata)$Data.Type <-  c(attributes(Metadata)$Data.Type,"SamplesAnnot")
        attributes(Metadata)$Export <- c(attributes(Metadata)$Export,"Yes")
        attributes(Metadata)$Cleaned = c(attributes(Metadata)$Cleaned,"Yes")}





    } else {

      cl_rolled <- clcl %>%

        # create groups by name
        group_by(samplesID) %>%

        dplyr::summarise(across(everything(), ~paste0(unique(na.omit(.x)), collapse = ";")))
      cl_rolled <- as.data.frame(cl_rolled)


      cl_rolled <-  cl_rolled[,c("samplesID", cc)]

      cl_rolled[cl_rolled==""] <- NA
      cl_rolled[cl_rolled=="NA"] <- NA


      if(exportname%in%names(Metadata)){
        if(force.replace==F){stop("set force.replace==T to subset object.")}
        message("Subsetting object.")


        if(FilterSamples==T){
          message("Selecting only Samples present in both Count and SamplesAnnot.")


          zz = which(attributes(Metadata)$Data.Type=="Count")[1]

          # Diviser chaque élément de la colonne en un vecteur de sous-chaînes
          substrings <- cl_rolled$samplesID
          # Vérifier si chaque valeur du vecteur est présente dans chaque vecteur de sous-chaînes
          est_present <- sapply(substrings, function(x) any(colnames(Metadata[[zz]]) %in% x))

          if(all(est_present==F)){
            substrings <- cl_rolled$patientsID
            est_present <- sapply(substrings, function(x) any(colnames(Metadata[[zz]]) %in% x))}

          # Récupérer les index de position correspondants
          indices <- which(est_present)
          cl_rolled <- cl_rolled[indices,]
        }
        rownames(cl_rolled) = cl_rolled$samplesID

        Metadata[[exportname]] <- cl_rolled

        t= which(str_detect(names(Metadata),exportname))

        attributes(Metadata)$Data.Type[t] <-  c("SamplesAnnot")
        attributes(Metadata)$Export[t] <- c("Yes")
        attributes(Metadata)$Cleaned[t] = c("Yes")


      }else {

        l <- length(Metadata)
        Metadata[[l+1]] <- cl_rolled
        names(Metadata)[l+1] <- exportname
        attributes(Metadata)$Data.Type <-  c(attributes(Metadata)$Data.Type,"SamplesAnnot")
        attributes(Metadata)$Export <- c(attributes(Metadata)$Export,"Yes")
        attributes(Metadata)$Cleaned = c(attributes(Metadata)$Cleaned,"Yes")}
    }

  #  file.show(paste0(Processpath(Metadata),"/Samples.CleanedProcess.txt"))

  } else if(type=="Patients")
  {
    message("------------------------")
    message("Patients clinic cleaning")


    NBP <- which(attributes(Metadata)$Data.Type=="Clinic" & attributes(Metadata)$Cleaned == "No")

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
    colnames(clinic) <-gsub(" ", "_",colnames(clinic))
    colnames(clinic) <-gsub("[:]", "_",colnames(clinic))
    colnames(clinic) <-gsub("[(]", "",colnames(clinic))
    colnames(clinic) <-gsub("[)]", "",colnames(clinic))
    colnames(clinic) <-gsub("[,]", "",colnames(clinic))
    colnames(clinic) = gsub("[[:punct:]]","-", colnames(clinic))
    LexicClinic <-  lapply(LexicClinic,function(x) gsub("[.]", "_",x))
    LexicClinic <-  lapply(LexicClinic,function(x) gsub("[:]", "_",x))
    LexicClinic <-  lapply(LexicClinic,function(x) gsub(" ", "_",x))
    LexicClinic <-  lapply(LexicClinic,function(x) gsub("[(]", "",x))
    LexicClinic <-  lapply(LexicClinic,function(x) gsub("[)]", "",x))
    LexicClinic <-  lapply(LexicClinic,function(x) gsub("[,]", "",x))
    LexicClinic = lapply(LexicClinic,function(x) gsub("[[:punct:]]","-",x))

    if(file.exists(paste0(Processpath(Metadata),"/Patients.CleanedProcess.txt"))){ file.remove(paste0(Processpath(Metadata),"/Patients.CleanedProcess.txt"))}

    cat("Patients.CleanedProcess" , file=paste0(Processpath(Metadata),"/Patients.CleanedProcess.txt"),sep="\n",append = T)
    cat("Raw.Clinic colnames Origin,Clean.Called" , file=paste0(Processpath(Metadata),"/Patients.CleanedProcess.txt"),sep="\n",append = T)


    for (i in 1:ncol(clinic)) {
      pat <- toupper(colnames(clinic)[i])
      col <- grep(paste("\\b",pat, "\\b",sep=""), LexicClinic)


    #      if(length(col)==0){ col <- grep(pat, LexicClinic)}

      cat(paste(pat,",", names(LexicClinic)[col]) , file=paste0(Processpath(Metadata),"/Patients.CleanedProcess.txt"),sep="\n",append = T)

      if(!length(col)==0){



        substrings <-   na.omit(clcl[,col])
        retrivedcols = na.omit(colnames(clcl)[apply(clcl,1, function(x)  x %in% substrings  )])

        if(all(is.na(clcl[,col] ))){
          clcl[,col] <- clinic[,i]} else {message(paste(names(LexicClinic)[col],
                                                        " already entered previously.\nCheck logs ('PatientLog()') and modify Lexic to match your need.\n"))}

        if(length(retrivedcols)>1){
          message(paste(colnames(clcl)[col], "values are duplicated in a previous column:" ,
                        retrivedcols[retrivedcols!=colnames(clcl)[col]],"\n"))
        }

      }
    }

    cc <- names(LexicClinic)
    cc <- cc[!cc%in%("patientsID")]

    if(all(is.na(clcl$patientsID))){
      message("No PatientsID found in raw clinical data. Using SamplesID instead")
      clcl$patientsID <- clcl$samplesID
    }





    if(length(which(duplicated(clcl$patientsID)))==0) {


      if(!"samplesID" %in%colnames(clcl) ){
        clcl$samplesID = NA
      }

      clinic2 <-  clcl[,c("patientsID", cc)]
      clinic2 <- as.data.frame(clinic2)
      clinic2[clinic2==""] <- NA
      clinic2[clinic2=="NA"] <- NA


      NBS <- which(attributes(Metadata)$Data.Type=="SamplesAnnot" & attributes(Metadata)$Cleaned == "No")

      if(length(NBS)>0){

        if(all(is.na(clinic2$samplesID))){} else {  psID =  clinic2$samplesID[which(clinic2$samplesID %in% Metadata[[NBS]]$samplesID)]




        SID = Metadata[[NBS]]$samplesID
        PID = clinic2$patientsID[which(clinic2$patientsID %in% Metadata[[NBS]][,"patientsID"])]

        if(length(psID) <= length(SID)){
          if(length(psID)>length(PID)) {
            message("More concordance by SamplesID than PatientsID")
            message("Checking Clinical data and Samples Annotation concordance by SamplesID")
            suma = summary(clinic2$samplesID %in% Metadata[[NBS]]$samplesID)
            print(suma)

            message(paste("Keeping all patientsID in", exportname))

          } else {

            message("Checking Clinical data and Samples Annotation concordance by PatientsID.")
            NBS = NBS[1]
            suma = summary(clinic2$patientsID %in% Metadata[[NBS]][,"patientsID"])
            print(suma)

            message(paste("Keeping all patientsID in", exportname))




            }

        } else {

          message("Checking Clinical data and Samples Annotation concordance by PatientsID.")
          NBS = NBS[1]
          suma = summary(clinic2$patientsID %in% Metadata[[NBS]][,"patientsID"])
          print(suma)

          message(paste("Keeping all patientsID in", exportname))

          }
      }}

      rownames(clinic2) = clinic2$patientsID

      if(FilterSamples ==T){stop("Data type is not 'Patients',you can't substet it by PatientsID. Try FilterSamples = T.")}
      if(FilterPatients==T){

      if(attributes(Metadata)$Omics.type != "Single.Cell"){

        if(!all(is.na(clinic2$samplesID))){
        message("Selecting only Patient present in both Count and clinical data.")

        zz = which(attributes(Metadata)$Data.Type=="Count")[1]

        # Diviser chaque élément de la colonne en un vecteur de sous-chaînes
        substrings <- clinic2$patientsID
        # Vérifier si chaque valeur du vecteur est présente dans chaque vecteur de sous-chaînes
        est_present <- sapply(substrings, function(x) any(colnames(Metadata[[zz]]) %in% x))

        if(all(est_present==F)){
        substrings <- clinic2$samplesID

        substrings <- clinic2$samplesID
        MultiSamples = all(na.omit(str_detect(substrings, ";"))==F)

        if(MultiSamples==F){

          substrings = na.omit(unlist(str_split(substrings, ";")))

          NBS <- which(attributes(Metadata)$Data.Type=="SamplesAnnot" & attributes(Metadata)$Cleaned == "No")
          NBSc <- which(attributes(Metadata)$Data.Type=="SamplesAnnot" & attributes(Metadata)$Cleaned == "Yes")

          if(length(NBSc)>0){
          filter.pat= na.omit(Metadata[[NBSc[1]]][substrings,"patientsID"])
          }
          if(length(NBS)>0 & length(NBSc)==0){
            filter.pat=  na.omit(Metadata[[NBS[1]]][substrings,"patientsID"])
           }

          } else {

        est_present <- sapply(substrings, function(x) any(colnames(Metadata[[zz]]) %in% x))}





        }

        # Récupérer les index de position correspondants
        if(exists("filter.pat")){indices=unique(filter.pat)
         }else{
        indices <- which(est_present)}
        clinic2 <- clinic2[indices,]
     } else { message("No filtering doable.")}

      }else {

        zz = which(attributes(Metadata)$Data.Type=="Count")[1]

        substrings <- clinic2$patientsID
        names(substrings) = substrings
        for(z in substrings){
        if("TRUE" %in% str_detect(pattern = paste0(z,"-"), colnames(Metadata[[zz]]))){substrings[z] = T }else { substrings[z] = F}
          }
        clinic2 <- clinic2[which(substrings==T),]

      }



        }





      if(exportname%in%names(Metadata)){
        if(force.replace==F){stop("set force.replace==T to subset object.")}
        message("Subsetting object.")

        Metadata[[exportname]] <- clinic2

        t= which(str_detect(names(Metadata),exportname))

        attributes(Metadata)$Data.Type[t] <-  c("Clinic")
        attributes(Metadata)$Export[t] <- c("Yes")
        attributes(Metadata)$Cleaned[t] = c("Yes")


      }else {

        l <- length(Metadata)
        Metadata[[l+1]] <- clinic2
        names(Metadata)[l+1] <- exportname
        attributes(Metadata)$Data.Type <-  c(attributes(Metadata)$Data.Type,"Clinic")
        attributes(Metadata)$Export <- c(attributes(Metadata)$Export,"Yes")
        attributes(Metadata)$Cleaned = c(attributes(Metadata)$Cleaned,"Yes")}



    } else {




      cl_rolled <- clcl %>%

        # create groups by name
        group_by(patientsID) %>%

        dplyr::summarise(across(everything(), ~paste0(unique(na.omit(.x)), collapse = ";")))

      isNA <- which(is.na( cl_rolled$patientsID))

      if(length(isNA)>0){ cl_rolled <- as.data.frame(cl_rolled[-isNA,])  }

      cl_rolled <- as.data.frame(cl_rolled)



      cl_rolled <-  cl_rolled[,c("patientsID", cc)]
      cl_rolled[cl_rolled==""] <- NA
      cl_rolled[cl_rolled=="NA"] <- NA

      NBS <- which(attributes(Metadata)$Data.Type=="SamplesAnnot" & attributes(Metadata)$Cleaned == "No")

      if(length(NBS)>0){


        psID =  cl_rolled$samplesID[which(cl_rolled$samplesID %in% Metadata[[NBS]]$samplesID)]
        SID = Metadata[[NBS]]$samplesID
        PID = cl_rolled$patientsID[which(cl_rolled$patientsID %in% Metadata[[NBS]][,"patientsID"])]

        if(length(psID) <= length(SID)){
          if(length(psID)>length(PID)) {
            message("More concordance by SamplesID than PatientsID")
            message("Checking Clinical data and Samples Annotation concordance by SamplesID")
            suma = summary(cl_rolled$samplesID %in% Metadata[[NBS]]$samplesID)
            print(suma)


            message(paste("Keeping all samplesID in", exportname))

          } else {

            message("Checking Clinical data and Samples Annotation concordance by PatientsID.")
            NBS = NBS[1]
            suma = summary(cl_rolled$patientsID %in% Metadata[[NBS]][,"patientsID"])
            print(suma)

            message(paste("Keeping all patientsID in", exportname))



          }

        } else {

          message("Checking Clinical data and Samples Annotation concordance by PatientsID.")
          NBS = NBS[1]
          suma = summary(cl_rolled$patientsID %in% Metadata[[NBS]][,"patientsID"])
          print(suma)

          message(paste("Keeping all patientsID in", exportname))


        }

      }

      rownames(cl_rolled) = cl_rolled$patientsID

      if(FilterPatients==T){
        message("Selecting only Patient present in both Count and clinical data.")
        if(attributes(Metadata)$Omics.type != "Single.Cell"){
        zz = which(attributes(Metadata)$Data.Type=="Count")[1]

        # Diviser chaque élément de la colonne en un vecteur de sous-chaînes
        substrings <- cl_rolled$patientsID
        # Vérifier si chaque valeur du vecteur est présente dans chaque vecteur de sous-chaînes
        est_present <- sapply(substrings, function(x) any(colnames(Metadata[[zz]]) %in% x))

        if(all(est_present==F)){
          substrings <- cl_rolled$samplesID
          est_present <- sapply(substrings, function(x) any(colnames(Metadata[[zz]]) %in% x))

          MultiSamples = all(na.omit(str_detect(substrings, ";"))==F)

          if(MultiSamples==F){

            substrings = na.omit(unlist(str_split(substrings, ";")))

            NBS <- which(attributes(Metadata)$Data.Type=="SamplesAnnot" & attributes(Metadata)$Cleaned == "No")
            NBSc <- which(attributes(Metadata)$Data.Type=="SamplesAnnot" & attributes(Metadata)$Cleaned == "Yes")

            if(length(NBSc)>0){
              filter.pat= na.omit(Metadata[[NBSc[1]]][substrings,"patientsID"])}
            if(length(NBS)>0 & length(NBSc)==0){
              filter.pat=  na.omit(Metadata[[NBS[1]]][substrings,"patientsID"])}

          } else {

            est_present <- sapply(substrings, function(x) any(colnames(Metadata[[zz]]) %in% x))}





        }

        # Récupérer les index de position correspondants
        if(exists("filter.pat")){indices=unique(filter.pat)}else{
          indices <- which(est_present)}
        cl_rolled <- cl_rolled[indices,]
          }

        else {

          zz = which(attributes(Metadata)$Data.Type=="Count")[1]

          substrings <- cl_rolled$patientsID
          names(substrings) = substrings
          for(z in substrings){
            if("TRUE" %in% str_detect(pattern = paste0(z,"-"), colnames(Metadata[[zz]]))){substrings[z] = T }else { substrings[z] = F}
          }
          cl_rolled <- cl_rolled[which(substrings==T),]

        }





      }




      if(exportname%in%names(Metadata)){
        if(force.replace==F){stop("set force.replace==T to subset object.")}
        message("Subsetting object.")

        Metadata[[exportname]] <- cl_rolled

        t= which(str_detect(names(Metadata),exportname))

        attributes(Metadata)$Data.Type[t] <-  c("Clinic")
        attributes(Metadata)$Export[t] <- c("Yes")
        attributes(Metadata)$Cleaned[t] = c("Yes")


      }else {

        l <- length(Metadata)
        Metadata[[l+1]] <- cl_rolled
        names(Metadata)[l+1] <- exportname
        attributes(Metadata)$Data.Type <-  c(attributes(Metadata)$Data.Type,"Clinic")
        attributes(Metadata)$Export <- c(attributes(Metadata)$Export,"Yes")
        attributes(Metadata)$Cleaned = c(attributes(Metadata)$Cleaned,"Yes")}





    }

    #file.show(paste0(Processpath(Metadata),"/Patients.CleanedProcess.txt"))


  } else if(type == "Cells"){

    message("------------------------")
    message("Cells annotation cleaning")
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
    colnames(clinic) <-gsub(" ", "_",colnames(clinic))
    colnames(clinic) <-gsub("[:]", "_",colnames(clinic))
    colnames(clinic) <-gsub("[(]", "",colnames(clinic))
    colnames(clinic) <-gsub("[)]", "",colnames(clinic))
    colnames(clinic) <-gsub("[,]", "",colnames(clinic))
    LexicClinic <-  lapply(LexicClinic,function(x) gsub("[.]", "_",x))
    LexicClinic <-  lapply(LexicClinic,function(x) gsub("[:]", "_",x))
    LexicClinic <-  lapply(LexicClinic,function(x) gsub(" ", "_",x))
    LexicClinic <-  lapply(LexicClinic,function(x) gsub("[(]", "",x))
    LexicClinic <-  lapply(LexicClinic,function(x) gsub("[)]", "",x))
    LexicClinic <-  lapply(LexicClinic,function(x) gsub("[,]", "",x))

    if(file.exists(paste0(Processpath(Metadata),"/CellAnnotation.CleanedProcess.txt"))){ file.remove(paste0(Processpath(Metadata),"/CellAnnotation.CleanedProcess.txt"))}

    cat("CellAnnotation.CleanedProcess" , file=paste0(Processpath(Metadata),"/CellAnnotation.CleanedProcess.txt"),sep="\n",append = T)
    cat("Raw.CellAnnotation colnames Origin,Clean.Called" , file=paste0(Processpath(Metadata),"/CellAnnotation.CleanedProcess.txt"),sep="\n",append = T)

    for (i in 1:ncol(clinic)) {

      pat <- toupper(colnames(clinic)[i])
      col <- grep(paste("\\b",pat, "\\b",sep=""), LexicClinic)


    #      if(length(col)==0){ col <- grep(pat, LexicClinic)}

      cat(paste(pat,",", names(LexicClinic)[col]) , file=paste0(Processpath(Metadata),"/CellAnnotation.CleanedProcess.txt"),sep="\n",append = T)

      if(!length(col)==0){



        substrings <-   na.omit(clcl[,col])
        retrivedcols = na.omit(colnames(clcl)[apply(clcl,1, function(x)  x %in% substrings  )])

        if(all(is.na(clcl[,col] ))){
          clcl[,col] <- clinic[,i]} else {message(paste(names(LexicClinic)[col],
                                                        " already entered previously.\nCheck logs ('CellLog()') and modify Lexic to match your need.\n"))}

        if(length(retrivedcols)>1){
          message(paste(colnames(clcl)[col], "values are duplicated in a previous column:" ,
                        retrivedcols[retrivedcols!=colnames(clcl)[col]],"\n"))
        }

      }

      }


    clcl$CellsBarcode = clinic$CellsBarcode
    cc <- names(LexicClinic)
    cc <- cc[!cc%in%("CellsBarcode")]

    if(all(is.na(clcl$CellsBarcode))){
      stop("No Cells found in raw Cells annotation.")

    }



      clinic2 <-  clcl[,c("CellsBarcode", cc)]
      clinic2 <- as.data.frame(clinic2)
      clinic2[clinic2==""] <- NA
      clinic2[clinic2=="NA"] <- NA


        message("Selecting only Cells present in both Count and Cells annotation data.")

        zz = which(attributes(Metadata)$Data.Type=="Count")[1]

        clinic2 <- clinic2[which(clinic2$CellsBarcode%in%colnames(Metadata[[zz]])),]

        rownames(clinic2) = clinic2$CellsBarcode


      if(exportname%in%names(Metadata)){
        if(force.replace==F){stop("set force.replace==T to subset object.")}
        message("Subsetting object.")

        Metadata[[exportname]] <- clinic2

        t= which(str_detect(names(Metadata),exportname))

        attributes(Metadata)$Data.Type[t] <-  c("CellsAnnot")
        attributes(Metadata)$Export[t] <- c("Yes")
        attributes(Metadata)$Cleaned[t] = c("Yes")


      }else {

        l <- length(Metadata)
        Metadata[[l+1]] <- clinic2
        names(Metadata)[l+1] <- exportname
        attributes(Metadata)$Data.Type <-  c(attributes(Metadata)$Data.Type,"CellsAnnot")
        attributes(Metadata)$Export <- c(attributes(Metadata)$Export,"Yes")
        attributes(Metadata)$Cleaned = c(attributes(Metadata)$Cleaned,"Yes")}





  } else{
    stop("Choose type = c('Samples', 'Patients', 'Cells')")}





  return(Metadata)

}



