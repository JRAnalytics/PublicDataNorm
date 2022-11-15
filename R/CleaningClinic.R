
#' CleaningClinic function : Cleaning clinical data to map samples and patients in Sequençing data
#'
#' @param Metadata Metadata object
#' @param type "Sample Pheno" or "Patients' clinical data" for building clean clinical data from raw clinical data.
#' @param list.files.path file path to find lexique of colnames
#' @importFrom utils menu
#' @return a data frame of Samples pheno or patients clinical data. If Sample ID and Patients ID are the sames, so Samples.pheno and Patient_clinic are the same data frame
#' @export
#'
#' @examples "none"
CleaningClinic <- function(Metadata, type = c("Sample Pheno", "Patients' clinical data"), list.files.path,
                           Organism , Organ, Disease, Compartment, SampleType, TisseConservation, TissueSampling, PatientSampling, SamplePathologicalState, SequencingRun
                           ){



  cnameClinic <- ColNameClinic(list.files.path)

  if(all(str_detect(names(Metadata),"clinic"))==T){stop("No clinical data in meta object")}

  clinic <- Metadata[which(str_detect(names(Metadata),"clinic"))][[1]]


  clcl <-  data.frame(matrix(nrow = nrow(clinic), ncol = length(cnameClinic)))
  colnames(clcl) = names(cnameClinic)
  cnameClinic <- lapply(cnameClinic, toupper)




for (i in 1:ncol(clinic)) {
  pat <- toupper(colnames(clinic)[i])
  col <- grep(paste("\\b",pat,"\\b", sep=""), cnameClinic)

  if(!length(col)==0){

    clcl[,col] <- clinic[,i]

  }
}




  if(type=="Sample Pheno"){

    if(length(which(duplicated(clcl$Samples)))==0) {

      clinic2 <-  clcl
      rownames(clinic2) <- clinic2$Samples

      Metadata$Sample.pheno <- clinic2[colnames(Metadata[which(str_detect(names(Metadata), c("Raw", "matrix")))][[1]]),]

    } else {

      cl_rolled <- clcl %>%

      # create groups by name
      group_by(Samples) %>%

      summarise(across(everything(), ~paste0(unique(na.omit(.x)), collapse = ";")))
      cl_rolled <- as.data.frame(cl_rolled)
      rownames(cl_rolled) <- cl_rolled$Samples

      Metadata$Sample.pheno <- cl_rolled
      }



    } else if(type=="Patients' clinical data")
      {

      if(length(which(duplicated(clcl$Patients)))==0) {
        clinic2 <-  clcl
        rownames(clinic2) <- clinic2$Patients

        Metadata$Patient.clinic <- clinic2

        } else {




          cl_rolled <- clcl %>%

            # create groups by name
            group_by(Patients) %>%

            summarise(across(everything(), ~paste0(unique(na.omit(.x)), collapse = ";")))

          isNA <- which(is.na( cl_rolled$Patients))

          if(length(isNA)>0){ cl_rolled <- as.data.frame(cl_rolled[-isNA,])  }



          rownames(cl_rolled) <- cl_rolled$Patients

          Metadata$Patient.clinic <- cl_rolled

          }


      } else {
  stop("Choose type = 'Sample Pheno' or 'Patients' clinical data' ")}



  return(Metadata)

}

