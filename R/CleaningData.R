#' CleaningData from Metaobject
#'
#' @param Metadata  a Metadata Object to fill
#' @param PatientsLexic  Lexic to use for cleaning. Created from CreateLexic function .
#' @param PatientsAnnotToClean a character string of name of clinical data to clean.
#' @param PatientsExportname  name to apply in Metadata object list
#' @param SamplesLexic  Lexic to use for cleaning. Created from CreateLexic function .
#' @param SamplesAnnotToClean a character string of name of Samples annotation data to clean.
#' @param SamplesExportname  name to apply in Metadata object list
#' @param FilterSP default F, if T, keep only retrieved samples in SamplesAnnotation file
#' @param force.replace set as F. T : replace an already object with the same name
#' @param keep.all.column default F, if T, copy all column from clinic.
#' @param FilterGenes Filter for genes found in geneannotation file and rownames of matrices.
#' @importFrom utils menu
#' @import dplyr
#' @return a meataobject
#' @export
#'
#' @examples "non"
CleaningData = function(Metadata = NULL,
                        PatientsLexic = NULL,
                        PatientsAnnotToClean = NULL,
                        PatientsExportname = NULL,
                        SamplesLexic = NULL,
                        SamplesAnnotToClean = NULL,
                        SamplesExportname = NULL,
                        FilterSP =  F,
                        force.replace = F,
                        keep.all.column = F,
                        FilterGenes = F){



  if(is.null(Metadata)){stop("No Metadata found.")}
  if(is.null(PatientsLexic)&is.null(SamplesLexic)){stop("a PatientsLexic or SamplesLexic is mandatory for data cleaning")}


  if(!is.null(SamplesLexic)){
      if(is.null(SamplesAnnotToClean)){stop("SamplesAnnotToClean has to be specify.")}

    SamAnnotRaw = which(attributes(Metadata)$Data.Type=="SamplesAnnot" & attributes(Metadata)$Cleaned == "No")
    ClinicRaw = which(attributes(Metadata)$Data.Type=="Clinic" & attributes(Metadata)$Cleaned == "No")
    if(length(SamAnnotRaw)!=0){


    Metadata <- CleaningClinic(Metadata = Metadata,
                             Lexic = SamplesLexic,
                             type = "Samples",
                             ClinicToClean = SamplesAnnotToClean,
                             exportname = SamplesExportname,
                             FilterSamples =  FilterSP,
                             force.replace = force.replace )

  if(keep.all.column==T){
    Metadata <- CleaningClinic(Metadata = Metadata,
                               Lexic = SamplesLexic,
                               type = "Samples",
                               ClinicToClean = SamplesAnnotToClean,
                               exportname = paste0(SamplesExportname,".fullCol"),
                               FilterSamples =  FilterSP,
                               force.replace = force.replace,
                               keep.all.column =keep.all.column )}

    }
  }


  if(!is.null(PatientsLexic)){
    if(is.null(PatientsAnnotToClean)){stop("PatientsAnnotToClean has to be specify.")}

    ClinicRaw = which(attributes(Metadata)$Data.Type=="Clinic" & attributes(Metadata)$Cleaned == "No")
    SamAnnotRaw = which(attributes(Metadata)$Data.Type=="SamplesAnnot" & attributes(Metadata)$Cleaned == "No")
    if(length(ClinicRaw)!=0){

      if(length(SamAnnotRaw)!=0){

       SID =  Metadata[[SamAnnotRaw[1]]]$SamplesID

       if(is.null( Metadata[[SamAnnotRaw[1]]]$PatientsID)){message(paste("No PatientsID column found in", names(Metadata)[SamAnnotRaw[1]], "\nCleaning may be sub optimal."))}

        PID =  Metadata[[ClinicRaw[1]]]$PatientsID

        if(is.null( Metadata[[ClinicRaw[1]]]$SamplesID)){message(paste("No SamplesID column found in", names(Metadata)[ClinicRaw[1]], "\nCleaning may be sub optimal."))}

          if(length(which(PID%in%SID))==length(SID)){
            message("All PatientsID are found as SamplesID")
            PatientsLexic <- AddKeyLexic(lexic = PatientsLexic, key ="SamplesID"  ,value = "SamplesID" )

            Metadata[[ClinicRaw[1]]]$SamplesID =  Metadata[[ClinicRaw[1]]]$PatientsID
            PatientsLexic <- AddKeyLexic(lexic = PatientsLexic, key ="SamplesID"  ,value = "SamplesID" )

          } else {

            if(length(which(PID%in%SID))<length(SID)){
              message("Not all PatientsID are found as SamplesID")
              Metadata[[ClinicRaw[1]]]= Metadata[[ClinicRaw[1]]][which(PID%in%SID),]

              Metadata[[ClinicRaw[1]]]$SamplesID =  Metadata[[ClinicRaw[1]]]$PatientsID
              PatientsLexic <- AddKeyLexic(lexic = PatientsLexic, key ="SamplesID"  ,value = "SamplesID" )}



           else {

            if(length(which(PID%in%SID))==0){
              message(paste("No PatientsID are found as SamplesID. Try setting 'PatientID' in", names(Metadata)[ClinicRaw[1]]))


            }

        }}

    Metadata <- CleaningClinic(Metadata = Metadata,
                               Lexic = PatientsLexic,
                               type = "Patients",
                               ClinicToClean = PatientsAnnotToClean,
                               exportname = PatientsExportname,
                               FilterPatients =  FilterSP,
                               FilterSamples = F,
                               force.replace = force.replace )

    if(keep.all.column==T){
      Metadata <- CleaningClinic(Metadata = Metadata,
                                 Lexic = PatientsLexic,
                                 type = "Patients",
                                 ClinicToClean = PatientsAnnotToClean,
                                 exportname = paste0(PatientsExportname,".fullCol"),
                                 FilterPatients =FilterSP,
                                 FilterSamples = F,
                                 force.replace = force.replace,
                                 keep.all.column =keep.all.column )}
      } else {



        Metadata <- CleaningClinic(Metadata = Metadata,
                                   Lexic = PatientsLexic,
                                   type = "Patients",
                                   ClinicToClean = PatientsAnnotToClean,
                                   exportname = PatientsExportname,
                                   FilterPatients =  FilterSP,
                                   FilterSamples = F,
                                   force.replace = force.replace )

        if(keep.all.column==T){
          Metadata <- CleaningClinic(Metadata = Metadata,
                                     Lexic = PatientsLexic,
                                     type = "Patients",
                                     ClinicToClean = PatientsAnnotToClean,
                                     exportname = paste0(PatientsExportname,".fullCol"),
                                     FilterPatients =FilterSP,
                                     FilterSamples = F,
                                     force.replace = force.replace,
                                     keep.all.column =keep.all.column )}










        }


  }}

if(FilterGenes == T){


  g =  which(attributes(Metadata)$Data.Type=="geneAnnot")[1]
  m <- which(attributes(Metadata)$Data.Type=="Count")
  object = Metadata[[g]]
  if(!is.na(g)){geneAnnot = as.matrix(Metadata[[g]])}

  for (i in m){

    gene = rownames(Metadata[[i]])
    sel = which(geneAnnot %in% gene)
    colT = which(apply(geneAnnot, 2, function(x) which(x %in% geneAnnot[sel[1]]))>0)[1]

    object <- object[object[,colT]%in%gene,]

    if(length(object[,colT])!=length(unique(object[,colT]))){object <- object[!duplicated(object[,colT]),]}

    Metadata[[i]] = Metadata[[i]][rownames(object),]

    attributes(Metadata)$Cleaned[i] = "Yes"

    }


    Metadata[[g]] = object


}





  return(Metadata)

}
