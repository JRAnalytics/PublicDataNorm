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
  if(is.null(PatientsLexic)&is.null(SamplesLexic)){stop("A PatientsLexic or SamplesLexic is mandatory for data cleaning")}


  if(!is.null(SamplesLexic) & attributes(Metadata)$Omics.type!="Single.Cell"){
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

 if(!all(is.null(PatientsExportname),is.null(PatientsLexic),is.null(PatientsAnnotToClean))){
    if(is.null(PatientsAnnotToClean) & attributes(Metadata)$Omics.type!="Single.Cell"   ){

      RP = which(attributes(Metadata)$Data.Type=="Clinic" &attributes(Metadata)$Cleaned=="No")
      message("Creating a Patient clinical table from Samples annotation")
      if(length(RP)>0){stop("Raw Patients clinical data found in Metaobject, Specify PatientsAnnotToClean")}
      if(is.null(PatientsExportname)){stop("Specify PatientsExportname.")}
      if(is.null(PatientsLexic)){stop("PatientsLexic is mandatory.")}


      Metadata <- CleaningClinic(Metadata = Metadata,
                                 Lexic = PatientsLexic,
                                 type = "Patients",
                                 ClinicToClean = SamplesAnnotToClean,
                                 exportname = PatientsExportname,
                                 FilterPatients =  FilterSP,
                                 FilterSamples = F,
                                 force.replace = force.replace,
                                 CleanFromOtherType = T)



      if(keep.all.column==T){
      Metadata <- CleaningClinic(Metadata = Metadata,
                                 Lexic = PatientsLexic,
                                 type = "Patients",
                                 ClinicToClean = SamplesAnnotToClean,
                                 exportname = paste0(PatientsExportname,".fullCol"),
                                 FilterPatients =  FilterSP,
                                 FilterSamples = F,
                                 force.replace = force.replace,
                                 CleanFromOtherType = T,
                                 keep.all.column = T)

      }





      } else{

      if(is.null(PatientsLexic)){stop("A PatientsLexic is mandatory for data cleaning")}


    ClinicRaw = which(attributes(Metadata)$Data.Type=="Clinic" & attributes(Metadata)$Cleaned == "No")
    SamAnnotRaw = which(attributes(Metadata)$Data.Type=="SamplesAnnot" & attributes(Metadata)$Cleaned == "No")
    if(length(ClinicRaw)!=0){

      if(length(SamAnnotRaw)!=0){

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


  }}}

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

    Metadata[[i]] = Metadata[[i]][object[,colT],]

    attributes(Metadata)$Cleaned[i] = "Yes"

    }


    Metadata[[g]] = object


}

  if(attributes(Metadata)$Omics.type=="Single.Cell"){

    if(is.null(SamplesLexic)){stop("You can set SamplesLexic with a Lexic for Cells Annotation cleaning.")}
    if(!"SamplesAnnot"%in%attributes(Metadata)$Data.Type & !is.null(SamplesLexic)){

      SamAnnotRaw = which(attributes(Metadata)$Data.Type=="CellsAnnot" & attributes(Metadata)$Cleaned == "No")
      pp = which(attributes(Metadata)$Data.Type=="Clinic" &attributes(Metadata)$Cleaned=="Yes")

      Metadata <- CleaningClinic(Metadata = Metadata,
                                 Lexic = SamplesLexic,
                                 type = "Cells",
                                 ClinicToClean = names(Metadata)[SamAnnotRaw],
                                 exportname = paste("Cells.Annotation"),
                                 FilterPatients =FilterSP,
                                 FilterSamples = F,
                                 CleanFromOtherType = F,
                                 force.replace = force.replace,
                                 keep.all.column =F )

      cellID = Metadata$Cells.Annotation$CellsBarcode

      if(length(pp)>0){
        clinic = Metadata[[pp]]
        for (z in clinic$PatientsID){

          Metadata$Cells.Annotation$PatientsID=ifelse(str_detect(pattern = paste0(z,"_"), cellID),z,Metadata$Cells.Annotation$PatientsID)
        }

      }


      if(keep.all.column==T){
      Metadata <- CleaningClinic(Metadata = Metadata,
                                 Lexic = SamplesLexic,
                                 type = "Cells",
                                 ClinicToClean = names(Metadata)[SamAnnotRaw],
                                 exportname = "Cells.Annotation.Full",
                                 FilterPatients =FilterSP,
                                 FilterSamples = F,
                                 CleanFromOtherType = F,
                                 force.replace = force.replace,
                                 keep.all.column =keep.all.column )


      cellID = Metadata$Cells.Annotation$CellsBarcode

      if(length(pp)>0){
        clinic = Metadata[[pp]]
        for (z in clinic$PatientsID){
          Metadata$Cells.Annotation.Full$PatientsID=ifelse(str_detect(pattern = paste0(z,"_"), cellID),z,Metadata$Cells.Annotation.Full$PatientsID)
        }

      }



      }


    }

    ClinicRaw = which(attributes(Metadata)$Data.Type=="Clinic" &attributes(Metadata)$Cleaned=="No")
    SampleAnnotRaw = which(attributes(Metadata)$Data.Type=="SamplesAnnot" &attributes(Metadata)$Cleaned=="No")

    if(length(ClinicRaw)<1 & !is.null(PatientsLexic)){


      SamAnnotRaw = which(attributes(Metadata)$Data.Type=="SamplesAnnot" & attributes(Metadata)$Cleaned == "No")


      Metadata <- CleaningClinic(Metadata = Metadata,
                                 Lexic = PatientsLexic,
                                 type = "Patients",
                                 ClinicToClean = names(Metadata)[SamAnnotRaw],
                                 exportname = PatientsExportname,
                                 FilterPatients =FilterSP,
                                 FilterSamples = FilterSP,
                                 CleanFromOtherType = T,
                                 force.replace = force.replace,
                                 keep.all.column =F )}



    if(length(ClinicRaw)>0 & !is.null(PatientsLexic)){

      CellNAnotRaw = which(attributes(Metadata)$Data.Type=="CellsAnnot" & attributes(Metadata)$Cleaned == "No")


      Metadata <- CleaningClinic(Metadata = Metadata,
                                 Lexic = SamplesLexic,
                                 type = "Cells",
                                 ClinicToClean = names(Metadata)[CellNAnotRaw],
                                 exportname = paste("Cells.Annotation"),
                                 FilterPatients =FilterSP,
                                 FilterSamples = FilterSP,
                                 CleanFromOtherType = F,
                                 force.replace = force.replace,
                                 keep.all.column =F )


      Metadata <- CleaningClinic(Metadata = Metadata,
                                 Lexic = PatientsLexic,
                                 type = "Patients",
                                 ClinicToClean = names(Metadata)[ClinicRaw],
                                 exportname = PatientsExportname,
                                 FilterPatients =FilterSP,
                                 FilterSamples = F,
                                 CleanFromOtherType = F,
                                 force.replace = force.replace,
                                 keep.all.column =F )
      if(keep.all.column==T){

        Metadata <- CleaningClinic(Metadata = Metadata,
                                   Lexic = PatientsLexic,
                                   type = "Patients",
                                   ClinicToClean = names(Metadata)[ClinicRaw],
                                   exportname = paste0(PatientsExportname,".Full"),
                                   FilterPatients =FilterSP,
                                   FilterSamples = F,
                                   CleanFromOtherType = F,
                                   force.replace = force.replace,
                                   keep.all.column =T )



      }

    }

    if(length(SampleAnnotRaw)>0 & !is.null(SamplesLexic)){

      CellNAnotRaw = which(attributes(Metadata)$Data.Type=="CellsAnnot" & attributes(Metadata)$Cleaned == "No")


      Metadata <- CleaningClinic(Metadata = Metadata,
                                 Lexic = SamplesLexic,
                                 type = "Cells",
                                 ClinicToClean = names(Metadata)[CellNAnotRaw],
                                 exportname = paste("Cells.Annotation"),
                                 FilterPatients =FilterSP,
                                 FilterSamples = FilterSP,
                                 CleanFromOtherType = F,
                                 force.replace = force.replace,
                                 keep.all.column =F )


      Metadata <- CleaningClinic(Metadata = Metadata,
                                 Lexic = SamplesLexic,
                                 type = "Samples",
                                 ClinicToClean = names(Metadata)[SampleAnnotRaw],
                                 exportname = SamplesExportname,
                                 FilterPatients =F,
                                 FilterSamples = FilterSP,
                                 CleanFromOtherType = F,
                                 force.replace = force.replace,
                                 keep.all.column =F )

      if(keep.all.column==T){
        Metadata <- CleaningClinic(Metadata = Metadata,
                                   Lexic = SamplesLexic,
                                   type = "Samples",
                                   ClinicToClean = names(Metadata)[SampleAnnotRaw],
                                   exportname = paste0(SamplesExportname,".Full"),
                                   FilterPatients =F,
                                   FilterSamples = FilterSP,
                                   CleanFromOtherType = F,
                                   force.replace = force.replace,
                                   keep.all.column =T )

      }

    }

    }



  return(Metadata)

}
