#' CheckLocalDatabase
#'
#' @param Local.Data.base.Path path to parent directory of the database.
#' @param Normalization.Method Defaul = NA
#' @param Sequencing.DeepenessPerMillionReads Defaul = NA
#' @param Technology Microarray, RNAseq, Single cell etc
#' @param Platform Brand of technology
#' @param Run.spec Deppenes, reads, etc etc
#' @param First.Author Defaul = NA
#' @param DOI doi of article if exist. Defaul = NA
#' @param project project name
#' @param Comment specify a comment for this export
#' @param User User who export the cleaned Meta object
#' @return a .text tab delimited database summary
#' @import data.table
#' @export
#'
#' @examples "non"
CheckLocalDatabase <- function(Local.Data.base.Path,
                               Normalization.Method = NA,
                               Technology = NA,
                               Platform = NA,
                               Run.spec=NA ,
                               First.Author = NA,
                               DOI = NA,
                               project = NA,
                               Comment = NA,
                               User= NA){

  Databasename = "DataBaseSummary.txt"


   lf <- list.files(Local.Data.base.Path)

   if(is.null(attributes(Meta)$Version)){ Version = "V1"} else {   Version <- attributes(Meta)$Version }


  if(length(names(Meta)[(str_detect(names(Meta), "Sample"))])==0){ Nsamples=0 } else{
      Nsamples <- nrow(Meta[[which(str_detect(names(Meta), "Sample"))]])}

  if(length(names(Meta)[(str_detect(names(Meta), "Patient"))])==0){ Npatient = Nsamples}else{
    NBP <- which(attributes(Meta)$Data.Type=="Patient.Clinical.data" & attributes(Meta)$Raw.data=="No")
      Npatient <- nrow(Meta[[NBP]])}

  if(length(which(attributes(Meta)$Data.Type=="Expression.Matrix" & attributes(Meta)$Raw.data=="Yes" ))==0){ RawGenes = 0} else {

    NB.raw.mat<- which(attributes(Meta)$Data.Type=="Expression.Matrix" & attributes(Meta)$Raw.data=="Yes" )
    RawGenes <- nrow(Meta[[NB.raw.mat[1]]])}

  if(length(which(attributes(Meta)$Data.Type=="Expression.Matrix" & attributes(Meta)$Raw.data=="No" ))==0){NormGenes=0}else{
    NB.norm.mat <- which(attributes(Meta)$Data.Type=="Expression.Matrix" & attributes(Meta)$Raw.data=="No")
      NormGenes <- nrow(Meta[[NB.norm.mat[1]]])}

   NBS <- which(attributes(Meta)$Data.Type=="Samples.Clinical.data" & attributes(Meta)$Raw.data=="No")
   if(all(is.na(Meta[[NBS]]$SamplePathologicalState))){tumor <- nrow(Meta[[NBS]]) } else {
     tumor <- length(which(str_detect(toupper(Meta[[NBS]]$SamplePathologicalState),"TUM|PRIMARY")))
     if(tumor==0 & !all(is.na(Meta[[NBS]]$SamplePathologicalState))) {

       normal <- length(which(str_detect(toupper(Meta[[NBS]]$SamplePathologicalState),"NORM|HEAL")))
       meta <- length(which(str_detect(toupper(Meta[[NBS]]$SamplePathologicalState),"MET")))
       na <- length(which(is.na(Meta[[NBS]]$SamplePathologicalState)))
       tumor = nrow(Meta[[NBS]])-normal-meta-na
     }

   }

       if(all(is.na(Meta[[NBS]]$SamplePathologicalState))){normal <- 0 } else {
       normal<- length(which(str_detect(toupper(Meta[[NBS]]$SamplePathologicalState),"NORM|HEAL")))

       if(normal==0 & !all(is.na(Meta[[NBS]]$SamplePathologicalState))) {

         tumor <- length(which(str_detect(toupper(Meta[[NBS]]$SamplePathologicalState),"TUM|PRIMARY")))
         meta <- length(which(str_detect(toupper(Meta[[NBS]]$SamplePathologicalState),"MET")))
         na <- length(which(is.na(Meta[[NBS]]$SamplePathologicalState)))
         normal = nrow(Meta[[NBS]])-tumor-meta-na
       }
       }

    NBP <- which(attributes(Meta)$Data.Type=="Patient.Clinical.data" & attributes(Meta)$Raw.data=="No")

  if(length(NBP)>0){ TTT<- Meta[[NBP]] %>% subset(TreatmentInfo=="Yes")%>%nrow()
    } else { TTT <- 0 }


  if( TTT==0  ){
    TTTinfo <- "No"
    TTTtype <- NA
  } else {
    TTTinfo <- "Yes"
    TTTtype <- paste(unique(Meta[[NBP]][,c("Treatment.AdjType"  ,    "Treatment.NeoAdjType"  , "Treatment.RT"    )]),collapse = ",")
    TTTtype <- paste(unique(Meta[[NBP]][,c("Treatment.AdjType"  ,    "Treatment.NeoAdjType" ,  "Treatment.RT"    )]),collapse = ",")
  }

  if(length(NBP)>0){
  if(all(is.na(Meta[[NBP]]$OSdelay))){ OSinfo <- "No" } else { OSinfo <- "Yes" }
  if(all(is.na(Meta[[NBP]]$PFSdelay ))){ PFSinfo <- "No" } else { PFSinfo <- "Yes" } }

  else {
    OSinfo <- "No"
  PFSinfo <- "No" }



  dt <- data.frame("Project" = project,
                   "Version" = Version,
                   "Date.of.Data.Norm" = Sys.Date(),
                   "N.Patients" = Npatient,
                   "N.Samples" = Nsamples,
                   "N.TumoralSamples" = tumor,
                   "N.NormalSamples" = normal,
                   "N.Metastasis" = meta,
                   "Overall.Survival" =OSinfo ,
                   "Progression.Free.Survival" = PFSinfo,
                   "Treatment.Information" = TTTinfo,
                   "Treatment.Type" = TTTtype,
                   "N.RawGenes" = RawGenes,
                   "N.NormalizedGenes"=NormGenes,
                   "Normalization.Method" = Normalization.Method,
                   "Technology" = Technology,
                   "Platform" = Platform,
                   "Run.spec"=Run.spec ,
                   "First.Author" = First.Author,
                   "Article.DOI" = DOI,
                   "Comment" = Comment,
                   "User"= User)





  message("Adding to data bases : ")
  print(dt)
  fp <- paste(c(Local.Data.base.Path,Databasename), collapse = "/")


  if(file.exists(fp)){


    x <-  as.data.frame(data.table::fread(fp))


    if(x$Project==project & x$Version[x$Project==project]==attributes(Meta)$Version){

      message(paste(project,"already existing in database. Reactualising database"))

      row <- max(which(x$Project==project & x$Version[x$Project==project]==attributes(Meta)$Version))

      x[row,] <- dt

    } else { x <- rbind(x,dt) }


    x <- x[order(x$Project,x$Version,decreasing = F),]

    write.table(x,fp,row.names = F, sep = "\t",dec = "." )
  } else {

    write.table(dt,fp,row.names = F, sep = "\t",dec = "." )
  }









}
