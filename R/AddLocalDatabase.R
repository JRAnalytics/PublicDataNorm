#' AddLocalDatabase
#'
#' @param Metadata Metadata object
#' @param Normalization.Method Defaul = NA
#' @param Technology Microarray, RNAseq, Single cell etc
#' @param Platform Brand of technology
#' @param Run.spec Deppenes, reads, etc etc
#' @param First.Author Defaul = NA
#' @param DOI doi of article if exist. Defaul = NA
#' @param project project name
#' @param Comment specify a comment for this export
#' @param User User who export the cleaned Metadata object
#' @return a .text tab delimited database summary
#' @import data.table
#' @export
#'
#' @examples "non"
AddLocalDatabase <- function(Metadata,
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

  Local.Data.base.Path <- attributes(Metadata)$File.path$Parent

  lf <- list.files(Local.Data.base.Path)

  if(is.null(attributes(Metadata)$Version)){ Version = "V1"} else {   Version <- attributes(Metadata)$Version }

  NBS <- which(attributes(Metadata)$Data.Type=="SamplesAnnot" & attributes(Metadata)$Export=="Yes")
  NBS = NBS[1]

  NBP <- which(attributes(Metadata)$Data.Type=="Clinic" & attributes(Metadata)$Export=="Yes")
  NBP = NBP[1]




  if(is.na(NBS)){ Nsamples=0 } else{
    Nsamples <- nrow(Metadata[[NBS[1]]])}

  if(is.na(NBP)){ Npatient = Nsamples}else{
    Npatient <- nrow(Metadata[[NBP[1]]])}

  if(length(which(attributes(Metadata)$Data.Type=="Count" & attributes(Metadata)$Export=="Yes" ))==0){ RawGenes = 0} else {

    NB.raw.mat<- which(attributes(Metadata)$Data.Type=="Count" & attributes(Metadata)$Export=="Yes" )
    RawGenes <- nrow(Metadata[[NB.raw.mat[1]]])}

  if(length(which(attributes(Metadata)$Data.Type=="Count" & attributes(Metadata)$Export=="No" ))==0){NormGenes=0}else{
    NB.norm.mat <- which(attributes(Metadata)$Data.Type=="Count" & attributes(Metadata)$Export=="No")
    NormGenes <- nrow(Metadata[[NB.norm.mat[1]]])}


  if(!is.na(NBS)){
  if(all(is.na(Metadata[[NBS[1]]]$SamplePathologicalState))){tumor <- nrow(Metadata[[NBS[1]]]) } else {
    tumor <- length(which(str_detect(toupper(Metadata[[NBS[1]]]$SamplePathologicalState),"TUM|PRIMARY|CARCINO")))
    if(tumor==0 & !all(is.na(Metadata[[NBS[1]]]$SamplePathologicalState))) {

      normal <- length(which(str_detect(toupper(Metadata[[NBS[1]]]$SamplePathologicalState),"NORM|HEAL")))
      met <- length(which(str_detect(toupper(Metadata[[NBS[1]]]$SamplePathologicalState),"MET")))
      na <- length(which(is.na(Metadata[[NBS[1]]]$SamplePathologicalState)))
      tumor = nrow(Metadata[[NBS[1]]])-normal-met-na
    }

  }

  if(all(is.na(Metadata[[NBS[1]]]$SamplePathologicalState))){normal <- 0 } else {
    normal<- length(which(str_detect(toupper(Metadata[[NBS[1]]]$SamplePathologicalState),"NORM|HEAL")))

    if(normal==0 & !all(is.na(Metadata[[NBS[1]]]$SamplePathologicalState))) {


      met <- length(which(str_detect(toupper(Metadata[[NBS[1]]]$SamplePathologicalState),"MET")))
      na <- length(which(is.na(Metadata[[NBS[1]]]$SamplePathologicalState)))
      normal = nrow(Metadata[[NBS[1]]])-tumor-met-na
    }
  }}else {
    met <- NA
    na <- NA
    normal = NA
    tumor = NA

  }






  if(!is.na(NBS)){

    if(!is.null(Metadata[[NBS[1]]][,"HadTreatment"])) {

      TTT <- length(which(str_detect(toupper(Metadata[[NBS[1]]][,"HadTreatment"]),"YES|OUI|TRUE|1")))

    } else { TTT <- 0   }

  } else { TTT <- 0 }


  if( TTT==0  ){
    TTTinfo <- "No"

  } else {
    TTTinfo <- "Yes"
  }

  if(!is.na(NBP)){
    if(all(is.na(Metadata[[NBP[1]]]$OSdelay))){ OSinfo <- "No" } else { OSinfo <- "Yes" }
    if(all(is.na(Metadata[[NBP[1]]]$PFSdelay ))){ PFSinfo <- "No" } else { PFSinfo <- "Yes" } } else {
    OSinfo <- "No"
    PFSinfo <- "No" }

  if(!is.na(NBS)){
  met <- length(which(str_detect(toupper(Metadata[[NBS[1]]]$SamplePathologicalState),"MET")))}else{met = 0}

  print("Ok")

  dt <- data.frame("Project" = project,
                   "Version" = Version,
                   "Date.of.Data.Norm" = format(Sys.Date(),format = "%d/%m/%Y" ),
                   "N.Patients" = Npatient,
                   "N.Samples" = Nsamples,
                   "N.TumoralSamples" = tumor,
                   "N.NormalSamples" = normal,
                   "N.Metastasis" = met,
                   "Overall.Survival" = OSinfo ,
                   "Progression.Free.Survival" = PFSinfo,
                   "Treatment.Information" = TTTinfo,
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

  fp <- paste(c(Local.Data.base.Path,Databasename), collapse = "/")


  if(file.exists(fp)){



    x <-  as.data.frame(data.table::fread(fp))
    x <- x[order(x$Project,x$Version,decreasing = F),]

    if(length(x$Project[x$Project==project])!=0){



      if(!all((x$Project[x$Project==project]==project & x$Version[x$Project==project]==Version)==F)){

        message(paste(project,"already existing in database. Reactualising database"))

        proj <- which(x$Project==project)
        row <- which(x$Version[proj]==Version)


        x[proj[row],] <- dt

        print(dt)

      } else {
           x <- rbind(x,dt)
        print(dt)
      }}

    if(length(x$Project[x$Project==project])==0){ x <- rbind(x,dt)
    print(dt)

    }

    LF <- list.files(Local.Data.base.Path$Project.VerifiedDataset)
    if(length(LF)!=0){
      df <- file.info(list.files(Local.Data.base.Path$Project.VerifiedDataset, full.names = T))
      df$Filenames <- unlist(lapply(str_split(rownames(df),paste0(project,"/")),"[[",2))
      filename2 <- unlist(lapply(str_split( df$Filenames ,".csv"),"[[",1))
      version <- unique(str_extract(filename2,"V[0-9]*"))

    }

    if(length(version)==1){
      if(is.na(version)){version = "V1"}}

    if(!all(is.na(version))){
      if(!all(x$Version[x$Project==project]%in%version)){

        proj <- which(x$Project==project)

        if(length(unique(x[proj,]$Project))!=1){stop("Error in actualising DataBaseSummary.txt")}

        outV <- which(!x[proj,]$Version%in%version)


        message(paste("Project",project   ,"version",x[proj[outV],]$Version," is missing in 04VerifiedDataset. Removing from DataBaseSummary.txt\n"))

        x <- x[-proj[outV],]
      }}


    x <- x[order(x$Project,x$Version,decreasing = F),]
    write.table(x,fp,row.names = F, sep = "\t",dec = "." )
  } else {
    print(dt)
    write.table(dt,fp,row.names = F, sep = "\t",dec = "." )
  }







}
