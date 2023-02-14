#' CheckLocalDatabase
#'
#' @param Local.Data.base.Path path to parent directory of the database.
#' @param Databasename name of the database.txt summary file
#' @param Normalization.Method Defaul = NA
#' @param Sequençing.DeepenessPerMillionReads Defaul = NA
#' @param Sequençing.LengthPer.pb Defaul = NA
#' @param Sequençing.Platform Defaul = NA
#' @param First.Author Defaul = NA
#' @param DOI doi of article if exist. Defaul = NA
#' @param project project name
#' @return a .text tab delimited database summary
#' @import data.table
#' @export
#'
#' @examples "non"
CheckLocalDatabase <- function(Local.Data.base.Path,
                               Databasename = "DataBaseSummary.txt",
                               Normalization.Method = NA,
                               Sequençing.DeepenessPerMillionReads= NA,
                               Sequençing.LengthPer.pb = NA,
                               Sequençing.Run = NA,
                               Sequençing.Platform = NA ,
                               First.Author = NA,
                               DOI = NA,
                              project = NA){



   lf <- list.files(Local.Data.base.Path)



  if(length(names(Meta)[(str_detect(names(Meta), "Sample"))])==0){ Nsamples=0 } else{
      Nsamples <- nrow(Meta[[which(str_detect(names(Meta), "Sample"))]])}

  if(length(names(Meta)[(str_detect(names(Meta), "Patient"))])==0){ Npatient = Nsamples}else{
      Npatient <- nrow(Meta[[which(str_detect(names(Meta), "Patient"))]])}

  if(length(names(Meta)[(str_detect(names(Meta), "Raw.matrix"))])==0){ RawGenes = 0}else{
    RawGenes <- nrow(Meta[[which(str_detect(names(Meta), "Raw.matrix"))]])}

  if(length(names(Meta)[(str_detect(names(Meta), "Normalized"))])==0){NormGenes=0}else{

      NormGenes <- nrow(Meta[[which(str_detect(names(Meta), "Normalized"))]])}

  tumor <- Meta[[which(str_detect(names(Meta), "Sample"))]] %>% subset(SamplePathologicalState=="Tumor")%>%nrow()
  normal<- Meta[[which(str_detect(names(Meta), "Sample"))]] %>% subset(SamplePathologicalState=="Normal")%>%nrow()


  if(length(which(str_detect(names(Meta), "Patient")))>0){ TTT<- Meta[[which(str_detect(names(Meta), "Patient"))]] %>% subset(TreatmentInfo=="Yes")%>%nrow()
    } else { TTT <- 0 }


  if( TTT==0  ){
    TTTinfo <- "No"
    TTTtype <- NA
  } else {
    TTTinfo <- "Yes"
    TTTtype <- paste(unique(Meta[[which(str_detect(names(Meta), "Patient"))]][,c("Treatment.AdjType"  ,    "Treatment.NeoAdjType"  , "Treatment.RT"    )]),collapse = ",")
    TTTtype <- paste(unique(Meta[[which(str_detect(names(Meta), "Patient"))]][,c("Treatment.AdjType"  ,    "Treatment.NeoAdjType" ,  "Treatment.RT"    )]),collapse = ",")
  }

  if(length(which(str_detect(names(Meta), "Patient")))>0){
  if(all(is.na(Meta[[which(str_detect(names(Meta), "Patient"))]]$OSdelay))){ OSinfo <- "No" } else { OSinfo <- "Yes" }
  if(all(is.na(Meta[[which(str_detect(names(Meta), "Patient"))]]$PFSdelay ))){ PFSinfo <- "No" } else { PFSinfo <- "Yes" } }

  else {
    OSinfo <- "No"
  PFSinfo <- "No" }



  dt <- data.frame("Project" = project,
                   "Date.of.Data.Norm" = Sys.Date(),
                   "N.Patients" = Npatient,
                   "N.Samples" = Nsamples,
                   "N.TumoralSamples" = tumor,
                   "N.NormalSamples" = normal,
                   "Overall.Survival" =OSinfo ,
                   "Progression.Free.Survival" = PFSinfo,
                   "Treatment.Information" = TTTinfo,
                   "Treatment.Type" = TTTtype,
                   "N.RawGenes" = RawGenes,
                   "N.NormalizedGenes"=NormGenes,
                   "Normalization.Method" = Normalization.Method,
                   "Sequençing.DeepenessPerMillionReads" = Sequençing.DeepenessPerMillionReads,
                   "Sequençing.LengthPer.pb" = Sequençing.LengthPer.pb,
                   "Sequençing.Run" = Sequençing.Run,
                   "Sequençing.Platform" =Sequençing.Platform ,
                   "First.Author" = First.Author,
                   "Article.DOI" = DOI,
                   row.names = project)





  message("Adding to data bases : ")
  print(dt)
  fp <- paste(c(Local.Data.base.Path,Databasename), collapse = "/")


  if(file.exists(fp)){


    x <-  as.data.frame(data.table::fread(fp))
    rownames(x) <- x$Project

    if(project%in%rownames(x)){

      message(paste(project,"already existing in database. Reactualising database"))
      x[project,] <- dt

    } else { x <- rbind(x,dt) }


    write.table(x,fp,row.names = F, sep = "\t",dec = "." )
  } else {
    write.table(dt,fp,row.names = F, sep = "\t",dec = "." )
  }









}
