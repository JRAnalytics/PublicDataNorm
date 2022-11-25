#' CheckLocalDatabase
#'
#' @param Local.Data.base.Path path to parent directory of the database.
#' @param Databasename name of the database.txt summary file
#' @return a .text tab delimited database summary
#' @import data.table
#' @export
#'
#' @examples "non"
CheckLocalDatabase <- function(Local.Database.Path, Databasename = "DataBaseSummary.txt"){



  lf <- list.files(Local.Data.base.Path)
  lf <- lf[!lf%in%"Lexic"]


  Npatient <- nrow(Meta[[which(str_detect(names(Meta), "Patient"))]])
  Nsamples <- nrow(Meta[[which(str_detect(names(Meta), "Sample"))]])
  RawGenes <- nrow(Meta[[which(str_detect(names(Meta), "RawCount"))]])
  NormGenes <- nrow(Meta[[which(str_detect(names(Meta), "Normalized"))]])

  tumor <- Meta[[which(str_detect(names(Meta), "Sample"))]] %>% subset(SamplePathologicalState=="Tumor")%>%nrow()
  normal<- Meta[[which(str_detect(names(Meta), "Sample"))]] %>% subset(SamplePathologicalState=="Normal")%>%nrow()


  TTT<- Meta[[which(str_detect(names(Meta), "Patient"))]] %>% subset(TreatmentInfo=="Yes")%>%nrow()

  colNamesPatientsClinic <- colnames(Meta[[which(str_detect(names(Meta), "Patient"))]])

  if( TTT==0  ){
    TTTinfo <- "No"
    TTTtype <- NA
  } else {
    TTTinfo <- "Yes"
    TTTtype <- paste(unique(Meta[[which(str_detect(names(Meta), "Patient"))]][,c("Treatment.AdjType"  ,    "Treatment.NeoAdjType"  , "Treatment.RT"    )]),collapse = ",")
    TTTtype <- paste(unique(Meta[[which(str_detect(names(Meta), "Patient"))]][,c("Treatment.AdjType"  ,    "Treatment.NeoAdjType" ,  "Treatment.RT"    )]),collapse = ",")
  }

  if(all(is.na(Meta[[which(str_detect(names(Meta), "Patient"))]]$OSdelay))){ OSinfo <- "No" } else { OSinfo <- "Yes" }
  if(all(is.na(Meta[[which(str_detect(names(Meta), "Patient"))]]$PFSdelay ))){ PFSinfo <- "No" } else { PFSinfo <- "Yes" }





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
                   "Normalization.Method" = "UpperQuartile",
                   "Sequençing.DeepenessPerMillionReads" = 30,
                   "Sequençing.LengthPer.pb" = 100,
                   "Sequençing.Run" = "Single-end",
                   "Sequençing.Platform" ="Illumina HiSeq 2000" ,
                   "First.Author" = "Carlo Maurer",
                   "Article.DOI" = "http://dx.doi.org/10.1136/gutjnl-2018-317706",
                   row.names = project)
  dt

  fp <- paste(c(Local.Data.base.Path,Databasename), collapse = "/")


  if(file.exists(fp)){


    x <-  fread(fp)

    x <- rbind(x,dt)

    write.table(x,fp,row.names = F, sep = "\t",dec = "." )
  } else {
    write.table(dt,fp,row.names = F, sep = "\t",dec = "." )
  }









}
