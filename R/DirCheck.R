
#' DirCheck Create project Dir in Script,n Precessed, Raw Clinic, Raw genomic, and references if does not exist.
#'
#' @param project a character string : project name
#' @importFrom rstudioapi getSourceEditorContext
#' @return a file path list of dir project in ~parent/dir/project file
#' @export
#'
#' @examples DirCheck("Maurer")
#'
#'
#'
DirCheck <- function(project){


  project <- project

  rstudioapi::documentSave(paste0(project,".R"))

  path <- rstudioapi::getSourceEditorContext()$path
  NScript <- which(unlist(strsplit(path, "/"))=="Script")

  file.path.parent <- paste(unlist(strsplit(path, "/"))[1:(NScript-1)],collapse = "/")

  file.path.Script <- paste(c(unlist(strsplit(path, "/"))[1:(NScript)],project),collapse = "/")
  file.path.Raw.genomic <- paste(c(unlist(strsplit(path, "/"))[1:(NScript-1)],"RawGenomic",project),collapse = "/")
  file.path.Raw.clinic <- paste(c(unlist(strsplit(path, "/"))[1:(NScript-1)],"RawClinic",project),collapse = "/")
  file.path.Processed <- paste(c(unlist(strsplit(path, "/"))[1:(NScript-1)],"Processed",project),collapse = "/")
  file.path.References <- paste(c(unlist(strsplit(path, "/"))[1:(NScript-1)],"References",project),collapse = "/")

  if(!dir.exists(file.path.Script)){dir.create(file.path.Script)}
  if(!dir.exists(file.path.Raw.genomic)){dir.create(file.path.Raw.genomic)}
  if(!dir.exists(file.path.Raw.clinic)){dir.create(file.path.Raw.clinic)}
  if(!dir.exists(file.path.Processed)){dir.create(file.path.Processed)}
  if(!dir.exists(file.path.References)){dir.create(file.path.References)}

  list.files.path <- list( "Parent"= file.path.parent ,
                           "Script" = file.path.Script,
                           "RawGenomic" = file.path.Raw.genomic,
                           "RawClinic" = file.path.Raw.clinic,
                           "Processed" = file.path.Processed,
                           "References"= file.path.References)

  return(list.files.path)


  }
