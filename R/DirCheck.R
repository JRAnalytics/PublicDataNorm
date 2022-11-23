#' DirCheck Create project Dir in Script,n Precessed, Raw Clinic, Raw genomic, and references if does not exist.
#'
#' @param project a character string : project name
#' @param path dir path of Parent directory
#' @importFrom rstudioapi getSourceEditorContext
#' @return a file path list of dir project in ~parent/dir/project file
#' @export
#'
#' @examples "none"
#'
#'
#'
DirCheck <- function(project,path){

  file.path.parent <- paste(unlist(strsplit(path, "/")),collapse = "/")
  file.path.Project <- paste(c(unlist(strsplit(path, "/")),project),collapse = "/")
  file.path.Script <- paste(c(unlist(strsplit(path, "/")),project,"Script"),collapse = "/")
  file.path.RawDataDump <- paste(c(unlist(strsplit(path, "/")),project,"RawDataDump"),collapse = "/")
  file.path.RawPhenoAnnotation <- paste(c(unlist(strsplit(path, "/")),project,"RawPhenoAnnotation"),collapse = "/")
  file.path.References <- paste(c(unlist(strsplit(path, "/")),project,"References"),collapse = "/")
  file.path.PipelineDump <- paste(c(unlist(strsplit(path, "/")),project,"PipelineDump"),collapse = "/")
  file.path.VerifiedDataSet <- paste(c(unlist(strsplit(path, "/")),project,"VerifiedDataSet"),collapse = "/")

  if(!dir.exists(file.path.Project)){dir.create(file.path.Project)}
  if(!dir.exists(file.path.Script)){dir.create(file.path.Script)}
  if(!dir.exists(file.path.RawDataDump)){dir.create(file.path.RawDataDump)}
  if(!dir.exists(file.path.RawPhenoAnnotation)){dir.create(file.path.RawPhenoAnnotation)}
  if(!dir.exists(file.path.References)){dir.create(file.path.References)}
  if(!dir.exists(file.path.PipelineDump)){dir.create(file.path.PipelineDump)}
  if(!dir.exists(file.path.VerifiedDataSet)){dir.create(file.path.VerifiedDataSet)}


  list.files.path <- list( "Parent"= file.path.parent,
                           "Project" = file.path.Project,
                           "Script" = file.path.Script,
                           "RawDataDump" = file.path.RawDataDump,
                           "RawPhenoAnnotation" = file.path.RawPhenoAnnotation,
                           "References"= file.path.References,
                           "PipelineDump" = file.path.PipelineDump,
                           "VerifiedDataSet" = file.path.VerifiedDataSet)

  message("Creating Project directories")

  for (i in list.files.path) {

    message(paste(i, "------------dir exsist--------------",dir.exists(i)))

  }





  return(list.files.path)


  }
