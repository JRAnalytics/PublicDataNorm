#' DirCheck Create project Dir in Processes,n Precessed, Raw Clinic, Raw genomic, and references if does not exist.
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
  # Directories in "parent" Dir. Refereces, Processes, RawData and VerifiedDataSet

       file.path.References <- paste(c(unlist(strsplit(path, "/")),"03References"),collapse = "/")

      file.path.Processes <- paste(c(unlist(strsplit(path, "/")),"02Processes"),collapse = "/")
            file.path.Project.Processes <- paste(c(file.path.Processes,project),collapse = "/")

      file.path.RawData <- paste(c(unlist(strsplit(path, "/")),"01RawData"),collapse = "/")
            file.path.Project.RawData <- paste(c(file.path.RawData,project),collapse = "/")

      file.path.VerifiedDataSet <- paste(c(unlist(strsplit(path, "/")),"04VerifiedDataSet"),collapse = "/")
            file.path.Project.VerifiedDataset <- paste(c(file.path.VerifiedDataSet,project),collapse = "/")

  if(!dir.exists(file.path.Processes)){dir.create(file.path.Processes)}
  if(!dir.exists(file.path.RawData)){dir.create(file.path.RawData)}
  if(!dir.exists(file.path.Project.RawData)){dir.create(file.path.Project.RawData)}
  if(!dir.exists(file.path.References)){dir.create(file.path.References)}
  if(!dir.exists(file.path.VerifiedDataSet)){dir.create(file.path.VerifiedDataSet)}
  if(!dir.exists(file.path.Project.VerifiedDataset)){dir.create(file.path.Project.VerifiedDataset)}
  if(!dir.exists(file.path.Project.Processes)){dir.create(file.path.Project.Processes)}




  list.files.path <- list( "Parent"= file.path.parent,
                           "RawData" = file.path.RawData,
                           "Project.RawData" = file.path.Project.RawData,
                           "Processes" = file.path.Processes,
                           "Project.Processes" = file.path.Project.Processes,
                           "References"= file.path.References,
                           "VerifiedDataSet" = file.path.VerifiedDataSet,
                           "Propject.VerifiedDataset" = file.path.Project.VerifiedDataset)

  message("Creating Project directories")

  for (i in list.files.path) {

    message(paste(i, "------------dir exsist--------------",dir.exists(i)))

  }





  return(list.files.path)


  }
