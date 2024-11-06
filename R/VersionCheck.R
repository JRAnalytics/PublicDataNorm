#' VersionCheck : Asks DataBaseSumary.txt file for the oldest project's version.
#'
#' @param Metadata the Metaobject of the project
#' @param path file path to DataBaseSumary.txt
#'
#' @return Return a Versionned metadata object
#' @export
#'
#' @examples "non"
VersionCheck = function(Metadata, path = NULL){

  if(is.null(Metadata)){ stop("A Metadata object created by CreateDataset(project,path) is mandatory")}
  if(is.null(path)){
    path <- attributes(Metadata)$File.path$Parent
    }

  Listfiles = list.files(path)
  project =  attributes(Metadata)$Project


  if(!"DataBaseSummary.txt" %in%  Listfiles){stop("No DataBaseSummary.txt found in file path. Change to the database path.")} else{


  if(is.null(  attributes(Metadata)$Version)){ message("The Metaobject is not 'Versionned', checking in the DataBaseSummary.txt if project exists.")}

    fp <- paste(c(path,"DataBaseSummary.txt"), collapse = "/")

      x <-  as.data.frame(data.table::fread(fp))
      x <- x[order(x$Project,x$Version,decreasing = F),]

      if(length(x$Project[x$Project==project])!=0){

          message(paste(project,"already existing in database. Extracting oldest version"))

          proj <- which(x$Project==project)

          oldversion = unique(x[proj,]$Version)
          oldversion=as.numeric(gsub("V", "",oldversion))
          oldversion = max(oldversion)



          newversion = oldversion+1
          message(paste0(project," Version is: V",oldversion,".\nIncrementing project to V",newversion))

      } else {

        message(paste(project,"is not in the database. Conragts to create the first version of the cleaned project."))

        newversion = 1

          }


      attributes(Metadata)$Version = paste0("V",newversion)


    return(Metadata)
  }
}
