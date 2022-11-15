#' removefilesrecursively remove donwloaded files from the TCGA
#'
#' @param files list files in dir to remove
#' @param remove.files.recursively function
#' @return nothing
#' @examples   "none"
#' @export
#'
#'
removefilesrecursively <- function(files){
  files2rm <- dirname(files)
  unlink(files2rm,recursive = TRUE)
  files2rm <- dirname(files2rm) # data category
  if(length(list.files(files2rm)) == 0) remove.files.recursively(files2rm)
}
