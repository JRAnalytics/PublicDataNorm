#' getProteinCodingMatrix
#'
#' @param Metadata Meta data from PublicDataNorm
#' @importFrom rlang .data
#' @return a Meta object
#' @export
#'
#' @examples "none"
getProteinCodingMatrix <- function(Metadata){


  Meta <- Metadata
  x <- subset(Meta$geneAnnotation, gene_type == "protein_coding")



  if(is.null(Meta$RawCount.proteincoding.matrix)){
  if(!all(str_detect(names(Meta), "Raw", negate = FALSE)=="FALSE")){

    if(all(rownames(Meta[[which(str_detect(names(Meta), "Raw", negate = FALSE)=="TRUE")]])%in%x$gene_names)==F) {message("Raw Matrice already with gene protein coding names as rownames")}else {

    message("Building RawCount.proteincoding matrice expression")

    y <- Meta[[which(str_detect(names(Meta), "Raw", negate = FALSE)=="TRUE")]][rownames(x),]
    y$gene_names <- x$gene_names

    y <- y %>% group_by(gene_names) %>% summarise_each(funs(round(mean(.))))
    y <- as.data.frame(y)
    rownames(y) <- y$gene_names
    y$gene_names <- NULL

    Meta$RawCount.proteincoding.matrix <- as.data.frame(y)

  }}else{message("No Raw data matrice in Meta oject existing. \n Add a matrice with AddExpressionMatrix function")}
    }else{message("Raw.proteincoding in Meta object already existing. \n Not computing code.")}

  if(is.null(Meta$TPM.proteincoding.matrix)){
  if(!all(str_detect(names(Meta), "TPM", negate = FALSE)=="FALSE")){


    if(all(rownames(Meta[[which(str_detect(names(Meta), "TPM", negate = FALSE)=="TRUE")]])%in%x$gene_names)==F) {message("TPM Matrice already with gene protein coding names as rownames")} else {


    message("Building TPM.proteincoding matrice expression")
    y <- Meta[[which(str_detect(names(Meta), "TPM", negate = FALSE)=="TRUE")]][rownames(x),]
    y$gene_names <- x$gene_names
    y <- y %>% group_by(gene_names) %>% summarise_each(funs(mean))
    y <- as.data.frame(y)
    rownames(y) <- y$gene_names
    y$gene_names <- NULL

    Meta$TPM.proteincoding.matrix <- as.data.frame(y)
  }}else{message("No TPM data matrice in Meta oject existing. \n Add a matrice with AddExpressionMatrix function")}
    }else{message("TPM.proteincoding in Meta object already existing. \n Not computing code.")}

  if(is.null(Meta$FKPM.proteincoding.matrix)){
  if(!all(str_detect(names(Meta), "FKPM", negate = FALSE)=="FALSE")){

    if(all(rownames(Meta[[which(str_detect(names(Meta), "FKPM", negate = FALSE)=="TRUE")]])%in%x$gene_names)==F) {message("FKPM Matrice already with gene protein coding names as rownames")} else {

    message("Building FKPM.proteincoding matrice expression")

    y <- Meta[[which(str_detect(names(Meta), "FKPM", negate = FALSE)=="TRUE")]][rownames(x),]
    y$gene_names <- x$gene_names
    y <- y %>% group_by(gene_names) %>% summarise_each(funs(mean))
    y <- as.data.frame(y)
    rownames(y) <- y$gene_names
    y$gene_names <- NULL
    Meta$FKPM.proteincoding.matrix <- as.data.frame(y)
  }}else{message("No FKPM data matrice in Meta oject existing. \n Add a matrice with AddExpressionMatrix function")}
    }else{message("FKPM.proteincoding in Meta object already existing. \n Not computing code.")}




  return(Meta)
}



