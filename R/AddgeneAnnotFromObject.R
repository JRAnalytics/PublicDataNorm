#' AddgeneAnnotFromObject to Meta
#'
#' @param Meta  Meta object to add geneAnnotation
#' @param object object to add
#' @param force.replace set as F. T : replace an already object with the same name
#' @import stringr
#' @import AnnotationDbi
#' @import data.table
#' @return return a Meta list of dataframe with geneAnntotation
#' @export
#'
#' @examples "none"
#'
AddgeneAnnotFromObject <- function(Meta ,object, force.replace=F){

  Meta <- Meta

  if(!is.null(Meta$geneAnnotation)){

    message("geneAnnotation already loaded.")

    if(force.replace==F){stop("Set force.replace==T to subset object.")}
    message("Subsetting object.")}


  if(is.null(object)){stop("Error in geneAnnotation function : no object found")}
  if(!all(class(object)%in%c("data.frame","matrix","array"))){stop("Object is not of a data.frame or a matrix class object.")}


  zz <- which(attributes(Meta)$Data.Type=="Count")[1]

  if(length(zz)!=1){ stop("attributes(Meta)$Data.Type=='Count' & attributes(Meta)$Export=='Yes',line 32, more than 1 object are detected.")}


  gene <- rownames(Meta[[zz]])
  geneAnnot = as.matrix(object)

  if(!summary(gene%in%geneAnnot)["TRUE"]>0){
    stop("No genes ('rownames(Meta[[1]])') found in object.")} else {
      sel = which(geneAnnot %in% gene)
      col = which(apply(geneAnnot, 2, function(x) which(x %in% geneAnnot[sel[1]]))>0)[1]}


  if(all(str_detect(gene, "ENSG000"))){

    message("Data matrice row names are ENSEMBL gene names.")

    if(length(gene)!=length(unique(gene))) { Meta[[zz]] <- Meta[[zz]][-which(duplicated(gene)),] }
    if(length(gene)==length(unique(gene))){ message("No rownames of raw matrix are duplicated")}

    gene <- unlist(lapply(str_split(rownames(Meta[[zz]]),"[.]"),"[[",1))
    rownames(Meta[[zz]]) <- gene

    geneAnnot =as.data.frame(geneAnnot)

    geneAnnot <- geneAnnot[geneAnnot[,col]%in%gene,]

    if(nrow(geneAnnot)==0){ stop("if(nrow(geneAnnot)==0), line 43, N rows of geneAnnotation is 0")}

    geneAnnot <- geneAnnot[order(geneAnnot[,col],decreasing = F),]

    if(is.null(Meta$geneAnnotation)){ attributes(Meta)$Data.Type <-  c(attributes(Meta)$Data.Type, "geneAnnot")
    attributes(Meta)$Export <- c(attributes(Meta)$Export,"Yes") }

    Meta$geneAnnotation <- geneAnnot


  } else {

    gene <- unlist(lapply(str_split(rownames(Meta[[zz]]),"[|]"),"[[",1))
    print(message("Data matrice row names are already in gene names."))

    geneAnnot =as.data.frame(geneAnnot)

    geneAnnot <- geneAnnot[order(geneAnnot[,col],decreasing = F),]

    if(is.null(Meta$geneAnnotation)){ attributes(Meta)$Data.Type <-  c(attributes(Meta)$Data.Type, "geneAnnot")
    attributes(Meta)$Export <- c(attributes(Meta)$Export,"Yes") }

    Meta$geneAnnotation <- geneAnnot



  }


  return(Meta)
}
