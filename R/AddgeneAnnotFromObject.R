#' AddgeneAnnotFromObject to Metadata
#'
#' @param Metadata  Metadata object to add geneAnnotation
#' @param object object to add
#' @param force.replace set as F. T : replace an already object with the same name
#' @import stringr
#' @import AnnotationDbi
#' @import data.table
#' @return return a Metadata list of dataframe with geneAnntotation
#' @export
#'
#' @examples "none"
#'
AddgeneAnnotFromObject <- function(Metadata ,object, force.replace=F){

  Metadata <- Metadata

  if(!is.null(Metadata$geneAnnotation)){

    message("geneAnnotation already loaded.")

    if(force.replace==F){stop("Set force.replace==T to subset object.")}
    message("Subsetting object.")}


  if(is.null(object)){stop("Error in geneAnnotation function : no object found")}
  if(!all(class(object)%in%c("data.frame","matrix","array"))){stop("Object is not of a data.frame or a matrix class object.")}


  zz <- which(attributes(Metadata)$Data.Type=="Count")[1]

  if(length(zz)!=1){ stop("attributes(Metadata)$Data.Type=='Count' & attributes(Metadata)$Export=='Yes',line 32, more than 1 object are detected.")}


  gene <- rownames(Metadata[[zz]])
  geneAnnot = as.matrix(object)

  if(!summary(gene%in%geneAnnot)["TRUE"]>0){
    stop("No genes ('rownames(Metadata[[1]])') found in object.")} else {
      sel = which(geneAnnot %in% gene)
      col = which(apply(geneAnnot, 2, function(x) which(x %in% geneAnnot[sel[1]]))>0)[1]}


  if(all(str_detect(gene, "ENSG000"))){

    message("Data matrice row names are ENSEMBL gene names.")

    if(length(gene)!=length(unique(gene))) { Metadata[[zz]] <- Metadata[[zz]][-which(duplicated(gene)),] }
    if(length(gene)==length(unique(gene))){ message("No rownames of raw matrix are duplicated")}

    gene <- unlist(lapply(str_split(rownames(Metadata[[zz]]),"[.]"),"[[",1))
    rownames(Metadata[[zz]]) <- gene

    geneAnnot =as.data.frame(geneAnnot)

    geneAnnot <- geneAnnot[geneAnnot[,col]%in%gene,]

    if(nrow(geneAnnot)==0){ stop("if(nrow(geneAnnot)==0), line 43, N rows of geneAnnotation is 0")}

    geneAnnot <- geneAnnot[order(geneAnnot[,col],decreasing = F),]

    if(is.null(Metadata$geneAnnotation)){ attributes(Metadata)$Data.Type <-  c(attributes(Metadata)$Data.Type, "geneAnnot")
    attributes(Metadata)$Export <- c(attributes(Metadata)$Export,"Yes") }

    Metadata$geneAnnotation <- geneAnnot


  } else {

    gene <- unlist(lapply(str_split(rownames(Metadata[[zz]]),"[|]"),"[[",1))
    print(message("Data matrice row names are already in gene names."))

    geneAnnot =as.data.frame(geneAnnot)

    geneAnnot <- geneAnnot[order(geneAnnot[,col],decreasing = F),]

    if(is.null(Metadata$geneAnnotation)){ attributes(Metadata)$Data.Type <-  c(attributes(Metadata)$Data.Type, "geneAnnot")
    attributes(Metadata)$Export <- c(attributes(Metadata)$Export,"Yes") }

    Metadata$geneAnnotation <- geneAnnot



  }


  return(Metadata)
}
