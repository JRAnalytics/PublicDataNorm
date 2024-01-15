#' AddgeneAnnotGtf to Metadata
#'
#' @param Metadata  Metadata object to add geneAnnotation
#' @param gtf.file.dir set xorking directory to path to gtf.file
#' @param gtf.files  gtf file to load
#' @param force.replace set as F. T : replace an already object with the same name
#' @import stringr
#' @import AnnotationDbi
#' @import data.table
#' @return return a Metadata list of dataframe with geneAnntotation
#' @export
#'
#' @examples "none"
#'
AddgeneAnnotGtf <- function(Metadata ,gtf.file.dir, gtf.files, force.replace=F){

  Metadata <- Metadata

  path <- file.path(gtf.file.dir)

  if(!is.null(Metadata$geneAnnotation)){

    message("geneAnnotation already loaded.")

    if(force.replace==F){stop("Set force.replace==T to replace sub-object.")}
    message("Replacing sub-object.")}

  if(!str_detect(gtf.files, ".gtf")){stop("A '.gtf' file is needed for AddgeneAnnot function.\n If an R object already exist with annotations, use AddObjetToMeta() function.") }
  if(!file.exists(paste0(path,"/",gtf.files))){ stop("File does not exist in directories.")}


  geneAnnot <- suppressMessages({.geneAnnotation(gtf.files = paste0(path,"/",gtf.files) ,saverds = F)})


  if(is.null(geneAnnot)){stop("Error in geneAnnotation function : is.null(geneAnnot==TRUE")}

  if(is.null(rownames(geneAnnot))){stop("Error in geneAnnotation function : no rownames in geneAnnot. But is not Null")}

  zz <- which(attributes(Metadata)$Data.Type=="Count")[1]

  if(length(zz)!=1){ stop("attributes(Metadata)$Data.Type=='Count' & attributes(Metadata)$Export=='Yes',line 27, more than 1 object are detected.")}

  ###rajouté ecrasement

  gene <- rownames(Metadata[[zz]])

  geneAnnot = as.matrix(geneAnnot)

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

    geneAnnot = as.data.frame(geneAnnot)
    geneAnnot <- geneAnnot[geneAnnot[,col]%in%gene,]

    if(nrow(geneAnnot)==0){ stop("if(nrow(geneAnnot)==0), line 43, N rows of geneAnnotation is 0")}

    geneAnnot <- geneAnnot[order(geneAnnot$GeneSymbol,decreasing = F),]

    if(is.null(Metadata$geneAnnotation)){ attributes(Metadata)$Data.Type <-  c(attributes(Metadata)$Data.Type, "geneAnnot")
    attributes(Metadata)$Export <- c(attributes(Metadata)$Export,"No") }

    Metadata$geneAnnotation <- geneAnnot


  } else {

    gene <- unlist(lapply(str_split(rownames(Metadata[[zz]]),"[|]"),"[[",1))
    print(message("Data matrice row names are already in gene names."))


    geneAnnot = as.data.frame(geneAnnot)
    geneAnnot <- geneAnnot[geneAnnot[,col]%in%gene,]

    geneAnnot <- geneAnnot[order(geneAnnot$GeneSymbol,decreasing = F),]

    if(is.null(Metadata$geneAnnotation)){ attributes(Metadata)$Data.Type <-  c(attributes(Metadata)$Data.Type, "geneAnnot")
    attributes(Metadata)$Export <- c(attributes(Metadata)$Export,"Yes") }

    Metadata$geneAnnotation <- geneAnnot



  }



  return(Metadata)
}
