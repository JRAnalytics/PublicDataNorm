#' AddgeneAnnot to Meta
#'
#' @param Meta  Meta object to add geneAnnotation
#' @param gtf.file.dir set xorking directory to path to gtf.file
#' @param gtf.files  gtf file to load
#' @param force.replace set as F. T : replace an already object with the same name
#' @import stringr
#' @import AnnotationDbi
#' @import data.table
#' @return return a Meta list of dataframe with geneAnntotation
#' @export
#'
#' @examples "none"
#'
AddgeneAnnot <- function(Meta ,gtf.file.dir, gtf.files, force.replace=F){

  Meta <- Meta

  path <- file.path(gtf.file.dir)

  if(!is.null(Meta$geneAnnotation)){

    message("geneAnnotation already loaded.")

    if(force.replace==F){stop("Set force.replace==T to subset object.")}
    message("Subsetting object.")}

  if(!str_detect(gtf.files, ".gtf")){stop("A '.gtf' file is needed for AddgeneAnnot function.\n If an R object already exist with annotations, use AddObjetToMeta() function.") }
  if(!file.exists(paste0(path,"/",gtf.files))){ stop("File does not exist in directories.")}


  geneAnnot <- geneAnnotation(gtf.files = paste0(path,"/",gtf.files) ,saverds = F)


  if(is.null(geneAnnot)){stop("Error in geneAnnotation function : is.null(geneAnnot==TRUE")}

  if(is.null(rownames(geneAnnot))){stop("Error in geneAnnotation function : no rownames in geneAnnot. But is not Null")}

  zz <- which(attributes(Meta)$Data.Type=="Expression.Matrix")[1]

  if(length(zz)!=1){ stop("attributes(Meta)$Data.Type=='Expression.Matrix' & attributes(Meta)$Raw.data=='Yes',line 27, more than 1 object are detected.")}

  ###rajouté ecrasement

   gene <- rownames(Meta[[zz]])

  if(all(str_detect(gene, "ENSG000"))){

    print(message("Data matrice row names are ENSEMBL gene names."))

    if(length(gene)!=length(unique(gene))) { Meta[[zz]] <- Meta[[zz]][-which(duplicated(gene)),] }
    if(length(gene)==length(unique(gene))){ message("No rownames of raw matrix are duplicated")}

  gene <- unlist(lapply(str_split(rownames(Meta[[zz]]),"[.]"),"[[",1))
  rownames(Meta[[zz]]) <- gene

  geneAnnot <- geneAnnot[rownames(geneAnnot)%in%gene,]

  if(nrow(geneAnnot)==0){ stop("if(nrow(geneAnnot)==0), line 43, N rows of geneAnnotation is 0")}

  geneAnnot <- geneAnnot[order(geneAnnot$GeneSymbol,decreasing = F),]

  if(is.null(Meta$geneAnnotation)){ attributes(Meta)$Data.Type <-  c(attributes(Meta)$Data.Type, "geneAnnotation.file")
  attributes(Meta)$Raw.data <- c(attributes(Meta)$Raw.data,"No") }

  Meta$geneAnnotation <- geneAnnot


  } else {

    gene <- unlist(lapply(str_split(rownames(Meta[[zz]]),"[|]"),"[[",1))
  print(message("Data matrice row names are already in gene names."))

    geneAnnot <- geneAnnot[order(geneAnnot$GeneSymbol,decreasing = F),]

    if(is.null(Meta$geneAnnotation)){ attributes(Meta)$Data.Type <-  c(attributes(Meta)$Data.Type, "geneAnnotation.file")
    attributes(Meta)$Raw.data <- c(attributes(Meta)$Raw.data,"No") }

    Meta$geneAnnotation <- geneAnnot



  }

  if(file.exists("Readme.txt")){

    tme <- Sys.Date()
    tme <- format(tme, format="%B %d %Y")

    name <- "geneAnnotation"

    sp <- data.frame("Type"="---------:" ,"Description"="---------")
    mod <- data.frame("Type"="geneAnnotation added the:" ,"Description"=tme)

    dt <- rbind(sp,mod,sp)


     if(str_detect(name, c("geneAnnotation"))==T){
      nr <- nrow(Meta[[name]])
      nc <- ncol(Meta[[name]])

      if(str_detect(name, "Annotation")){ Assay ="Gene annotation "}

      if(gtf.files =="gencode.v19.chr_patch_hapl_scaff.annotation.gtf") {Version ="Genome of reference hg19 ; GRCh37.p13 ; https://www.gencodegenes.org/human/release_19.html" }
      if(gtf.files =="gencode.v40.annotation.gtf") {Version ="Genome of reference hg38; V40 ; GRCh37 ; https://www.gencodegenes.org/human/release_40lift37.html" }
      if(gtf.files =="gencode.v41.chr_patch_hapl_scaff.annotation.gtf") {Version ="Genome of reference hg38 ; V41 ; GRCh37 ; https://www.gencodegenes.org/human/release_40lift37.html" }
      if(gtf.files =="gencode.v33.annotation.gtf") {Version ="Genome of reference hg38, V33 ; GRCh38.p13  ;https://www.gencodegenes.org/human/release_33.html" }



      ltest <- data.frame("Type"=c("File: " , "Class: " ,"Dimension: ", "Assay: ","Version:" ,"Rownames: ","Colnames: "),"Description" = c(
        paste0(name,".csv"),
        class(Meta[[name]]),
        paste(nr,"x",nc),
        Assay,
        Version,
        paste(rownames(Meta[[name]])[1],"...",rownames(Meta[[name]])[nrow(Meta[[name]])]),
        paste(colnames(Meta[[name]])[2],"...",colnames(Meta[[name]])[ncol(Meta[[name]])])
      ))



      dt <- rbind(dt, ltest,sp)

      write.table(dt,"Readme.txt",row.names = F,quote = FALSE, append = T, col.names=FALSE)
      file.show("Readme.txt")
      closeAllConnections()}}



  return(Meta)
}
