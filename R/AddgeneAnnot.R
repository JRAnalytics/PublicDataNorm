#' AddgeneAnnot to Meta
#'
#' @param Meta  Meta object to add geneAnnotation
#' @param gtf.file.dir set xorking directory to path to gtf.file
#' @param gtf.files  gtf file to load
#' @import stringr
#' @import AnnotationDbi
#' @import data.table
#' @return return a Meta list of dataframe with geneAnntotation
#' @export
#'
#' @examples "none"
#'
AddgeneAnnot <- function(Meta ,gtf.file.dir, gtf.files){

  Meta <- Meta

  path <- file.path(gtf.file.dir)


  geneAnnot <- geneAnnotation(gtf.files = paste0(path,"/",gtf.files) ,saverds = F)

  zz <- which(str_detect(names(Meta), c("Raw", "matrix")))

  gene <- rownames(Meta[[zz]])

  if(all(str_detect(gene, "ENSG000"))){

    message("Data matrice row names are ENSEMBL gene names.")

  Meta[[zz]] <- Meta[[zz]][-which(duplicated(gene)),]
  gene <- unlist(lapply(str_split(rownames(Meta[[zz]]),"[.]"),"[[",1))
  rownames(Meta[[zz]]) <- gene

  geneAnnot <- geneAnnot[rownames(geneAnnot)%in%gene,]
  geneAnnot <- geneAnnot[order(geneAnnot$GeneSymbol,decreasing = F),]
  Meta$geneAnnotation <- geneAnnot
  } else {

    gene <- unlist(lapply(str_split(rownames(Meta[[zz]]),"[|]"),"[[",1))
  message("Data matrice row names are already in gene names.")

    geneAnnot <- geneAnnot[order(geneAnnot$GeneSymbol,decreasing = F),]



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
