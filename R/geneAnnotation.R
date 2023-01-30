#' geneAnnotation function to form a genes annotated data frame for data claening rnaseq data.
#'
#' @param gtf.files GTF file to load, a character if wd is already set, or an object of GTF file already load
#' @param file character to names the geneAnnotation rds file
#' @param saverds if TRUE save a ".rds" file with the file field
#' @import data.table
#' @import AnnotationDbi

#' @importFrom org.Hs.eg.db org.Hs.eg.db
#' @return return a data.frame with ensembl Id as rownames and attributes for each genes.
#'
#' @examples
#' "none"
#'
#' @export
geneAnnotation <- function(gtf.files = "gencode.v40.annotation.gtf", saverds=TRUE,file= ".rds" ){


#loading of the gtf files
  if(class(gtf.files)=="character") {



      genes <- fread(gtf.files)
      setnames(genes, names(genes), c("chr","source","type","start","end","score","strand","phase","attributes") )
  } else {

     if(!all(colnames(gtf.files)==c("chr","source","type","start","end","score","strand","phase","attributes")))
      stop("GTF file has not the rigth structure as c('chr','source','type','start','end','score','strand','phase','attributes') ")}



  split_gene <- as.data.frame(subset(genes, subset = c(type=="gene" | type== "transcript") ))

  if(nrow(split_gene)==0){ stop("Fun geneAnnotation : failed in subseting .gtf file")}

  split_gene$EnsemblID <- unlist(lapply(split_gene $attributes, extract_attributes, "gene_id"))
  split_gene$EnsemblID <- unlist(lapply(strsplit(split_gene$EnsemblID, "[.]"), "[[", 1))

  if(length(which(duplicated(split_gene $EnsemblID)))>0){split_gene  <- split_gene[-which(duplicated(split_gene$EnsemblID)),] }

  rownames(split_gene ) <- split_gene$EnsemblID

  split_gene$GeneSymbol <- unlist(lapply(split_gene $attributes, extract_attributes, "gene_name"))
  split_gene$GeneType <- unlist(lapply(split_gene $attributes, extract_attributes, "gene_type"))

  if(!all(is.na(unlist(lapply(split_gene $attributes, extract_attributes, "hgnc_id"))))) { split_gene $GeneType <- unlist(lapply(split_gene $attributes, extract_attributes, "gene_type"))}

  if(str_detect(gtf.files, "19")) {site <-"http://genome.ucsc.edu/cgi-bin/hgTracks?db=hg19&position=" }
  if(str_detect(gtf.files, "v33")) {site <-"http://genome.ucsc.edu/cgi-bin/hgTracks?db=hg33&position=" }
  if(str_detect(gtf.files, "v40")) {site <-"http://genome.ucsc.edu/cgi-bin/hgTracks?db=hg40&position=" }
  if(str_detect(gtf.files, "v41")) {site <-"http://genome.ucsc.edu/cgi-bin/hgTracks?db=hg41&position=" }
  if(str_detect(gtf.files, "v42")) {site <-"http://genome.ucsc.edu/cgi-bin/hgTracks?db=hg42&position=" }

  split_gene$UCSC_link <- paste0(site,split_gene$chr,"%3A",split_gene$start ,"%2D",split_gene$end)


  zz <- which(unlist(lapply(strsplit(split_gene$chr,"chr"), length))>1)
  split_gene$Ensemble_link <- rep(NA,nrow(split_gene))
  split_gene$Ensemble_link[zz] <- paste0("http://www.ensembl.org/Homo_sapiens/Gene/Summary?g=",split_gene$EnsemblID,";r=",lapply(strsplit(split_gene$chr,"chr")[zz],"[[",2),":",split_gene$start ,"-",split_gene$end)

  split_gene [which(is.na(split_gene$HGNC)),"HGNC"] <- ":-"
  split_gene [which((split_gene $HGNC)==":-"),"HGNC"] <- "-"

  split_gene$Entrez.id <- mapIds(org.Hs.eg.db, split_gene$GeneSymbol, 'ENTREZID', 'SYMBOL')
  split_gene$Entrez_geneID_link <- paste0("https://www.ncbi.nlm.nih.gov/gene/?term=",split_gene$Entrez.id)
  split_gene[which(is.na(split_gene$Entrez.id)),"Entrez_geneID_link"] <- "-"


if(nrow(split_gene)==0){ stop("Fun geneAnnotation return a DF with 0 rows")}

return(split_gene)

if(saverds==TRUE) {
saveRDS(split_gene,file)}

}
