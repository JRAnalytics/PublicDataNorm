#' AddObjetToMeta
#'
#' @param Metadata  a Metadata Object to fill
#' @param object an object to add
#' @param Omics.type one of these categories : "RNAseq", "Single.Cell", "Microarray", "Spatial".
#' @param name name to apply in Metadata objet list
#' @param geneFilter Default F, select only genes found in both the geneAnnot and the Expression Matrix.
#' @param SamplesFilter Default F, select only samples found in both the clinical data and the Expression Matrix.
#' @param Data.type attribute c("Count","Clinic", "SamplesAnnot","geneAnnot","CellAnnot" )
#' @param Export data to export after cleaning c("Yes", "No")
#' @import matrixStats
#'
#' @return a Metadata Object
#' @export
#'
#' @examples "none"
AddObjetToMeta <- function(Metadata,
                           object,
                           Omics.type = c("RNAseq", "Single.Cell", "Microarray", "Spatial") ,
                           name,
                           geneFilter = F,
                           SamplesFilter = F,
                           Data.type = c("Count","Clinic", "SamplesAnnot","geneAnnot","CellAnnot" ),
                           Export = c("Yes", "No")){


  if(!is.list(Metadata)) { stop("Metadata object should be a list")}
  if(is.null(object)) { stop("Object is null")}
  if(is.null(name)) { stop("name should be specified, as character")}
  if(!is.character(name)) { stop("name should be specified as character")}
  if(is.null(Data.type)) { stop("Data.type should be specified")}
  if(is.null(Export)) { stop("Export should be specified")}
  if(!Data.type%in%c("Count","Clinic", "SamplesAnnot","geneAnnot","CellAnnot" )) { stop("Data.type should be specified from this values c('Count','Clinic', 'SamplesAnnot','geneAnnot','CellAnnot' ) ")}
  if(!Export%in%c("Yes", "No")) { stop("Export should be specified 'Yes' or 'No'")}

  if(is.null(Omics.type)){ stop("Omics.type must be from one of this character c(RNAseq, Single.Cell, Microarray, Spatial)")}
  if(!inherits(Omics.type, "character")){ stop("Omics.type is not a character string")}
  if(!Omics.type%in% c("RNAseq", "Single.Cell", "Microarray", "Spatial") ){ stop("Omics.type must be from one of this character c(RNAseq, Single.Cell, Microarray, Spatial)")}



  l <- length(Metadata)



  if(Data.type=="geneAnnot") {

    zz <- which(attributes(Metadata)$Data.Type=="Count")[1]
    if(length(zz)==0){stop("No Count found in Metadata object")}
    gene <-  rownames(Metadata[[zz]])
    if(length(gene)==length(unique(gene))){ message("No rownames of Count matrix are duplicated")}
    val <- gene[which(gene%in%as.matrix(object))[1]]
    colT <- which(matrixStats::colAnys(as.matrix(object),value = val))

    if(length(colT)>1){ colT <- colT[1]}



      if(all(str_detect(gene, "ENSG")==T)) { message("Data matrice row names are as ENSEMBL.")
        if(length(colT)==0){stop("No genes as ENSEMBL found in geneAnnot.")}}

      if(all(str_detect(gene, "ILMN_")==T)) { message("Data matrice row names are as Illumina Bead Array Probes")
        if(length(colT)==0){stop("No genes as Illumina Bead Array Probes 'ILMN_' found in geneAnnot.")}}

      if(all(is.numeric(as.numeric(gene)))==T)  { message("Data matrice row names are as ENTREZ gene id")
        if(length(colT)==0){stop("No genes as ENTREZ.ID found in geneAnnot.")}}

      if(all(str_detect(gene, "_at")==T)) { message("Data matrice row names are as Illumina Microarray Probes")
        if(length(colT)==0){stop("No genes as ILLUMINA '_at' found in geneAnnot.")}}

      if("ACTB"%in%gene){ message("Data matrice row names are in GeneSymbols.")
        if(length(colT)==0){stop("No genes as GeneSymbols found in geneAnnot.")}}


      message("Found genes :")
      print(summary(gene%in%as.matrix(object)))

      message("Selecting retrieved genes from Count matrix")

      if(geneFilter==T){

        object <- object[object[,colT]%in%gene,]

      if(length(object[,colT])!=length(unique(object[,colT]))){object <- object[!duplicated(object[,colT]),]} }


      }

  if(Data.type=="Clinic"|Data.type== "SamplesAnnot"){

    zz <- which(attributes(Metadata)$Data.Type=="Count")[1]
    samples <- colnames(Metadata[[zz]])
    val <- samples[which(samples%in%as.matrix(object))[1]]
    colT <- which(matrixStats::colAnys(as.matrix(object),value = val))

    if(length(colT)>1){ colT <- colT[1]}

    message("Found Samples :")
    print(summary(samples%in%as.matrix(object)))
    message("Unfound Samples :")
    print(samples[!samples%in%as.matrix(object)])

    if(SamplesFilter==T){
      message("Selecting only Samples present in both Count and clinical.data.")

      # Diviser chaque élément de la colonne en un vecteur de sous-chaînes
      substrings <- strsplit(as.data.frame(object)[,colT], ";")
      # Vérifier si chaque valeur du vecteur est présente dans chaque vecteur de sous-chaînes
      est_present <- sapply(substrings, function(x) any(colnames(Metadata[[zz]]) %in% x))

      # Récupérer les index de position correspondants
      indices <- which(est_present)
      object <- object[indices,]


    }


    }


  Metadata[[l+1]] <- object
  names( Metadata )[l+1] <- name

  attributes(Metadata)$Data.Type[l+1] <- Data.type

  attributes(Metadata)$Export[l+1] <- Export
  attributes(Metadata)$Cleaned[l+1] = "No"


return(Metadata)

}
