#' AddObjetToMeta
#'
#' @param Meta  a Meta Object to fill
#' @param object an object to add
#' @param name name to apply in Meta objet list
#' @param geneFilter Default F, select only genes found in both the geneAnnotation.file and the Expression Matrix.
#' @param SamplesFilter Default F, select only samples found in both the clinical data and the Expression Matrix.
#' @param Data.type attribute c("Expression.Matrix","Patient.Clinical.data", "Samples.Clinical.data","geneAnnotation.file" )
#' @param Raw attribute c("Yes", "No")
#' @import matrixStats
#'
#' @return a Meta Object
#' @export
#'
#' @examples "none
AddObjetToMeta <- function(Meta, object, name ,
                           geneFilter = F,
                           SamplesFilter = F,
                           Data.type = c("Expression.Matrix","Patient.Clinical.data", "Samples.Clinical.data","geneAnnotation.file" ),
                           Raw = c("Yes", "No")){


  if(!is.list(Meta)) { stop("Meta object should be a list")}
  if(is.null(object)) { stop("Object is null")}
  if(is.null(name)) { stop("name should be specified, as character")}
  if(!is.character(name)) { stop("name should be specified as character")}
  if(is.null(Data.type)) { stop("Data.type should be specified")}
  if(is.null(Raw)) { stop("Raw should be specified")}
  if(!Data.type%in%c("Expression.Matrix","Patient.Clinical.data", "Samples.Clinical.data","geneAnnotation.file" )) { stop("Data.type should be specified from this values c('Expression.Matrix','Patient.Clinical.data', 'Samples.Clinical.data','geneAnnotation.file' ")}
  if(!Raw%in%c("Yes", "No")) { stop("Raw should be specified 'Yes' or 'No'")}


  l <- length(Meta)



  if(Data.type=="geneAnnotation.file") {

    zz <- which(attributes(Meta)$Data.Type=="Expression.Matrix")[1]
    if(length(zz)==0){stop("No Expression.Matrix found in Meta object")}
    gene <-  rownames(Meta[[zz]])
    if(length(gene)==length(unique(gene))){ message("No rownames of raw matrix are duplicated")}
    val <- gene[which(gene%in%as.matrix(object))[1]]
    colT <- which(matrixStats::colAnys(as.matrix(object),value = val))

    if(length(colT)>1){ colT <- colT[1]}



      if(all(str_detect(gene, "ENSG")==T)) { message("Data matrice row names are as ENSEMBL.")
        if(length(colT)==0){stop("No genes as ENSEMBL found in geneAnnotation.file.")}}

      if(all(str_detect(gene, "ILMN_")==T)) { message("Data matrice row names are as Illumina Bead Array Probes")
        if(length(colT)==0){stop("No genes as Illumina Bead Array Probes 'ILMN_' found in geneAnnotation.file.")}}

      if(all(is.numeric(as.numeric(gene)))==T)  { message("Data matrice row names are as ENTREZ gene id")
        if(length(colT)==0){stop("No genes as ENTREZ.ID found in geneAnnotation.file.")}}

      if(all(str_detect(gene, "_at")==T)) { message("Data matrice row names are as Illumina Microarray Probes")
        if(length(colT)==0){stop("No genes as ILLUMINA '_at' found in geneAnnotation.file.")}}

      if("ACTB"%in%gene){ message("Data matrice row names are in GeneSymbols.")
        if(length(colT)==0){stop("No genes as GeneSymbols found in geneAnnotation.file.")}}


      message("Found genes :")
      print(summary(gene%in%as.matrix(object)))

      message("Selecting retrieved genes from Raw expression matrix")

      if(geneFilter==T){

        object <- object[object[,colT]%in%gene,]

      if(length(object[,colT])!=length(unique(object[,colT]))){object <- object[!duplicated(object[,colT]),]} }


      }

  if(Data.type=="Patient.Clinical.data"|Data.type== "Samples.Clinical.data"){

    zz <- which(attributes(Meta)$Data.Type=="Expression.Matrix")[1]
    samples <- colnames(Meta[[zz]])
    val <- samples[which(samples%in%as.matrix(object))[1]]
    colT <- which(matrixStats::colAnys(as.matrix(object),value = val))

    if(length(colT)>1){ colT <- colT[1]}

    message("Found Samples :")
    print(summary(samples%in%as.matrix(object)))
    message("Unfound Samples :")
    print(samples[!samples%in%as.matrix(object)])

    if(SamplesFilter==T){
      message("Selecting only Samples present in both Expression.matrix and clinical.data.")
      object <- object[object[,colT]%in%samples,]
      Meta[[zz]] <- Meta[[zz]][samples%in%object[,colT],]


    }}


  Meta[[l+1]] <- object
  names( Meta )[l+1] <- name

  attributes(Meta)$Data.Type[l+1] <- Data.type

  attributes(Meta)$Raw.data[l+1] <- Raw



return(Meta)

}
