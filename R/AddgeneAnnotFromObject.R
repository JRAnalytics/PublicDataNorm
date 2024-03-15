#' AddgeneAnnotFromObject to Metadata
#'
#' @param Metadata  Metadata object to add geneAnnotation
#' @param object object to add
#' @param geneAnnotIDcolumn column object to fetch geneIDs
#' @param force.replace set as F. T : replace an already object with the same name
#' @param Filter.Genes default F, if T, keep only retrieved genes in Count matrix
#' @import stringr
#' @import AnnotationDbi
#' @import data.table
#' @return return a Metadata list of dataframe with geneAnntotation
#' @export
#'
#' @examples "none"
#'
AddgeneAnnotFromObject <- function(Metadata=NULL ,object=NULL, geneAnnotIDcolumn= NULL, Filter.Genes= F, force.replace=F){


  if(is.null(Metadata)){stop("A Metadata object is needed.")}
  if(!is.null(Metadata$geneAnnotation)){

    message("geneAnnotation already loaded.")

    if(force.replace==F){stop("Set force.replace==T to subset object.")}
    message("Subsetting object.")}


  if(is.null(object)){stop("Error in geneAnnotation function : no object found")}
  if(!all(class(object)%in%c("data.frame","matrix","array"))){stop("Object is not of a data.frame or a matrix class object.")}

  if(is.null(geneAnnotIDcolumn)){stop("geneAnnotIDcolumn must be refered to a colname from object to add.")}


  zz <- which(attributes(Metadata)$Data.Type=="Count")[1]
  if(length(zz)==0){stop("No Count found in Metadata object")}
  gene <-  rownames(Metadata[[zz]])

  if(length(gene)==length(unique(gene))){ message("No rownames of Count matrix are duplicated")}


  colT <- geneAnnotIDcolumn


  if(all(str_detect(gene, "ENSG")==T)) { message("Data matrice row names are as ENSEMBL.")
    }

  if(all(str_detect(gene, "ILMN_")==T)) { message("Data matrice row names are as Illumina Bead Array Probes")
   }

  if(is.numeric(na.omit(as.numeric(gene))) & length(na.omit(as.numeric(gene)))!=0)  { message("Data matrice row names are as ENTREZ gene id")
  }

  if(all(str_detect(gene, "_at")==T)) { message("Data matrice row names are as Illumina Microarray Probes")
    }

  if("ACTB"%in%gene){ message("Data matrice row names are in GeneSymbols.")
   }


  message("Found genes :")
  print(summary(gene%in%object[,colT]))

  message("Selecting retrieved genes from Count matrix")

  if(Filter.Genes==T){

    object <- object[object[,colT]%in%gene,]

    if(length(object[,colT])!=length(unique(object[,colT]))){object <- object[!duplicated(object[,colT]),]} }

    if(is.null(Metadata$geneAnnotation)){ attributes(Metadata)$Data.Type <-  c(attributes(Metadata)$Data.Type, "geneAnnot")
    attributes(Metadata)$Export <- c(attributes(Metadata)$Export,"Yes")
    attributes(Metadata)$Cleaned = c(attributes(Metadata)$Cleaned,"Yes")} else {

      }

  Metadata$geneAnnotation = object



  return(Metadata)
}
