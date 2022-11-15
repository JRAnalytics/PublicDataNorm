#' SurvivalMeta Sruvial OS or PFS analysis from MetaData
#'
#' @param Metadata a metaData object
#' @param target gene ID
#' @param Group.sep median ou quantil 0.25/0.75 separation
#' @param type OS of PFS
#' @import survival
#' @importFrom stats median
#' @importFrom stats quantile
#' @import ggplot2
#' @import survminer
#'
#' @return a survival plot and gene expression plus grouping in Patient's clinical data sheet.
#' @export
#'
#' @examples "non"
SurvivalMeta <- function(Metadata, target, Group.sep = c("median", "Q1 et 3"), type){



  clinic <- Metadata[[which(str_detect(names(Metadata),"Patient"))]]
  coln <- colnames(Metadata[[which(str_detect(names(Metadata),"Patient"))]])
  patient <- rownames((Metadata[[which(str_detect(names(Metadata),"Patient"))]]))

  if(is.null(which(str_detect(names(Metadata),"Normalized")))) { stop("Need a Normalized expression data set for Survival analyses.")}
  DF <-  as.data.frame(Metadata[[which(str_detect(names(Metadata),"Normalized"))]])

  if(all(is.na(DF[target,]))){stop("Gene not retrived in Normalized dataset")}

  clinic$gene <- t(DF[target,rownames(clinic)])[,1]

  if(Group.sep == c("median")){

    clinic$Group <- ifelse(clinic$gene>median(clinic$gene), "High", "Low")
  }
  if(Group.sep == c("Q1 et 3")){

    clinic$Group <- ifelse(clinic$gene>quantile(clinic$gene,0.75), "High",ifelse(clinic$gene<quantile(clinic$gene,0.25), "Low", NA))

    }

if(type == "OS"){

  OS01 <- which(str_detect(colnames(clinic), "OS01"))
  OStime <- which(str_detect(colnames(clinic), "OS") & colnames(clinic)!="OS01" )

} else {  OS01 <- which(str_detect(colnames(clinic), "PFS01"))
OStime <- which(str_detect(colnames(clinic), "PFS") & colnames(clinic)!="PFS01" ) }

  if(length(which(str_detect(unique(clinic[,OS01]), "Alive")))>1) {

    clinic[,OS01] <-  gsub(clinic[,OS01], pattern = "Alive",replacement =  0)
   clinic[,OS01] <-  gsub(clinic[,OS01], pattern = "Dead",replacement =  1) ; clinic[,OS01] <- as.numeric(clinic[,OS01])   }

  surv_object <- Surv(time = (clinic[,OStime]), event = (clinic[,OS01]))
  fit1 <- survfit(surv_object ~ Group, data = clinic )
  print(ggsurvplot(fit1, data = clinic, pval = TRUE, legend.title = paste(target,"median genes epxression"),
             legend.labs = c("High", "Low"), legend = "bottom"))

  if(length(which(str_detect(colnames(Metadata[[which(str_detect(names(Metadata),"Patient"))]]),target)))<1){

    Metadata[[which(str_detect(names(Metadata),"Patient"))]] <-  clinic

    colnames(Metadata[[which(str_detect(names(Metadata),"Patient"))]]) <-  c(coln,target ,paste(target, "_Group"))

    }



  return(Metadata)

}
