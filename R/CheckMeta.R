#' CheckMeta : Checking samples accors MetaData files
#'
#' @param MetaData a mMetaData list
#'
#' @return none
#' @export
#'
#' @examples "none"
CheckMeta <- function(MetaData) {

  if(is.null(MetaData)){stop("Need a MetaData List file")}
  if(!is.list(MetaData)){stop("Need a MetaData List file")}


  l <-length(names(MetaData))
  MetaDataN <- names(MetaData)


  m <- which(attributes(MetaData)$Data.Type=="Expression.Matrix")
  c <- which(attributes(MetaData)$Data.Type=="Patient.Clinical.data" |attributes(MetaData)$Data.Type=="Samples.Clinical.data")

  ColN <- colnames(MetaData[[m[1]]])

  message("-------------------------")
  message(paste("Checking colnames of matrices in MetaData from", names(MetaData)[m[1]]))
  message("-------------------------")
  for (i in m){



    if(all(ColN %in% colnames(MetaData[[i]]))==T) {   message(paste(MetaDataN[i]), " colnames : PASS") } else {  message(paste(MetaDataN[i]), " colnames : FAIL") }

  }
  message("-------------------------")
  message("Checking Matrices probes")
  message("-------------------------")


    for (i in m){



      if(all(str_detect(rownames(MetaData[[i]]), "ENSG")==T)) { message(paste(names(MetaData[i]), " gene probes as ENSEMBL"))

        suma <- summary(rownames(MetaData[[i]])%in%MetaData$geneAnnotation$EnsemblID)
        if(length(suma)==3) {  names(suma) <- c("Mode", "Gene not found", "Found")  }
        if(length(suma)==2 & names(suma)[2]=="TRUE") {  names(suma) <- c("Mode", "Found") }
        if(length(suma)==2 & names(suma)[2]=="FALSE") {  names(suma) <- c("Mode", "Gene not found") }
        print(suma)
      } else {

      if(all(rownames(MetaData[[i]])%in%MetaData$geneAnnotation$Entrez.id)==T)  { message(paste(names(MetaData[i]), " gene probes as ENTREZ gene id"))

        suma <- summary(rownames(MetaData[[i]])%in%MetaData$geneAnnotation$Entrez.id)
        if(length(suma)==3) {  names(suma) <- c("Mode", "Gene not found", "Found")  }
        if(length(suma)==2 & names(suma)[2]=="TRUE") {  names(suma) <- c("Mode", "Found") }
        if(length(suma)==2 & names(suma)[2]=="FALSE") {  names(suma) <- c("Mode", "Gene not found") }
        print(suma)
      } else {


       if(all(str_detect(rownames(MetaData[[i]]), "ILMN_")==T)) { message(paste(names(MetaData[i]), " gene probes as Illumina Bead Array Probes"))

        suma <- summary(rownames(MetaData[[i]])%in%rownames(MetaData$geneAnnotation))
        if(length(suma)==3) {  names(suma) <- c("Mode", "Gene not found", "Found")  }
        if(length(suma)==2 & names(suma)[2]=="TRUE") {  names(suma) <- c("Mode", "Found") }
        if(length(suma)==2 & names(suma)[2]=="FALSE") {  names(suma) <- c("Mode", "Gene not found") }
        print(suma)
      } else{


      if(all(str_detect(rownames(MetaData[[i]]), "_at")==T)) { message(paste(names(MetaData[i]), " gene probes as Illumina Microarray Probes"))

        suma <- summary(rownames(MetaData[[i]])%in%rownames(MetaData$geneAnnotation))
        if(length(suma)==3) {  names(suma) <- c("Mode", "Gene not found", "Found")  }
        if(length(suma)==2 & names(suma)[2]=="TRUE") {  names(suma) <- c("Mode", "Found") }
        if(length(suma)==2 & names(suma)[2]=="FALSE") {  names(suma) <- c("Mode", "Gene not found") }
        print(suma)
      } else {



      if(length(which(rownames(MetaData[[i]])%in%MetaData$geneAnnotation$GeneSymbol))>1)   {

        message(paste(names(MetaData[i]), " gene probes as genes Symbols"))

        suma <- summary(rownames(MetaData[[i]])%in%MetaData$geneAnnotation$GeneSymbol)
        if(length(suma)==3) {  names(suma) <- c("Mode", "Gene not found", "Found")  }
        if(length(suma)==2 & names(suma)[2]=="TRUE") {  names(suma) <- c("Mode", "Found") }
        if(length(suma)==2 & names(suma)[2]=="FALSE") {  names(suma) <- c("Mode", "Gene not found") }
        print(suma)
      } else {

      if(all(rownames(MetaData[[i]])%in%rownames(MetaData$geneAnnotation))) { message(paste(names(MetaData[i]), " gene probes manualy entered from published data."))

        suma <- summary(rownames(MetaData[[i]])%in%rownames(MetaData$geneAnnotation))
        if(length(suma)==3) {  names(suma) <- c("Mode", "Gene not found", "Found")  }
        if(length(suma)==2 & names(suma)[2]=="TRUE") {  names(suma) <- c("Mode", "Found") }
        if(length(suma)==2 & names(suma)[2]=="FALSE") {  names(suma) <- c("Mode", "Gene not found") }
        print(suma)
      }

    }}}}}}




  message("-------------------------")

  message(paste("Checking colnames of matrices in MetaData from", names(MetaData)[m[1]]), "in clinical data.")
  message("-------------------------")

  for (i in c){

    if(all(ColN%in%as.matrix(MetaData[[i]]))) {   message(paste(MetaDataN[i]), " colnames : PASS") } else {  message(paste(MetaDataN[i]), " rownames : FAIL") }

  }









}
