#' CheckMeta : Checking samples accors MetaData files
#'
#' @param MetaData a MetaData list
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


  m <- which(attributes(MetaData)$Data.Type=="Count")

  g =  which(attributes(MetaData)$Data.Type=="geneAnnot")[1]
  if(!is.na(g)){geneAnnot = as.matrix(MetaData[[g]])}


  if(attributes(MetaData)$Omics.type!="Single.Cell"){
  c <- which(attributes(MetaData)$Data.Type=="Clinic" | attributes(MetaData)$Data.Type=="SamplesAnnot")}

  if(attributes(MetaData)$Omics.type=="Single.Cell"){
    c <- which(attributes(MetaData)$Data.Type=="Clinic" & attributes(MetaData)$Export=="No")}

  ColN <- colnames(MetaData[[m[1]]])

  message("-------------------------")
  message(paste("Checking colnames of matrices in MetaData from", names(MetaData)[m[1]]))
  message("-------------------------")
  for (i in m){



    if(all(ColN %in% colnames(MetaData[[i]]))==T) {   message(paste(MetaDataN[i]), " colnames : PASS") } else {  message(paste(MetaDataN[i]), " colnames : FAIL") }

  }

  if(!is.na(g)){
  message("-------------------------")
  message("Checking Matrices probes")
  message("-------------------------")


    for (i in m){

      gene = rownames(MetaData[[i]])
      sel = which(geneAnnot %in% gene)
      col = which(apply(geneAnnot, 2, function(x) which(x %in% geneAnnot[sel[1]]))>0)[1]


      if(all(str_detect(gene, "ENSG")==T)) { message(paste(names(MetaData[i]), " gene probes as ENSEMBL"))


        suma <- summary(gene%in%geneAnnot[,col])
        if(length(suma)==3) {  names(suma) <- c("Mode", "Gene not found", "Found")  }
        if(length(suma)==2 & names(suma)[2]=="TRUE") {  names(suma) <- c("Mode", "Found") }
        if(length(suma)==2 & names(suma)[2]=="FALSE") {  names(suma) <- c("Mode", "Gene not found") }
        print(suma)
      } else {

      if(all(is.numeric(as.numeric(gene))))  { message(paste(names(MetaData[i]), " gene probes as ENTREZ gene id"))

        suma <- summary(gene%in%geneAnnot[,col])
        if(length(suma)==3) {  names(suma) <- c("Mode", "Gene not found", "Found")  }
        if(length(suma)==2 & names(suma)[2]=="TRUE") {  names(suma) <- c("Mode", "Found") }
        if(length(suma)==2 & names(suma)[2]=="FALSE") {  names(suma) <- c("Mode", "Gene not found") }
        print(suma)
      } else {


       if(all(str_detect(gene, "ILMN_")==T)) { message(paste(names(MetaData[i]), " gene probes as Illumina Bead Array Probes"))

         suma <- summary(gene%in%geneAnnot[,col])
        if(length(suma)==3) {  names(suma) <- c("Mode", "Gene not found", "Found")  }
        if(length(suma)==2 & names(suma)[2]=="TRUE") {  names(suma) <- c("Mode", "Found") }
        if(length(suma)==2 & names(suma)[2]=="FALSE") {  names(suma) <- c("Mode", "Gene not found") }
        print(suma)
      } else{


      if(all(str_detect(gene, "_at")==T)) { message(paste(names(MetaData[i]), " gene probes as Illumina Microarray Probes"))

        suma <- summary(gene%in%geneAnnot[,col])
        if(length(suma)==3) {  names(suma) <- c("Mode", "Gene not found", "Found")  }
        if(length(suma)==2 & names(suma)[2]=="TRUE") {  names(suma) <- c("Mode", "Found") }
        if(length(suma)==2 & names(suma)[2]=="FALSE") {  names(suma) <- c("Mode", "Gene not found") }
        print(suma)
      } else {



      if(length(which(gene%in%geneAnnot$GeneSymbol))>1)   {

        message(paste(names(MetaData[i]), " gene probes as genes Symbols"))

        suma <- summary(gene%in%geneAnnot[,col])
        if(length(suma)==3) {  names(suma) <- c("Mode", "Gene not found", "Found")  }
        if(length(suma)==2 & names(suma)[2]=="TRUE") {  names(suma) <- c("Mode", "Found") }
        if(length(suma)==2 & names(suma)[2]=="FALSE") {  names(suma) <- c("Mode", "Gene not found") }
        print(suma)
      } else {

      if(all(gene%in%geneAnnot)) { message(paste(names(MetaData[i]), " gene probes manualy entered from published data."))

        suma <- summary(gene%in%geneAnnot[,col])
        if(length(suma)==3) {  names(suma) <- c("Mode", "Gene not found", "Found")  }
        if(length(suma)==2 & names(suma)[2]=="TRUE") {  names(suma) <- c("Mode", "Found") }
        if(length(suma)==2 & names(suma)[2]=="FALSE") {  names(suma) <- c("Mode", "Gene not found") }
        print(suma)
      }

    }}}}}}}




  message("-------------------------")
if(length(c>0)){
  message(paste("Checking colnames of", names(MetaData)[m[1]] ,"in MetaData in clinical data."))
  message("-------------------------")

  for (i in c){

    if(attributes(MetaData)$Omics.type!="Single.Cell"){

    if(all(ColN%in%as.matrix(MetaData[[i]]))) {   message(paste(MetaDataN[i]), " : PASS") } else {  message(paste(MetaDataN[i]), " : FAIL") }
    }


    if(attributes(MetaData)$Omics.type=="Single.Cell"){

      message(paste("Patients from Single.Cell data:", names(MetaData)[i]))
      tot=0
      for (z in gene) {
        t = summary(str_detect(pattern = z, ColN))["TRUE"][1]

        if(is.na(as.numeric(t))){ t = 0}

        tot=tot+as.numeric(t)


      message(c(z," N= ",as.numeric(t)))

      }
      message("Total = " , tot, "\nAre all Patients found in Expression matrix ? ", tot/length(ColN)==1)
      message("-------------------------")



     p =  which(attributes(MetaData)$Data.Type=="SamplesAnnot" & attributes(MetaData)$Export=="No")

     if(!is.null(p)){

      message(paste("Samples from Single.Cell data:", names(MetaData)[p]))

       tot = as.numeric(summary(rownames(MetaData[[p]])%in%ColN)["TRUE"][1])

       message("Total = " , tot, "\nAre all Patients found in Expression matrix ? ", tot/length(ColN)==1)
       message("-------------------------")
   }

    }
  }}









}
