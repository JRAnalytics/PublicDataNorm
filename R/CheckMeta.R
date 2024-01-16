#' CheckMeta : Checking samples accors Metadata files
#'
#' @param Metadata a Metadata list
#'
#' @return none
#' @export
#'
#' @examples "none"
CheckMeta <- function(Metadata) {

  if(is.null(Metadata)){stop("Need a Metadata List file")}
  if(!is.list(Metadata)){stop("Need a Metadata List file")}


  l <-length(names(Metadata))
  MetaDataN <- names(Metadata)


  m <- which(attributes(Metadata)$Data.Type=="Count")

  g =  which(attributes(Metadata)$Data.Type=="geneAnnot")[1]
  if(!is.na(g)){geneAnnot = as.matrix(Metadata[[g]])}


  if(attributes(Metadata)$Omics.type!="Single.Cell"){
  c <- which(attributes(Metadata)$Data.Type=="Clinic" | attributes(Metadata)$Data.Type=="SamplesAnnot")}

  if(attributes(Metadata)$Omics.type=="Single.Cell"){
    c <- which(attributes(Metadata)$Data.Type=="Clinic"  & attributes(Metadata)$Export=="No")}

  ColN <- colnames(Metadata[[m[1]]])

  message("-------------------------")
  message(paste("Checking colnames of matrices in Metadata from", names(Metadata)[m[1]]))
  message("-------------------------")
  for (i in m){



    if(all(ColN %in% colnames(Metadata[[i]]))==T) {   message(paste(MetaDataN[i]), " colnames : PASS") } else {  message(paste(MetaDataN[i]), " colnames : FAIL") }

  }

  if(!is.na(g)){
  message("-------------------------")
  message("Checking Matrices probes")
  message("-------------------------")


    for (i in m){

      gene = rownames(Metadata[[i]])
      sel = which(geneAnnot %in% gene)
      col = which(apply(geneAnnot, 2, function(x) which(x %in% geneAnnot[sel[1]]))>0)[1]


      if(all(str_detect(gene, "ENSG")==T)) { message(paste(names(Metadata[i]), " gene probes as ENSEMBL"))


        suma <- summary(gene%in%geneAnnot[,col])
        if(length(suma)==3) {  names(suma) <- c("Mode", "Gene not found", "Found")  }
        if(length(suma)==2 & names(suma)[2]=="TRUE") {  names(suma) <- c("Mode", "Found") }
        if(length(suma)==2 & names(suma)[2]=="FALSE") {  names(suma) <- c("Mode", "Gene not found") }
        print(suma)
      } else {

      if(all(is.numeric(as.numeric(gene))))  { message(paste(names(Metadata[i]), " gene probes as ENTREZ gene id"))

        suma <- summary(gene%in%geneAnnot[,col])
        if(length(suma)==3) {  names(suma) <- c("Mode", "Gene not found", "Found")  }
        if(length(suma)==2 & names(suma)[2]=="TRUE") {  names(suma) <- c("Mode", "Found") }
        if(length(suma)==2 & names(suma)[2]=="FALSE") {  names(suma) <- c("Mode", "Gene not found") }
        print(suma)
      } else {


       if(all(str_detect(gene, "ILMN_")==T)) { message(paste(names(Metadata[i]), " gene probes as Illumina Bead Array Probes"))

         suma <- summary(gene%in%geneAnnot[,col])
        if(length(suma)==3) {  names(suma) <- c("Mode", "Gene not found", "Found")  }
        if(length(suma)==2 & names(suma)[2]=="TRUE") {  names(suma) <- c("Mode", "Found") }
        if(length(suma)==2 & names(suma)[2]=="FALSE") {  names(suma) <- c("Mode", "Gene not found") }
        print(suma)
      } else{


      if(all(str_detect(gene, "_at")==T)) { message(paste(names(Metadata[i]), " gene probes as Illumina Microarray Probes"))

        suma <- summary(gene%in%geneAnnot[,col])
        if(length(suma)==3) {  names(suma) <- c("Mode", "Gene not found", "Found")  }
        if(length(suma)==2 & names(suma)[2]=="TRUE") {  names(suma) <- c("Mode", "Found") }
        if(length(suma)==2 & names(suma)[2]=="FALSE") {  names(suma) <- c("Mode", "Gene not found") }
        print(suma)
      } else {



      if(length(which(gene%in%geneAnnot$GeneSymbol))>1)   {

        message(paste(names(Metadata[i]), " gene probes as genes Symbols"))

        suma <- summary(gene%in%geneAnnot[,col])
        if(length(suma)==3) {  names(suma) <- c("Mode", "Gene not found", "Found")  }
        if(length(suma)==2 & names(suma)[2]=="TRUE") {  names(suma) <- c("Mode", "Found") }
        if(length(suma)==2 & names(suma)[2]=="FALSE") {  names(suma) <- c("Mode", "Gene not found") }
        print(suma)
      } else {

      if(all(gene%in%geneAnnot)) { message(paste(names(Metadata[i]), " gene probes manualy entered from published data."))

        suma <- summary(gene%in%geneAnnot[,col])
        if(length(suma)==3) {  names(suma) <- c("Mode", "Gene not found", "Found")  }
        if(length(suma)==2 & names(suma)[2]=="TRUE") {  names(suma) <- c("Mode", "Found") }
        if(length(suma)==2 & names(suma)[2]=="FALSE") {  names(suma) <- c("Mode", "Gene not found") }
        print(suma)
      }

    }}}}}}}




  message("-------------------------")
if(length(c>0)){
  message(paste("Checking colnames of", names(Metadata)[m[1]] ,"in Metadata in clinical data."))
  message("-------------------------")

  for (i in c){

    if(attributes(Metadata)$Omics.type!="Single.Cell"){

    if(all(ColN%in%as.matrix(Metadata[[i]]))) {   message(paste(MetaDataN[i]), " : PASS") } else {  message(paste(MetaDataN[i]), " : FAIL") }
    }


    if(attributes(Metadata)$Omics.type=="Single.Cell"){

      message(paste("Patients from Single.Cell data:", names(Metadata)[i]))
      tot=0
      for (z in rownames(Metadata[[i]])) {
        t = summary(str_detect(pattern = z, ColN))["TRUE"][1]

        if(is.na(as.numeric(t))){ t = 0}

        tot=tot+as.numeric(t)


      message(c(z," N= ",as.numeric(t)))

      }
      message("Total = " , tot, "\nAre all Patients found in Expression matrix ? ", tot/length(ColN)==1)
      message("-------------------------")



     p =  which(attributes(Metadata)$Data.Type=="SamplesAnnot" & attributes(Metadata)$Export=="No")

     if(length(p)>0){

      message(paste("Samples from Single.Cell data:", names(Metadata)[p]))

       tot=0
       for (z in rownames(Metadata[[p]])) {
         t = summary(str_detect(pattern = z, ColN))["TRUE"][1]

         if(is.na(as.numeric(t))){ t = 0}

         tot=tot+as.numeric(t)


         message(c(z," N= ",as.numeric(t)))

       }
       message("Total = " , tot, "\nAre all Patients found in Expression matrix ? ", tot/length(ColN)==1)
       message("-------------------------")
   }

    }
  }}









}
