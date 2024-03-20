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
  if(is.null(attributes(Metadata)$Omics.type)){attributes(Metadata)$Omics.type="NotDefine"}
  if("CellsAnnot"%in%attributes(Metadata)$Data.Type){attributes(Metadata)$Omics.type="Single.Cell"}

  l <-length(names(Metadata))
  MetaDataN <- names(Metadata)


  m <- which(attributes(Metadata)$Data.Type=="Count")

  g =  which(attributes(Metadata)$Data.Type=="geneAnnot")[1]
  if(!is.na(g)){geneAnnot = as.matrix(Metadata[[g]])}


  if(attributes(Metadata)$Omics.type!="Single.Cell"){
  c <- which(c(attributes(Metadata)$Data.Type=="Clinic" | attributes(Metadata)$Data.Type=="SamplesAnnot")  & attributes(Metadata)$Cleaned=="No")

  c2 <- which(c(attributes(Metadata)$Data.Type=="Clinic" |attributes(Metadata)$Data.Type=="SamplesAnnot")  & attributes(Metadata)$Cleaned=="Yes")

  if(length(c2)>0){c=c2}

  }

  if(attributes(Metadata)$Omics.type=="Single.Cell"){ c <- which(attributes(Metadata)$Data.Type=="Clinic")}



  if(attributes(Metadata)$Data.Type[c[1]]=="SamplesAnnot"){sID <- Metadata[[c[1]]][,"SamplesID"] }
  if(attributes(Metadata)$Data.Type[c[1]]=="Clinic"){sID <- Metadata[[c[1]]][,"PatientsID"] }

  if(attributes(Metadata)$Data.Type[c[1]]=="Clinic" & attributes(Metadata)$Omics.type=="Single.Cell"){
    cellannot = which(attributes(Metadata)$Data.Type=="CellsAnnot"& attributes(Metadata)$Cleaned=="No")
    cellannot2 = which(attributes(Metadata)$Data.Type=="CellsAnnot"& attributes(Metadata)$Cleaned=="Yes")
    c <- which(attributes(Metadata)$Data.Type=="Clinic"  & attributes(Metadata)$Cleaned=="No")
    c2 <- which(attributes(Metadata)$Data.Type=="Clinic"  & attributes(Metadata)$Cleaned=="Yes")

    cellID <- Metadata[[cellannot[1]]][,"CellsBarcode"]
    if(length(c2>0)){     sID <- Metadata[[c2[1]]][,"PatientsID"];c=c2} else {     sID <- Metadata[[c[1]]][,"PatientsID"]}
    if(length(cellannot2>0)){ cellID <- Metadata[[cellannot2[1]]][,"CellsBarcode"];cellannot=cellannot2} else {   cellID <- Metadata[[cellannot[1]]][,"CellsBarcode"]}

    }



  if(attributes(Metadata)$Data.Type[c[1]]=="SamplesAnnot" & attributes(Metadata)$Omics.type=="Single.Cell" & !"Clinic" %in% attributes(Metadata)$Data.Type){
    cellannot = which(attributes(Metadata)$Data.Type=="CellsAnnot"& attributes(Metadata)$Cleaned=="No")
    cellannot2 = which(attributes(Metadata)$Data.Type=="CellsAnnot"& attributes(Metadata)$Cleaned=="Yes")
    c <- which(attributes(Metadata)$Data.Type=="SamplesAnnot"  & attributes(Metadata)$Cleaned=="No")
    c2 <- which(attributes(Metadata)$Data.Type=="SamplesAnnot"  & attributes(Metadata)$Cleaned=="Yes")

    cellID <- Metadata[[cellannot[1]]][,"CellsBarcode"]
    if(length(c2>0)){     sID <- Metadata[[c2[1]]][,"SamplesID"]} else {     sID <- Metadata[[c[1]]][,"SamplesID"]}
    if(length(cellannot2>0)){ cellID <- Metadata[[cellannot2[1]]][,"CellsBarcode"]} else {   cellID <- Metadata[[cellannot[1]]][,"CellsBarcode"]}




    }




  for (i in m){


if(attributes(Metadata)$Omics.type!="Single.Cell"){
  message("-------------------------")

  ccc = which(attributes(Metadata)$Cleaned=="Yes"& attributes(Metadata)$Data.Type!="geneAnnot")

  if(length(ccc)>0){  message(paste("Checking SamplesID in Cleaned Metadata sub-objects from", names(Metadata)[c2[1]]))}else {
    message(paste("Checking SamplesID in Metadata sub-objects from", names(Metadata)[c[1]])) }

  message("-------------------------")




    if(all(sID %in% colnames(Metadata[[i]]))==T) {   message(paste(MetaDataN[i]), " colnames : PASS") } else {
      message(paste(MetaDataN[i]), " colnames : FAIL")
      if(summary(sID %in% colnames(Metadata[[i]]))["TRUE"]==ncol(Metadata[[i]]) ){message(paste("All samples from", MetaDataN[i],"are found in Samples or clinical annotation file."))}
      message(paste("Samples not found in ", MetaDataN[i]," : "), paste0(na.omit(sID[!sID%in%colnames(Metadata[[i]])]),collapse = "; "))}
}


    if(attributes(Metadata)$Omics.type=="Single.Cell"){
      message("-------------------------")

      ccc = which(attributes(Metadata)$Cleaned=="Yes"& attributes(Metadata)$Data.Type!="geneAnnot")

      if(length(ccc)>0 ){  message(paste("Checking  Cells barcodes and Samples/Patients correspondances in cleaned Metadata sub-objects from", names(Metadata)[c2[1]]))}else {
      message(paste("Checking Cells barcodes and Samples/Patients correspondances in Metadata sub-objects from", names(Metadata)[c[1]]))}
      message("-------------------------")


      if(all(cellID %in% colnames(Metadata[[i]]))==T) {   message(paste("All Cells barcodes in",MetaDataN[i], "colnames : PASS")) } else {

          message(paste(MetaDataN[i]), " colnames : FAIL")

          if(length(!cellID%in%colnames(Metadata[[i]]))>10){
            mismatch = na.omit(cellID[!cellID%in%colnames(Metadata[[i]])])[1:10]
          ext = paste("\n showing ten of",length(cellID[!cellID%in%colnames(Metadata[[i]])]), "not found") }else {
            mismatch=na.omit(cellID[!cellID%in%colnames(Metadata[[i]])])
            ext = ""}
          message(paste("Cells barcodes not found in ", MetaDataN[i]," : "), paste0(mismatch,collapse = "; "),ext)}


    }

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

  ccc = which(attributes(Metadata)$Cleaned=="Yes"& attributes(Metadata)$Data.Type!="geneAnnot")

  if(length(ccc)>0){  message(paste("Checking Common Samples from", names(Metadata)[c2[1]] ,"in other Cleaned Samples or Patients annotations data."))}else {
  message(paste("Checking Common Samples from", names(Metadata)[c[1]] ,"in other Samples or Patients annotations data."))}
  message("-------------------------")

  for (i in c){

    if(attributes(Metadata)$Omics.type!="Single.Cell"){

    if(all(sID%in%as.matrix(Metadata[[i]]))) {   message(paste(MetaDataN[i]), " : PASS") } else {

      message(paste(MetaDataN[i]), " : FAIL")
      message(paste("Samples not found in ", MetaDataN[i]," : "), paste0(na.omit(sID[!sID%in%as.matrix(Metadata[[i]])]),collapse = "; "))


      }
    }


    if(attributes(Metadata)$Omics.type=="Single.Cell"){

      message(paste("Patients from Single.Cell data:", names(Metadata)[i]))
      tot=0
      for (z in rownames(Metadata[[i]])) {
        t = summary(str_detect(pattern = paste0(z,"_"), cellID))["TRUE"][1]

        if(is.na(as.numeric(t))){ t = 0}

        tot=tot+as.numeric(t)


      message(c(z," N= ",as.numeric(t)))

      }
      message("Total = " , tot, "\n Are all cells barcode associated to Patients found in clinical data ? ", tot/length(cellID)==1)
      message("-------------------------")



     p =  which(attributes(Metadata)$Data.Type=="SamplesAnnot" & attributes(Metadata)$Export=="No")

     if(length(p)>0){

      message(paste("Samples from Single.Cell data:", names(Metadata)[p]))

       tot=0
       for (z in rownames(Metadata[[p]])) {
         t = summary(str_detect(pattern = paste0(z,"_"), cellID))["TRUE"][1]

         if(is.na(as.numeric(t))){ t = 0}

         tot=tot+as.numeric(t)


         message(c(z," N= ",as.numeric(t)))

       }
       message("Total = " , tot, "\nAre all Cells barcode associated to Samples in SamplesAnnot ? ", tot/length(cellID)==1)
       message("-------------------------")
     }



     mm =  which(attributes(Metadata)$Data.Type=="Count")

     if(length(mm)>0){

       message(paste("Samples from Single.Cell data:", names(Metadata)[mm]))

       tot=0
       for (z in sID) {
         t = summary(str_detect(pattern = paste0(z,"_"), colnames(Metadata[[mm]])))["TRUE"][1]

         if(is.na(as.numeric(t))){ t = 0}

         tot=tot+as.numeric(t)


         message(c(z," N= ",as.numeric(t)))

       }
       message("Total = " , tot, "\nAre all cells barcode found in Matrix count ? ", tot/length(cellID)==1)
       message("-------------------------")
     }

    }
  }}









}
