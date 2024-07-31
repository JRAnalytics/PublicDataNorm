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
  c <- which(attributes(Metadata)$Data.Type=="Clinic" & attributes(Metadata)$Cleaned=="No")
  c2 <- which(attributes(Metadata)$Data.Type=="Clinic" & attributes(Metadata)$Cleaned=="Yes")
  if(length(c2)>0){c=c2}
  if(length(c)==0){stop("A Patients' Clinical data must be loaded.")}

  PpID <- unique(Metadata[[c[1]]][,"PatientsID"])
  PsID = Metadata[[c[1]]][,"SamplesID"]
  PsID = unique(unlist(strsplit(PsID, ";")))


  s <- which(attributes(Metadata)$Data.Type=="SamplesAnnot" & attributes(Metadata)$Cleaned=="No")
  s2 <- which(attributes(Metadata)$Data.Type=="SamplesAnnot" & attributes(Metadata)$Cleaned=="Yes")


  if(length(s2)>0){s=s2}

  SsID <- unique(Metadata[[s[1]]][,"SamplesID"])
  SpID <- unique(Metadata[[s[1]]][,"PatientsID"])

    if(length(s)!=0){
    if(attributes(Metadata)$Data.Type[s[1]]=="SamplesAnnot"){sID <- Metadata[[s[1]]][,"SamplesID"] }}
  if(length(c)!=0){
    if(attributes(Metadata)$Data.Type[c[1]]=="Clinic"){pID <- Metadata[[c[1]]][,"PatientsID"] }}

  }





  if(attributes(Metadata)$Omics.type=="Single.Cell"){

    c= NULL
    if("Clinic" %in%attributes(Metadata)$Data.Type){

      c <- which(attributes(Metadata)$Data.Type=="Clinic" & attributes(Metadata)$Cleaned=="No")
      c2 <- which(attributes(Metadata)$Data.Type=="Clinic" & attributes(Metadata)$Cleaned=="Yes")
      if(length(c2)>0){c=c2}
      if(length(c)==0){stop("A Patients' Clinical data must be loaded")}


      PpID <- unique(Metadata[[c[1]]][,"PatientsID"])
      PsID = Metadata[[c[1]]][,"SamplesID"]
      PsID = unique(unlist(strsplit(PsID, ";")))
    }




    if("SamplesAnnot" %in%attributes(Metadata)$Data.Type ){
      s <- which(attributes(Metadata)$Data.Type=="SamplesAnnot" & attributes(Metadata)$Cleaned=="No")
      s2 <- which(attributes(Metadata)$Data.Type=="SamplesAnnot" & attributes(Metadata)$Cleaned=="Yes")
      if(length(s2)>0){s=s2}

      SsID <- unique(Metadata[[s[1]]][,"SamplesID"])
      SpID <- unique(Metadata[[s[1]]][,"PatientsID"])
      }


    if("CellsAnnot" %in%attributes(Metadata)$Data.Type){
      cellannot <- which(attributes(Metadata)$Data.Type=="CellsAnnot" & attributes(Metadata)$Cleaned=="No")
      cellannot2 <- which(attributes(Metadata)$Data.Type=="CellsAnnot" & attributes(Metadata)$Cleaned=="Yes")

      if(length(cellannot2)>0){cellannot=cellannot2}

      cellID <- Metadata[[cellannot[1]]][,"CellsBarcode"]

      }






  }







if(attributes(Metadata)$Omics.type!="Single.Cell"){
  message("-------------------------")

  ccc = which(attributes(Metadata)$Cleaned=="Yes"& attributes(Metadata)$Data.Type!="geneAnnot")

  if(length(ccc)>0){  message(paste("Checking SamplesID in Cleaned Metadata sub-objects from", names(Metadata)[s2[1]]))}else {
    message(paste("Checking SamplesID in Metadata sub-objects from", names(Metadata)[s[1]])) }

  message("-------------------------")


  for (i in m){

    if(length(s)>0 & length(c)<1 | length(C)>0 & length(s)>0 ){
    if(all(SsID %in% colnames(Metadata[[i]]))==T) {   message(paste(MetaDataN[i]), " colnames : PASS") } else {
      message(paste(MetaDataN[i]), " colnames : FAIL")
      if(summary(PsID %in% colnames(Metadata[[i]]))["TRUE"]==ncol(Metadata[[i]]) ){message(paste("All samples from", MetaDataN[i],"are found in Samples or clinical annotation file."))}
      message(paste("Samples not found in ", MetaDataN[i]," : "), paste0(na.omit(SsID[!SsID%in%colnames(Metadata[[i]])]),collapse = "; "))}
    }
    if(length(C)>0 & length(s)<1  ){
      if(all(PsID %in% colnames(Metadata[[i]]))==T) {   message(paste(MetaDataN[i]), " colnames : PASS") } else {
        message(paste(MetaDataN[i]), " colnames : FAIL")
        if(summary(PsID %in% colnames(Metadata[[i]]))["TRUE"]==ncol(Metadata[[i]]) ){message(paste("All samples from", MetaDataN[i],"are found in Samples or clinical annotation file."))}
        message(paste("Samples not found in ", MetaDataN[i]," : "), paste0(na.omit(PsID[!PsID%in%colnames(Metadata[[i]])]),collapse = "; "))}
    }

    }
  }


    if(attributes(Metadata)$Omics.type=="Single.Cell"){
      message("-------------------------")

      ccc = which(attributes(Metadata)$Cleaned=="Yes"& attributes(Metadata)$Data.Type!="geneAnnot")

      if(length(ccc)>0 ){  message(paste("Checking  Cells barcodes and Samples/Patients correspondances in cleaned Metadata sub-objects from", names(Metadata)[c2[1]]))}else {
      message(paste("Checking Cells barcodes and Samples/Patients correspondances in Metadata sub-objects from", names(Metadata)[c[1]]))}
      message("-------------------------")


      for (i in m){
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


if(length(c)>0){ccc = which(attributes(Metadata)$Cleaned=="Yes"& attributes(Metadata)$Data.Type!="geneAnnot")}
mm =  which(attributes(Metadata)$Data.Type=="Count")


  if(length(ccc)>0){  message(paste("Checking Common Patients from", names(Metadata)[c[1]] ,"in other Cleaned Samples or Patients annotations data."))
    }else {
  message(paste("Checking Common Patients from", names(Metadata)[c[1]] ,"in other Samples or Patients annotations data."))}
  message("-------------------------")






    if(attributes(Metadata)$Omics.type!="Single.Cell"){
      for (i in c(c[-1],s)){
      target = unique(Metadata[[i[1]]][,"PatientsID"])

      if(length(which(PpID %in% as.matrix(Metadata[[i[1]]])))==length(PpID)){message(paste(MetaDataN[i]), " : PASS") }
      if(length(which(PpID %in% as.matrix(Metadata[[i[1]]])))<length(PpID)){
        message(paste(MetaDataN[i]), " : FAIL")
        message(paste("PatientsID not found in ", MetaDataN[i]," : "), paste0(na.omit(PpID[!PpID%in%target]),collapse = "; "))
        }


      }}



    if(attributes(Metadata)$Omics.type=="Single.Cell"){
      for (i in c(c,s)){
      message(paste0("PatientsID from '", names(Metadata)[i],"', in CellsAnnotation object"))
      tot=0
      for (z in unique(Metadata[[i]][,"PatientsID"])) {
        if(str_detect(pattern = paste0('[a-zA-Z]'),z)){ pattern = paste0(z,"-")}else {pattern=  paste0('[a-zA-Z]',z,"-")}
        t = summary(str_detect(pattern = pattern, cellID))["TRUE"][1]

        if(is.na(as.numeric(t))){ t = 0}

        tot=tot+as.numeric(t)



      }
      message("Total = " , tot,"/",length(cellID), "\n Passed Checkpoint? ", tot/length(cellID)==1)
      message("-------------------------")}


      if(length(ccc)>0){ p =  which(attributes(Metadata)$Data.Type=="SamplesAnnot" & attributes(Metadata)$Cleaned=="Yes")}else{
        p =  which(attributes(Metadata)$Data.Type=="SamplesAnnot" & attributes(Metadata)$Cleaned=="No")
      }

      if(length(p)>0){

        for (i in p){
          message(paste0("SamplesID from '", names(Metadata)[i],"', in CellsAnnotation object :"))

          tot=0
          for (z in Metadata[[i]][,"SamplesID"]) {
            if(str_detect(pattern = paste0('[a-zA-Z]'),z)){ pattern = paste0(z,"-")}else {  paste0('[a-zA-Z]',z,"-")}
            t = summary(str_detect(pattern = pattern, cellID))["TRUE"][1]

            if(is.na(as.numeric(t))){ t = 0}

            tot=tot+as.numeric(t)


          }
          message("Total = " , tot, "/",length(cellID), "\n Passed Checkpoint? ", tot/length(cellID)==1)
          message("-------------------------")
        }
      }



      }













}











