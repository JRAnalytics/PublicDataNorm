
#' CheckCleaning function for checking numerical values in Clean clinical data, as all numerical.
#'
#' @param Metadata
#'.
#' @export
#'
#' @examples "non"
CheckCleaning <- function(Metadata) {
for (j in which(c(attributes(Metadata)$Data.Type =="SamplesAnnot" | attributes(Metadata)$Data.Type =="Clinic") & attributes(Metadata)$Export =="Yes" )) {

  {res <- list()

  for (i in colnames(Metadata[[j]])) {

    s <- suppressWarnings(list(summary(as.factor(!is.na(as.numeric(as.character(Metadata[[j]][,i])))))))
    names(s) <- i
    res <- c(res, s)

  }


  missValue <- data.frame("Colnames" = as.character(),
                          "Row" = as.numeric(),
                          "Values" = as.character())
  for (z in names(res)){


    if(length(res[[z]][names(res[[z]])=="TRUE"])!=0){

      if(res[[z]][names(res[[z]])=="TRUE"]!=sum(res[[z]])){


        warning(paste("In numerical data:",z,  ", a part is not numeric."))



        row <- suppressWarnings(which(as.factor(!is.na(as.numeric(as.character(Metadata[[j]][,z]))))==F))


        res2 <- data.frame("Colnames" = z,
                           "Row" = row,
                           "Values" = Metadata[[j]][row,z])

        missValue <-  rbind(missValue,res2)


      }
    }}


  missValue[missValue=="NA"] <- NA

  if(nrow(missValue)>0){
    message(paste(names(Metadata[j])), " MissCalled values")

    return(missValue)
  }
  }
}
}
