#' Charac.to.Num  transform only numbers column from character to numeric
#'
#' @param clinic clinic data to perform
#' @import readr
#' @return
#' @export
#'
#' @examples
Charac.to.Num=function(clinic) {
  numbers_only <- function(x) { t = gsub("[.]", "", x)
  t = gsub("[-]", "", t )
  !grepl("\\D", t)}

  mixte <- function(x) grepl("\\D", x)
  letters_only <- function(x) !grepl("[^A-Za-z]", x)
  double_only <- function(x) grepl("^[:digit:]+$", x)

  for (i in 1:ncol(clinic)){


    if(all(numbers_only(na.omit(clinic[,i]))==T)& #F
       all(letters_only(na.omit(clinic[,i]))==F)& #F
       all(double_only(na.omit(clinic[,i] ))==F)){

      t = gsub("[.]", "", clinic[,i] )
      t = gsub("[-]", "", t )

      clinic[,i][clinic[,i]=="NA"] = NA

      if(length(na.omit(t))!=0){
        if(all(numbers_only(na.omit(t))== T)){

          clinic[,i][is.na( clinic[,i])] = "NA"
          clinic[,i] = as.numeric(parse_number(clinic[,i], na = "NA" ))
        }}}


    if(all(numbers_only(na.omit(clinic[,i]))==T)& #F
       all(letters_only(na.omit(clinic[,i]))==F)& #F
       all(double_only(na.omit(clinic[,i] ))==F)) {

      clinic[,i] = as.numeric(clinic[,i])
    }
  }

  return(clinic)
}
