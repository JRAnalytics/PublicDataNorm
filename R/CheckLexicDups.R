#' CheckLexicDups
#'
#' @param Lexic Lexic to check if duplicated values (Targeted column) are in the lexic items.
#'
#' @return nothing
#' @export
#'
#' @examples "non"
CheckLexicDups = function(Lexic) {



DupItems = duplicated(toupper(names(Lexic)))
  if(!all(DupItems==F)){

    DupItems = names(Lexic)[which(DupItems)]

    message(paste("Duplicated Items '",DupItems,"' in lexic"))}

dupInItems = unlist(lapply(Lexic, function(x){x = x %>% duplicated() %>% which()}))

if(length(dupInItems)>1){

  dupInItemsTargets = names(Lexic)[unlist(lapply(names(dupInItems), function(x) which(str_detect(x, names(Lexic)))))]
  dt = data.frame("Position" =dupInItems ,"LexicItem" =dupInItemsTargets, "Dup" = NA)

  for (i in rownames(dt)){

    dt[i,]$Dup = Lexic[[ dt[i,]$LexicItem ]][ dt[i,]$Position]
  }
  dt$DupDup = paste0(dt$LexicItem,dt$Dup)
  dt = dt[-which(duplicated(dt$DupDup)),]

  message(c(paste0("Duplicated column target(s) '",dt$Dup,"' in '", dt$LexicItem,"'.\n")))}

if(length(dupInItems)==1){
    dupInItemsTargets = Lexic[[names(dupInItems)]][dupInItems]
    names(dupInItemsTargets) = names(dupInItems)
    message(c(paste0("Duplicated column target(s) '",dupInItemsTargets,"' in '", names(dupInItems),"'.\n")))

    }





lexiceclate = unlist(lapply(Lexic, function(x){x = x}))

dupInLexic= lexiceclate[which(duplicated(lexiceclate))]

if(length(dupInItems)==1){
dupInLexic = dupInLexic[dupInLexic!=dupInItemsTargets]}

if(length(dupInItems)>1){dupInLexic = dupInLexic[dupInLexic!=dt$Dup]}


NamesdupInLexic = names(lexiceclate)[lexiceclate %in% dupInLexic]

OriItem = NamesdupInLexic[ !NamesdupInLexic %in%  names(dupInLexic)]
OriItem = names(Lexic)[unlist(lapply(OriItem, function(x) which(str_detect(x, names(Lexic)))))]
Dupinitem = names(Lexic)[unlist(lapply(names(dupInLexic), function(x) which(str_detect(x, names(Lexic)))))]



if(length(dupInLexic)>0){
  message(paste0("column target '",dupInLexic,"' from Lexic's '", OriItem, "' is duplicated in '", Dupinitem,"'\n" ))
  }


}



