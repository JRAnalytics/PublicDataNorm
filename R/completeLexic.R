
#' completeLexic
#'
#' @param liste_de_listes Lexic
#' @param noms_listes colnames(clinic)
#'
#' @return un lexic
#' @export
#'
#' @examples none
#'
completeLexic=function(liste_de_listes, noms_listes){
  for (nom_liste in noms_listes) {
    add = FALSE
    for (liste in liste_de_listes){

      key = liste[1]

      if (toupper(nom_liste) %in% toupper(liste)) {
        # Ajouter l'élément à la liste existante
        liste_de_listes[[key]] <- unique(c(liste_de_listes[[key]], nom_liste))
        add = TRUE
      }
    }

    if (!add) {
      # Créer une nouvelle liste avec l'élément
      liste_de_listes[[nom_liste]] <- unique(nom_liste)
    }

  }

  return(liste_de_listes)
}
