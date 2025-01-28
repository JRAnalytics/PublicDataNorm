# .onLoad <- function(libname, pkgname)
# {
#   library.dynam("PublicDataNorm", pkgname, libname)
# }

PDNStartupMessage <- function(){
  # Startup message obtained as
  # > figlet -f slant MCLUST
  msg <- c(paste0(
    "PublicDataNorm version ",
    packageVersion("PublicDataNorm")))

  ver = suppressMessages(check_github("JRAnalytics/PublicDataNorm"))

   if(ver$up_to_date==T){ msg = paste(msg, "\nPackage up to date")}else{

     msg = paste(msg, "\nJRAnalytics/PublicDataNorm repos version", ver$latest_version,
                 "\n please uptdate package using :\n",
                 "devtools::install_github('JRAnalytics/PublicDataNorm', upgrade = 'always')")
   }


  return(msg)
}



.onAttach <- function(lib, pkg){
  # startup message
  msg <- PDNStartupMessage()

  if(!interactive())
    msg[1] <- paste("Package 'PublicDataNorm' version", packageVersion("PublicDataNorm"))
  packageStartupMessage(msg)
  invisible()

}



check_github <- function(pkg) {
  check_github_gitlab(pkg, "github")
}


check_github_gitlab <- function(pkg, repo="github") {
  installed_version <- tryCatch(packageVersion(gsub(".*/", "", pkg)), error=function(e) NA)

  if(repo == "github") {
    url <- paste0("https://raw.githubusercontent.com/", pkg, "/master/DESCRIPTION")
  } else if (repo == "gitlab") {
    url <- paste0("https://gitlab.com/", pkg, "/raw/master/DESCRIPTION")
  } else {
    stop("only work with github and gitlab")
  }

  x <- readLines(url)
  remote_version <- gsub("Version:\\s*", "", x[grep('Version:', x)])

  res <- list(package = pkg,
              installed_version = installed_version,
              latest_version = remote_version,
              up_to_date = NA)

  if (is.na(installed_version)) {
    message(paste("##", pkg, "is not installed..."))
  } else {
    if (remote_version > installed_version) {
      msg <- paste("##", pkg, "is out of date...")
      message(msg)
      res$up_to_date <- FALSE
    } else if (remote_version == installed_version) {
      message("package is up-to-date devel version")
      res$up_to_date <- TRUE
    }
  }

  return(res)
}

