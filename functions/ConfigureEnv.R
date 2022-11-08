ConfigureEnv <- function() {
  ## Load file configuration, install and load packages needed to run the application
  
  start <- Sys.time()
  
  message('Configuring and checking environment...')
  
  libs <- c("ncdf4", "ggplot2", "gridExtra", "mapproj", "reliabilitydiag")
  
  # Installation of all packages specified in the config file ('config.yml')
  sapply(libs, function(x) 
    if (!(x %in% rownames(installed.packages()))) { 
      print(paste('Installing ', x, ' package...'))
      install.packages(x, clean=TRUE, verbose=FALSE, quiet=TRUE) 
    })
  
  sapply(libs, function(x) library(x, character.only=TRUE))
  
  # Source all functions used in the main program
  #sapply(list.files(pattern="[.]R$", path="functions/", full.names=TRUE), source);
  message('Environment configured sucessfully!')
  
  message("Execution time \"ConfigureEnv.R\": ", round(Sys.time() - start, 4), " secs")
  
}