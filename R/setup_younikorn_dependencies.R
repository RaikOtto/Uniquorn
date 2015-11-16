#' Install packages required by Younikorn
#' @export
prepare_younikorn_dependencies = function(){
  
  required_packages = c( "dplyr","readxl")
  source("https://bioconductor.org/biocLite.R")
  
  pkgTest <- function( x ){
    
    if ( ! require( x, character.only = TRUE ) ){
      
      install.packages( x, dep = T )
      biocLite( x )
      
      if( !require( x, character.only = T ) ) stop( 
        paste0(
          paste0(
            "Package ", 
            x
          ),
          " not found"
        )
      )
    }
  }
  
  lapply( required_packages, FUN = pkgTest )
}