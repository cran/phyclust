.First.lib <- function(lib, pkg)
{
  library(ape)
  library.dynam("phyclust", pkg, lib)
} # End of .First.lib()

.Last.lib <- function(libpath)
{
#  ret <- .dynLibs()
#  for(dynlib in ret){
#    if(dynlib[1] == "phyclust"){
#      phyclust.libpath <- gsub("/libs/phyclust\\..*", "", dynlib[2])
#      library.dynam.unload("phyclust", phyclust.libpath)
#    }
#  }

  library.dynam.unload("phyclust", libpath)
} # End of .Last.lib()
