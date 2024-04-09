# functions

get.dens.inits <- function(Areas, Counts) {

  Inits <- ((Counts + 0.5) / Areas)
  return(Inits)

}