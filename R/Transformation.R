
#' Transformation of age information
#'
#' @param x Vector of untransformed ages
#' @param method The transformation method to use. Supported values are: "identity"
#'
#' @return Vector with transformed ages
#' @export
#'
transformAge <- function(x,
                         method = "identity"
){
  method <- tolower(method)
  m <- match.arg(method,c("identity"))

  if(m == "identity"){
    return(.identityTransform(x))
  }

  warning(paste0("Method \"",method,"\" is not supported. Returning untransformed data!"))
  return(x)
}

#' Identity transformation
#'
#' @param x Vector of untransformed ages
#'
#' @return Vector with transformed ages
#'
.identityTransform <- function(x
){
  return(x)
}





