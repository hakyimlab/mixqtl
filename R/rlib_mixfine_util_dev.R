#' @importFrom emma emma.MLE.noX
#' @importFrom emma emma.MLE
calc_var = function(y, x, fixed_effect = NULL) {
  kinship = x %*% t(x)
  if(is.null(fixed_effect)) {
    res = emma:::emma.MLE.noX(y, kinship)
    return(res$ve)
  } else {
    res = emma::emma.MLE(y, fixed_effect, kinship)
    return(res$ve)
  }
}
