#' Combining total and allele-specific reads for building prediction model
#'
#' Given total and allele-specific read counts along with library size,
#' the two-step inference procedure is implemented to build prediction model.
#' Step 1: estimate variances of total and allele-specific count response respectively
#' along with intercept for total count response
#' Step 2: transform observation according to step 1 and build prediction model
#' using glmnet (Elastic net with alpha = alpha and nfold-fold cross-validation)
#'
#' @param geno1 genotype of haplotype 1 (dimension = N x P)
#' @param geno2 genotype of haplotype 2 (dimension = N x P)
#' @param y1 allele-specific read count for haplotype 1 (dimension = N x 1)
#' @param y2 allele-specific read count for haplotype 2 (dimension = N x 1)
#' @param ytotal total read count (dimension = N x 1)
#' @param lib_size library size (dimension = N x 1)
#' @param cov_offset predicted effect of covariates on total read count response term,
#' log(ytotal / lib_size / 2) (dimension = N x 1)
#' @param trc_cutoff total read count cutoff to exclude observations with ytotal lower than the cutoff
#' @param asc_cutoff allele-specific read count cutoff to exclude observations with y1 or y2 lower than asc_cutoff
#' @param weight_cap the maximum weight difference (in fold) is min(weight_cap, floor(sample_size / 10)). The ones exceeding the cutoff is capped. Set to NULL then no weight_cap applied.
#' @param asc_cap exclude observations with y1 or y2 higher than asc_cap
#' @param alpha alpha parameter in elastic net model of glmnet (lasso: alpha = 1; ridge: alpha = 0). (default = 0.5).
#' @param nfold number of fold for cross-validation to pick lambda parameter in glmnet. (default = 5).
#' @param nobs_asc_cutoff don't consider ASC if number of observations is smaller than nobs_asc_cutoff
#' @param ... Extra args for the last fit_glmnet_with_cv call
#'
#' @return prediction model
#'
#' @examples
#' mixpred(
#'   geno1 = matrix(sample(c(0, 0.5, 1), 200, replace = TRUE), ncol = 2),
#'   geno2 = matrix(sample(c(0, 0.5, 1), 200, replace = TRUE), ncol = 2),
#'   y1 = rpois(100, 100),
#'   y2 = rpois(100, 100),
#'   ytotal = rpois(100, 1000),
#'   lib_size = rpois(100, 10000),
#'   cov_offset = runif(100),
#'   trc_cutoff = 100,
#'   asc_cutoff = 50,
#'   weight_cap = 100,
#'   asc_cap = 1000
#' )
#'
#' @importFrom susieR susie
mixpred_old = function(geno1, geno2, y1, y2, ytotal, lib_size, cov_offset = NULL, trc_cutoff = 20, asc_cutoff = 5, weight_cap = 100, asc_cap = 5000, alpha = 0.5, nfold = 5, nobs_asc_cutoff = 1, ...) {
  # prepare X
  h1 = geno1
  h1[is.na(h1)] = 0.5
  h2 = geno2
  h2[is.na(h2)] = 0.5
  X = rbind(h1 - h2, (h1 + h2) / 2)

  # apply trc and asc cutoffs
  trc_pass_ind = !is.na(ytotal) & (ytotal >= trc_cutoff)
  asc_pass_ind = !is.na(y1) & !is.na(y2) & (y1 >= asc_cutoff) & (y2 >= asc_cutoff) & (y1 <= asc_cap) & (y2 <= asc_cap)

  # prepare log y/lib
  asc = y1 + y2
  trc = ytotal
  if(is.null(cov_offset)) {
    logy = c(
      log(y1 / y2),
      log(trc / 2 / lib_size)
    )
  } else {
    logy = c(
      log(y1 / y2),
      log(trc / 2 / lib_size) - cov_offset
    )
  }
  # intercept
  inpt = c(
    rep(0, length(asc)),
    rep(1, length(trc))
  )
  # prepare weight
  weights = c(
    1 / (1 / y1 + 1 / y2),
    trc  # it is not really used as weight at downstream
  )

  # merge to df
  df = data.frame(y = logy, inpt = inpt, weights = weights)
  pass_qc_ind = c(asc_pass_ind, trc_pass_ind)
  df = df[pass_qc_ind, , drop = F]
  X = X[pass_qc_ind, , drop = F]
  is_na = rowSums(is.na(df) | is.infinite(df$y)) > 0
  df = df[!is_na, , drop = F]
  X = X[!is_na, , drop = F]

  # applying weight cap on ASC rows
  if(!is.null(weight_cap)) {
    weights_asc = df$weights[df$inpt == 0]
    sample_size = sum(df$inpt == 0)
    weight_cap = min(weight_cap, floor(sample_size / 10))
    weight_cutoff = min(weights_asc) * weight_cap
    weights_asc[weights_asc > weight_cutoff] = weight_cutoff
    df$weights[df$inpt == 0] = weights_asc
  }

  # training
  if(sum(df$inpt == 0) >= nobs_asc_cutoff) {
    susie_data1 = approx_susie(X[df$inpt == 0, , drop = F], df$y[df$inpt == 0], w = df$weights[df$inpt == 0], intercept = FALSE)
    # impute effective y and X
    if(is.na(susie_data1$sigma)) {
      X1 = NULL
      y1 = NULL
      susie_data1 = NULL
    } else {
      X1 = diag(sqrt(df$weights[df$inpt == 0])) %*% X[df$inpt == 0, , drop = F] / susie_data1$sigma
      y1 = susie_data1$y
    }
  } else {
    X1 = NULL
    y1 = NULL
  }
  # susie_data1 = approx_susie(X[df$inpt == 0, , drop = F], df$y[df$inpt == 0], w = df$weights[df$inpt == 0], intercept = FALSE)
  susie_data2 = approx_susie(X[df$inpt == 1, , drop = F], df$y[df$inpt == 1], w = NULL, intercept = TRUE)

  if(is.na(susie_data2$sigma)) {
    message('Unexpected failure of fitting sigma for total counts. Too few observations? Quit!')
    quit()
  }

  # impute effective y and X
  # X1 = diag(sqrt(df$weights[df$inpt == 0])) %*% X[df$inpt == 0, , drop = F] / susie_data1$sigma
  X2 = X[df$inpt == 1, , drop = F] / susie_data2$sigma
  X = rbind(X1, X2)
  mod222 = susie(X2, susie_data2$y, standardize = F, intercept = T)
  y2 = susie_data2$y - mod222$intercept
  df = data.frame(y = c(y1, y2))


  # fit elastic net model with imputed y and X
  mod = fit_glmnet_with_cv(X, df$y, nfold = nfold, alpha = alpha, ...)  # intercept = F)
  
  mod
}
