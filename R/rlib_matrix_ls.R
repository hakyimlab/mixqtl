#' trcQTL for nominal pass in matrix form
#'
#' Perform trcQTL in matrix form which can be used for obtaining nominal p-value
#'
#' @param trc total read count (dimension = N x 1)
#' @param lib_size library size (dimension = N x 1)
#' @param x 'genotype' used in trcQTL (which is defined as (x1 + x2) / 2) (dimension = N x P)
#' @param cov effect of covariates on log(trc/2) (dimension = N x 1)
#' @param trc_cutoff total read count cutoff to exclude observations with trc lower than the cutoff
#'
#' @return a list of summary statistics
#'         beta_hat: estimated log aFC (dimension = 1 x P)
#'         beta_se: standard deviation of log aFC (dimension = 1 x P)
#'         sample_size: number of observations used in regression
#'
#' @examples
#' matrix_ls_trc(
#'   trc = rpois(100, 100),
#'   lib_size = rpois(100, 10000),
#'   x = matrix(sample(c(0, 0.5, 1), 200, replace = TRUE), ncol = 2),
#'   cov = rnorm(100),
#'   trc_cutoff = 20
#' )
#'
#' @export
#' @importFrom stats var
matrix_ls_trc = function(trc, lib_size, x, cov, trc_cutoff = 20) {
  trc_in = trc
  trc = log(trc / 2 / lib_size) - cov
  snp_id = 1 : ncol(x)
  output = list(beta_hat = matrix(NA, ncol = length(snp_id), nrow = 1), beta_se = matrix(NA, ncol = length(snp_id), nrow = 1))
  na_ind = is.na(trc) | trc_in < trc_cutoff
  x = x[!na_ind, , drop = F]
  trc = trc[!na_ind, drop = F]
  mono_ind = apply(x, 2, var) == 0
  x_filtered = x[, !mono_ind, drop = F]
  sample_size = nrow(x_filtered)
  intercept_mat = matrix(1, nrow = nrow(x_filtered), ncol = ncol(x_filtered))
  if(sample_size > 2 & ncol(x_filtered) > 0) {
    out = matrixQTL_two_dim(matrix(trc, ncol = 1), x_filtered, intercept_mat, matrix(sample_size, ncol = ncol(x_filtered), nrow = 1))
    output$beta_hat[, !mono_ind] = t(out$beta1_hat)
    output$beta_se[, !mono_ind] = t(out$beta1_se)
  }
  output$sample_size = sample_size
  return(output)
}

#' ascQTL for nominal pass in matrix form
#'
#' Perform ascQTL in matrix form which can be used for obtaining nominal p-value
#'
#' @param asc1 allele-specific read count for haplotype 1 (dimension = N x 1)
#' @param asc2 allele-specific read count for haplotype 2 (dimension = N x 1)
#' @param x 'genotype' used in ascQTL (which is defined as x1 - x2) (dimension = N x P)
#' @param asc_cutoff allele-specific read count cutoff to exclude observations with asc1 or asc2 lower than asc_cutoff
#' @param weight_cap the maximum weight difference (in fold) is min(weight_cap, floor(sample_size / 10)). The ones exceeding the cutoff is capped.
#' @param asc_cap exclude observations with asc1 or asc2 higher than asc_cap
#'
#' @return a list of summary statistics
#'         beta_hat: estimated log aFC (dimension = 1 x P)
#'         beta_se: standard deviation of log aFC (dimension = 1 x P)
#'         sample_size: number of observations used in regression
#'
#' @examples
#' matrix_ls_asc(
#'   asc1 = rpois(100, 100),
#'   asc2 = rpois(100, 100),
#'   x = matrix(sample(c(0, 0.5, 1), 200, replace = TRUE), ncol = 2),
#'   asc_cutoff = 20,
#'   weight_cap = 10,
#'   asc_cap = 1000
#' )
#'
#' @export
#' @importFrom stats var
matrix_ls_asc = function(asc1, asc2, x, asc_cutoff = 5, weight_cap = 100, asc_cap = 5000) {
  snp_id = 1 : ncol(x)
  output = list(beta_hat = matrix(NA, ncol = length(snp_id), nrow = 1), beta_se = matrix(NA, ncol = length(snp_id), nrow = 1))
  passed_ind = asc1 >= asc_cutoff & asc2 >= asc_cutoff & asc1 <= asc_cap & asc2 <= asc_cap
  x = x[passed_ind, , drop = F]
  asc1 = asc1[passed_ind, drop = F]
  asc2 = asc2[passed_ind, drop = F]
  asc = log(asc1 / asc2)
  mono_ind = apply(x, 2, var) == 0
  x_filtered = x[, !mono_ind, drop = F]
  sample_size = nrow(x_filtered)
  weights = harmonic_sum_(asc1, asc2)
  weight_cap = min(weight_cap, floor(sample_size / 10))
  weight_cutoff = min(weights) * weight_cap
  weights[weights > weight_cutoff] = weight_cutoff

  if(sample_size > 2 & ncol(x_filtered) > 0) {
    asc = diag(sqrt(weights)) %*% asc
    x_filtered = diag(sqrt(weights)) %*% x_filtered
    out = matrixQTL_one_dim(matrix(asc, ncol = 1), x_filtered, matrix(sample_size, ncol = ncol(x_filtered), nrow = 1))
    output$beta_hat[, !mono_ind] = t(out$beta_hat)
    output$beta_se[, !mono_ind] = t(out$beta_se)
  }
  output$sample_size = sample_size
  return(output)
}

#' Solve Y = Xb + e in matrix form
#'
#' For each column i in X, solve Y_k = X_i b_i + e as least squares problem and output estimated effect size and standard deviation
#'
#' @param Y response to regress against (dimension = N x K)
#' @param X P predictors to perform regression separately (dimension = N x P)
#' @param n sample size (dimension = P x K)
#'
#' @return a list of summary statistics
#'         beta_hat: estimated b, b_hat (dimension = K x P)
#'         beta_se: standard deviation of b_hat (dimension = K x P)
#'
#' @examples
#' matrixQTL_one_dim(
#'   Y = matrix(rnorm(300), ncol = 3),
#'   X = matrix(sample(c(0, 0.5, 1), 200, replace = TRUE), ncol = 2),
#'   n = matrix(100, ncol = 2, nrow = 3)
#' )
#'
#' @export
#' @import tensorA
matrixQTL_one_dim = function(Y, X, n) {
  # input: n -> i, p -> j, k -> k
  # Y: n by k
  # X: n by p
  # n: n, i.e. sample size
  # output:
  # for each y_k x_p pair run y ~ -1 + x
  # betahat: p by k
  # betase: p by k
  if(dim(Y)[1] != dim(X)[1] | dim(n)[2] != dim(X)[2]) {
    message('Wrong dimension')
    return(NULL)
  }
  X_tensor <- to.tensor(as.vector(X), c(i = dim(X)[1], j = dim(X)[2]))
  # Y_tensor <- to.tensor(as.vector(Y), c(i = dim(Y)[1], j = dim(Y)[2]))
  Y_tensor <- to.tensor(as.vector(Y), c(i = dim(Y)[1], k = dim(Y)[2]))
  n_tensor <- to.tensor(as.vector(n), c(j = dim(n)[2], k = dim(n)[1]))
  YtX_tensor = einstein.tensor(X_tensor, by = 'j', Y_tensor)
  # YtX_tensor = mul.tensor(X_tensor, i = 'i', Y = Y_tensor)
  XtX_tensor = einstein.tensor(X_tensor, by = 'j', X_tensor)
  hatbeta_tensor = YtX_tensor / XtX_tensor

  # R_tensor = Y_tensor - diagmul.tensor(X_tensor, D = hatbeta_tensor, j = 'j', i = 'j')
  # Rsq_tensor = einstein.tensor(R_tensor, by = c('j', 'k'), R_tensor)
  # new way to calculate Rsq_tensor
  YtY_tensor = einstein.tensor(Y_tensor, by = 'k', Y_tensor)
  Rsq_tensor = YtY_tensor - 2 * hatbeta_tensor * YtX_tensor + hatbeta_tensor^2 * XtX_tensor

  hatsigma_tensor = sqrt(Rsq_tensor / (n_tensor - 1))
  sebeta_tensor = hatsigma_tensor / sqrt(XtX_tensor)
  # if(dim(Y)[2] == 1) {
  #   return(data.frame(beta_hat = as.vector(hatbeta_tensor), beta_se = as.vector(sebeta_tensor)))
  # } else {
  return(list(beta_hat = to.matrix.tensor(hatbeta_tensor, 'j', 'k'), beta_se = to.matrix.tensor(sebeta_tensor, 'j', 'k')))
  # }
}

#' Solve Y = X1 b1 + X2 b2 + e in matrix form
#'
#' For each column i in X, solve Y_k = X_1i b_1i + X_2i b_2i + e as least squares problem and output estimated effect size and standard deviation
#'
#' @param Y response to regress against (dimension = N x K)
#' @param X1 P predictors (as the first predictor) to perform regression separately (dimension = N x P)
#' @param X2 P predictors (as the second predictor) to perform regression separately (dimension = N x P)
#' @param n sample size (dimension = P x K)
#'
#' @return a list of summary statistics
#'         beta1_hat: estimated b1, b1_hat (dimension = K x P)
#'         beta1_se: standard deviation of b1_hat (dimension = K x P)
#'         beta2_hat: estimated b2, b2_hat (dimension = K x P)
#'         beta2_se: standard deviation of b2_hat (dimension = K x P)
#'
#' @examples
#' matrixQTL_two_dim(
#'   Y = matrix(rnorm(100), ncol = 3),
#'   X1 = matrix(sample(c(0, 0.5, 1), 200, replace = TRUE), ncol = 2),
#'   X2 = matrix(sample(c(0, 0.5, 1), 200, replace = TRUE), ncol = 2),
#'   n = matrix(200, ncol = 2, nrow = 3)
#' )
#'
#' @export
#' @import tensorA
matrixQTL_two_dim = function(Y, X1, X2, n) {
  # input: n -> i, p -> j, k -> k
  # Y: n by k
  # X1: n by p
  # X2: n by p
  # n: n, i.e. sample size
  # output:
  # for each y_k x_p pair run y ~ -1 + x1 + x2
  # betahat1: p by k
  # betase1: p by k
  # betahat2: p by k
  # betase2: p by k
  if(dim(Y)[1] != dim(X1)[1] | dim(n)[2] != dim(X1)[2] | dim(Y)[1] != dim(X2)[1]) {
    message('Wrong dimension')
    return(NULL)
  }
  X1_tensor <- to.tensor(as.vector(X1), c(i = dim(X1)[1], j = dim(X1)[2]))
  X2_tensor <- to.tensor(as.vector(X2), c(i = dim(X2)[1], j = dim(X2)[2]))
  # Y_tensor <- to.tensor(as.vector(Y), c(i = dim(Y)[1], j = dim(Y)[2]))
  Y_tensor <- to.tensor(as.vector(Y), c(i = dim(Y)[1], k = dim(Y)[2]))
  n_tensor <- to.tensor(as.vector(n), c(j = dim(n)[2], k = dim(n)[1]))

  T1_tensor = einstein.tensor(X1_tensor, by = 'j', Y_tensor)
  T2_tensor = einstein.tensor(X2_tensor, by = 'j', Y_tensor)
  # T1_tensor = mul.tensor(X1_tensor, i = 'i', Y = Y_tensor)
  # T2_tensor = mul.tensor(X2_tensor, i = 'i', Y = Y_tensor)

  S11_tensor = einstein.tensor(X1_tensor, by = 'j', X1_tensor)
  S12_tensor = einstein.tensor(X1_tensor, by = 'j', X2_tensor)
  S22_tensor = einstein.tensor(X2_tensor, by = 'j', X2_tensor)

  Delta_tensor = abs(S11_tensor * S22_tensor - S12_tensor * S12_tensor)

  hatbeta1_tensor = (S22_tensor * T1_tensor - S12_tensor * T2_tensor) / Delta_tensor
  hatbeta2_tensor = (S11_tensor * T2_tensor - S12_tensor * T1_tensor) / Delta_tensor

  # R_tensor = Y_tensor - diagmul.tensor(X1_tensor, D = hatbeta1_tensor, j = 'j', i = 'j') - diagmul.tensor(X2_tensor, D = hatbeta2_tensor, j = 'j', i = 'j')
  #
  # Rsq_tensor = einstein.tensor(R_tensor, by = c('j', 'k'), R_tensor)
  # new way to calculate Rsq_tensor
  YtY_tensor = einstein.tensor(Y_tensor, by = 'k', Y_tensor)
  Rsq_tensor = YtY_tensor - 2 * hatbeta1_tensor * T1_tensor - 2 * hatbeta2_tensor * T2_tensor + 2 * hatbeta1_tensor * hatbeta2_tensor * S12_tensor + hatbeta1_tensor^2 * S11_tensor + hatbeta2_tensor^2 * S22_tensor

  hatsigma_tensor = sqrt(Rsq_tensor / (n_tensor - 2))

  sebeta1_tensor = hatsigma_tensor * sqrt(S22_tensor / Delta_tensor)
  sebeta2_tensor = hatsigma_tensor * sqrt(S11_tensor / Delta_tensor)

  # if(dim(Y)[2] == 1) {
  #   return(
  #     data.frame(
  #       beta1_hat = as.vector(hatbeta1_tensor),
  #       beta1_se = as.vector(sebeta1_tensor),
  #       beta2_hat = as.vector(hatbeta2_tensor),
  #       beta2_se = as.vector(sebeta2_tensor)
  #     )
  #   )
  # } else {
  return(
    list(
      beta1_hat = to.matrix.tensor(hatbeta1_tensor, 'j', 'k'),
      beta1_se = to.matrix.tensor(sebeta1_tensor, 'j', 'k'),
      beta2_hat = to.matrix.tensor(hatbeta2_tensor, 'j', 'k'),
      beta2_se = to.matrix.tensor(sebeta2_tensor, 'j', 'k')
    )
  )
  # }
}



