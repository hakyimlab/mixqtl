# approx_susie_no_intercept = function(x, y, w, only_weight = F) {
#   x_susie = mydiag(sqrt(w)) %*% x
#   y_susie = mydiag(sqrt(w)) %*% y
#   if(only_weight == T) {
#     return(list(x = x_susie, y = y_susie))
#   }
#   qr_x = qr(x)
#   x_tilde = qr.Q(qr_x)
#   bhat = solve(x_tilde, y)
#   mod = cv.glmnet(x, y, weights = w, standardize = F, intercept = F)
#   beta = as.vector(mod$glmnet.fit$beta[, which.min(mod$cvm)])
#   sigma = sqrt(mean((y - x %*% beta)^2))
#   list(x = x_susie / sigma, y = y_susie / sigma)
# }

#' @importFrom glmnet cv.glmnet
#' @importFrom stats coef residuals
approx_susie = function(x, y, w = NULL, intercept = T, only_weight = F) {
  if(!is.null(w)) {
    x_susie = mydiag(sqrt(w)) %*% x
    y_susie = mydiag(sqrt(w)) %*% y
    if(intercept == T) {
      intercept_susie = sqrt(w)
      mod = lm(y_susie ~ intercept_susie - 1)
      y_susie = residuals(mod)
    }
    if(only_weight == T) {
      return(list(x = x_susie, y = y_susie))
    }
    sigma = tryCatch(
      {
        mod = glmnet::cv.glmnet(x, y, weights = w, standardize = F, intercept = intercept, nfold = 4)
        # To stablize the estimation of sigma, we calculate MSE using best beta according to CV mean, CV low, CV up and take the max value 
        # beta1 = as.vector(coef(mod$glmnet.fit)[, which.min(mod$cvm)])
        # beta2 = as.vector(coef(mod$glmnet.fit)[, which.min(mod$cvlo)])
        # beta3 = as.vector(coef(mod$glmnet.fit)[, which.min(mod$cvup)])
        beta = as.vector(coef(mod$glmnet.fit)[, which.min(mod$cvm)])
        x_pred = cbind(rep(1, nrow(x)), x)
        sigma = sqrt(sum((y - x_pred %*% beta)^2) / sum(1 / w))
        # sigma1 = sqrt(mean((y - x_pred %*% beta1)^2 * w))
        # sigma2 = sqrt(mean((y - x_pred %*% beta2)^2 * w))
        # sigma3 = sqrt(mean((y - x_pred %*% beta3)^2 * w))
        # sigma = max(c(sigma1, sigma2, sigma3))
         
        # for further debugging
        #
        # mod = glmnet::cv.glmnet(x_susie, y_susie, standardize = F, intercept = intercept)
        # beta = as.vector(coef(mod$glmnet.fit)[, which.min(mod$cvm)])
        # x_pred = cbind(rep(1, nrow(x_susie)), x_susie)
        # sigma = sqrt(mean((y_susie - x_pred %*% beta)^2))
        sigma
      },
      error = function(cond) {
        message('cv.glmnet failed with error:')
        message(cond)
        message(' discard the fit')
        # message('use sigma = sqrt(sum(y_bar^2) / sum(1/w))')
        # mod = lm(y ~ 1, weights = w)
        # y_bar = residuals(mod)
        # sigma = sqrt(sum(y_bar^2) / sum(1 / w))
        sigma = NA
        return(sigma)
      } 
    )
    if(is.na(sigma)) {
      return(list(x = NA, y = NA, sigma = NA))
    }
    return(list(x = x_susie / sigma, y = y_susie / sigma, sigma = sigma))  # , mod = mod))
  } else {
    sigma = tryCatch(
      {
        mod = glmnet::cv.glmnet(x, y, standardize = F, intercept = intercept)
        beta = as.vector(coef(mod$glmnet.fit)[, which.min(mod$cvm)])
        x_pred = cbind(rep(1, nrow(x)), x)
        sigma = sqrt(mean((y - x_pred %*% beta)^2))
        sigma
      },
      error = function(cond) {
        message('cv.glmnet failed with error:')
        message(cond)
        message(' discard the fit')
        # message('use sigma = sqrt(mean(y_bar^2))')
        # mod = lm(y ~ 1)
        # y_bar = residuals(mod)
        # sigma = sqrt(mean(y_bar^2))
        sigma = NA
        return(sigma)
      } 
    )
    if(is.na(sigma)) {
      return(list(x = NA, y = NA, sigma = NA))
    }
    
    if(intercept == T) {
      y = scale(y, T, F)
    }
    return(list(x = x / sigma, y = y / sigma, sigma = sigma))  # , mod = mod))
  }
}

#' @importFrom susieR susie
run_susie_default = function(x, y, ...) {
  mod = tryCatch({
    susie(x, y, ...)
  }, error = function(e) {
    vars = data.frame(variable = 1 : ncol(x), variable_prob = 0, cs = -1)
    list(vars = vars, cs = NULL)
  })
  mod
}

mydiag = function(vals) {
  diag(matrix(vals)[, 1])
}

#' Wrapper for fitting elastic net model with nested cross-validation
#'
#' Build additive linear model to predict y from x using elastic net.
#' It fits the use case when the number of predictors is much larger than
#' the number of observations.
#'
#'
#'
#' @param x x (N by P, usually P >> N)
#' @param y y (N by 1)
#' @param nfold fold in cross-valiation to pick hyperparameter lambda (default = 5)
#' @param lambda_seq specify the sequence of hyperparameter lambda to search.
#' If don't have such lambda_seq, set it to NULL,
#' and the function will generate one internally (default = NULL)
#' @param alpha alpha parameter in glmnet (defualt = 0.5)
#' @param n_times how many rounds of partitions you'd like to perform to estimate
#' CV MSE for each lambda (it is specific for nested CV procedure; default = 3)
#' @param ... additional arguments passed to cv.glmnet
#'
#' @seealso \code{\link{do_elastic_net}} for details on nested CV procedure.
#'
fit_glmnet_with_cv = function(x, y, nfold = 5, lambda_seq = NULL, alpha = 0.5, n_times = 3, ...) {
  # set.seed(1)
  n_k_folds = nfold
  n_samples = nrow(x)

  ## adapted from https://github.com/hakyimlab/PredictDB_Pipeline_GTEx_v7/blob/master/model_training/scripts/gtex_v7_elasticnet.R
  cv_fold_ids <- matrix(nrow = n_samples, ncol = n_times)
  for (j in 1:n_times)
    cv_fold_ids[,j] = sample(1 : n_k_folds, n_samples, replace = TRUE)
  elnet_out = do_elastic_net(x, y, n_k_folds, cv_fold_ids, alpha, ...)
  ## END

  best_lam_ind = which.min(elnet_out$avg_cvm[-1]) + 1
  best_fit = elnet_out$cv_fit$glmnet.fit
  beta = best_fit$beta[, best_lam_ind]
  beta = c(best_fit$a0[best_lam_ind], beta)
  mod = list(beta = beta, fit = elnet_out, best_idx = best_lam_ind)
  mod
}

#' Adapted from https://github.com/hakyimlab/PredictDB_Pipeline_GTEx_v7/blob/master/model_training/scripts/gtex_v7_elasticnet.R
#'
#' Fit y ~ x + 1 with elastic net using nested cross-validation.
#' By "nested", we mean determine MSE using multiple partitions (K)
#' instead of only one in vanilla cross-validation.
#' The procedure is:
#' 1. determine lambda sequence and fit y ~ x at different lambda values;
#' 2. take in K copies of the partition of the input (one of the input, cv_fold_ids).
#' For instance, 5-fold partition if you'd like 5-fold CV;
#' 3. compute CV MSE at each lambda value for each fold in each partition;
#' 4. compute CV MSE for each lambda as the mean of CV MSE of all partitions;
#' 5. select the best model with the smallest CV MSE at that lambda.
#'
#' @param cis_gt x (genotype matrix, sample x variant)
#' @param expr_adj y (expression level vector)
#' @param n_folds cross-validation fold
#' @param cv_fold_ids K partitions with N-fold (sample x K, where each column is partition of samples into N parts represented by the index of each part)
#' @param alpha alpha parameter in glmnet
#' @param ... additional arguments passed to cv.glmnet
#'
#' @importFrom glmnet cv.glmnet
do_elastic_net <- function(cis_gt, expr_adj, n_folds, cv_fold_ids, alpha, ...) {
  #tryCatch({
  cis_gt = as.matrix(cis_gt)
  n_times = ncol(cv_fold_ids)
  fit = cv.glmnet(cis_gt, expr_adj, nfolds = n_folds, alpha = alpha, keep = TRUE, type.measure = 'mse', foldid = cv_fold_ids[,1], parallel = FALSE, ...)
  lambda_seq = fit$lambda
  cvms = matrix(nrow = length(lambda_seq), ncol = n_times)
  fits = list()
  fits[[1]] = fit
  cvms = matrix(nrow = 100, ncol = n_times)
  cvms[1:length(fit$cvm),1] = fit$cvm
  for (i in 2:(n_times)) {
    fit = cv.glmnet(cis_gt, expr_adj, lambda = lambda_seq, nfolds = n_folds, alpha = alpha, keep = FALSE, foldid = cv_fold_ids[,i], parallel = FALSE, ...)
    fits[[i]] <- fit
    cvms[1 : length(fit$cvm), i] <- fit$cvm
  }
  avg_cvm = rowMeans(cvms)
  best_lam_ind = which.min(avg_cvm)
  best_lambda = lambda_seq[best_lam_ind]
  out = list(cv_fit = fits[[1]], min_avg_cvm = min(avg_cvm, na.rm = T), best_lam_ind = best_lam_ind, best_lambda = best_lambda, lambda_seq = lambda_seq, avg_cvm = avg_cvm)
  out
}
## END

#' @importFrom emma emma.MLE.noX emma.MLE
calc_var = function(y, x, fixed_effect = NULL) {
  kinship = x %*% t(x)
  if(is.null(fixed_effect)) {
    res = emma.MLE.noX(y, kinship)
    return(res$ve)
  } else {
    res = emma.MLE(y, fixed_effect, kinship)
    return(res$ve)
  }
}
