% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rlib_meta.R
\name{meta_analyze}
\alias{meta_analyze}
\title{Meta-analyze trcQTL and ascQTL results}
\usage{
meta_analyze(trc, asc, n_cutoff = 15)
}
\arguments{
\item{trc}{list of effect size estimates and standard deviation (e.g. output of matrix_ls_trc)}

\item{asc}{list of effect size estimates and standard deviation (e.g. output of matrix_ls_asc)}

\item{n_cutoff}{threshold to approximate test statistic as z-score.}
}
\value{
a list of summary statistics and p-values obtained from trc, asc, and meta-analysis of the two
        the object within list is named by the method (trc, asc, and meta) and each of which includes
        pval (p-value), stat (test statistic), stat_type (z-value or t-value),
        bhat (effect size estimate), se (standard deviation of effect size),
        method (only for 'meta' object showing if the result is from 'trc', 'asc', or 'meta')
}
\description{
For each variant/gene pair, meta-analyze summary statistics from trcQTL and ascQTL.
If sample size is greater than or equal to n_cutoff,
approximate test statistic as z-score,
otherwise compute p-value from student-t distribution with d.o.f. = 1.
If both trcQTL and ascQTL have no fewer than n_cutoff samples, perform
inverse variance-based meta-analysis.
}
\examples{
meta_analyze(
  trc = list(
    beta_hat = matrix(rnorm(100), ncol = 20),
    beta_se = abs(matrix(rnorm(100), ncol = 20)),
    sample_size = 100
  ),
  asc = list(
    beta_hat = matrix(rnorm(100), ncol = 20),
    beta_se = abs(matrix(rnorm(100), ncol = 20)),
    sample_size = 100
  )
)

}
