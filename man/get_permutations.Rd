% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rlib_matrix_ls_with_mask.R
\name{get_permutations}
\alias{get_permutations}
\title{Generate permutation matrix}
\usage{
get_permutations(n, nindiv, seed = 1)
}
\arguments{
\item{n}{number of permutation (as K)}

\item{nindiv}{number of observations to permutate (as N)}

\item{seed}{set.seed(seed) (default = 1)}
}
\value{
N-by-K matrix where each column is permutation of 1:N
}
\description{
Given the number of observations and the number of permutations, generate the permutated indexes (from 1 to the total number of observations) at each column of output. The number of columns of the output equals to the number of permutation desired.
}
\examples{
get_permutations(n = 10, nindiv = 100, seed = 2019)

}
