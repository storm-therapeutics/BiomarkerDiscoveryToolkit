## Functions for biomarker discovery (e.g. for drug sensitivity in cell lines)
## Here: utility functions for correlations

#' Compute correlation coefficients between a vector and all rows/columns of a matrix
#'
#' Convenience function for calculating multiple correlation coefficients for a dataset.
#'
#' The length of `vec` must be compatible with the number of rows or columns in `mat` (depending on `margin`).
#' Multiple correlation coefficients per row/column can be calculated if multiple strings are passed to `methods`.
#'
#' @param vec Numeric vector
#' @param mat Numeric matrix
#' @param margin Whether to calculate correlations per row (1) or column (2) of `mat`
#' @param methods Which correlation coefficient(s) to calculate (see `method` argument of [cor()])
#' @return Matrix of correlation coefficients with one column per method used
#' @export compute.correlations
compute.correlations <- function(vec, mat, margin=2, methods=c("pearson", "spearman")) {
  res <- lapply(methods, function(m)
    apply(mat, margin, function(x) cor(x, vec, use="complete.obs", method=m)))
  res <- do.call(cbind, res)
  colnames(res) <- methods
  res
}


#' Generate null distributions of correlation coefficients by permutation
#'
#' This function repeatedly calls [compute.correlations()] after shuffling its first argument (`vec`).
#'
#' @param n Number of repeats
#' @param vec Numeric vector to permutate
#' @param ... Further parameters passed to [compute.correlations()]
#' @return List of length `n` with results of [compute.correlations()] after shuffling `vec`
get.null.correlations <- function(n, vec, ...) {
  ## TODO: parallelise
  lapply(1:n, function(i) {
    message(i, " ", appendLF=FALSE)
    shuffled <- sample(vec)
    compute.correlations(shuffled, ...)
  })
}


#' Estimate p-values based on a normal distribution
#'
#' Intended for estimating the significance of correlation coefficients, this function assumes both `x` and `-x` are equally "extreme" outcomes, hence the p-value is for a two-sided test.
#'
#' @param x Test statistic
#' @param sd Standard deviation of the normal distribution
#' @return P(|X| > |x|) for X ~ N(0, `sd`)
estimate.pvalue <- function(x, sd) {
  2 * (1 - pnorm(abs(x), sd=sd))
}


#' Calculate correlation coefficients and estimate their significance
#'
#' This function runs [compute.correlations()] to calculate correlation coefficients, then uses [get.null.correlations()] to generate corresponding null distributions and estimate p-values.
#'
#' In contrast to [cor.test()], the test is not for "any (non-zero) association", but for "as or more extreme association" (than observed).
#' If `return.null` is TRUE, the null distribution samples are returned as attribute "null.distributions" of the result, and their standard deviations (for each method used) as attribute "null.sd".
#'
#' @param vec Numeric vector
#' @param mat Numeric matrix
#' @param ... Further arguments passed to [compute.correlations()] and [get.null.correlations()]
#' @param n.null Number of repeats for the null distribution
#' @param return.null If `TRUE`, return details about the null distributions as attributes of the result
#'
#' @return Data frame with correlation coefficients and their p-values (in separate columns)
#' @export
correlations.with.pvalues <- function(vec, mat, ..., n.null=10, return.null=FALSE) {
  cors <- compute.correlations(vec, mat, ...)
  methods <- colnames(cors)
  if (n.null == 0) { # skip null distributions/p-values
    colnames(cors) <- paste0("cor.", methods)
    return(cors)
  }

  message("Generating null distributions: ", appendLF=FALSE)
  cors.null <- get.null.correlations(n.null, vec, mat, ...)
  ## standard deviations of null distributions:
  sd.null <- sapply(methods, function(method) {
    sds <- sapply(cors.null, function(x) sd(x[, method], na.rm=TRUE))
    mean(sds)
  })
  ## estimate p-values:
  cor.pvs <- sapply(methods, function(method) {
    pvalues <- estimate.pvalue(cors[, method], sd.null[method])
    pvalues
  })
  colnames(cor.pvs) <- paste0("p.", methods)
  colnames(cors) <- paste0("cor.", methods)
  res <- cbind(cors, cor.pvs)
  if (return.null) {
    attr(res, "null.distributions") <- cors.null
    attr(res, "null.sd") <- sd.null
  }
  res
}


#' Plot densities of correlation coefficients and accompanying null distributions
#'
#' @param cor.res Output of [correlations.with.pvalues()] with `return.null=TRUE`
#' @param method Results from which method to plot, if multiple were used (default: first)
#' @export plot.correlation.densities
plot.correlation.densities <- function(cor.res, method=NULL) {
  if (is.null(attr(cor.res, "null.distributions")))
    stop("No null distribution data available")
  if (is.null(method)) {
    colname <- colnames(cor.res)[1]
    method <- sub("^cor\\.", "", colname)
  } else {
    colname <- paste0("cor.", method)
  }
  real.cors <- cor.res[, colname]
  ## fit densities:
  real.density <- density(real.cors, na.rm=TRUE)
  null.densities <- lapply(attr(cor.res, "null.distributions"), function(x) density(x[, method], na.rm=TRUE))
  ## plot real and null distributions:
  ymax <- max(real.density$y, sapply(null.densities, function(d) max(d$y)))
  plot(real.density, col="red", lwd=2, xlim=c(-1, 1), ylim=c(0, ymax), xlab="correlation coefficient",
       main=paste("Distribution of", stringr::str_to_title(method), "correlations"))
  for (d in null.densities) {
    lines(d, col="grey")
  }
  ## plot normal approximation to null distributions:
  sd.null <- attr(cor.res, "null.sd")
  curve(dnorm(x, sd=sd.null), add=TRUE, col="blue")
  labels <- c("real data", paste0("random permutation (x", length(null.densities), ")"),
              paste0("normal approx. (sd=", format(sd.null[method], digits=3), ")"))
  legend("topleft", labels, lty=1, col=c("red", "grey", "blue"))
  ## show real correlation coefficients on x axis:
  points(real.cors, rep(-0.2, length(real.cors)), pch="|", cex=0.5)
  abline(v=0, lty=3)
}
