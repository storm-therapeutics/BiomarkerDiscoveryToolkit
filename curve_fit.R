## Functions for biomarker discovery (e.g. for drug sensitivity in cell lines)
## Here: utility functions for dose-response fits

library(dr4pl)

#' Compute the area under the curve (AUC) of a dose-response fit
#'
#' Calculates the integral under the dose-response curve within the dose range and normalizes by the dose range (so the result is between 0 and 1).
#' Dose points are expected to be on a logarithmic scale.
#'
#' @param dr4pl.fit Dose-response curve fit ([dr4pl()] result)
#' @param dose.range Dose range (optional - taken from `dr4pl.fit` if missing)
#' @param bounded Whether response values have to be between 0 and 1
#' @return AUC value
compute.AUC <- function(dr4pl.fit, dose.range, bounded=TRUE) {
  if (missing(dose.range))
    dose.range <- range(dr4pl.fit$data$Dose)
  response <- function(x) MeanResponse(coef(dr4pl.fit), 2^x) # undo the log
  if (bounded)
    curve <- function(x) pmax(pmin(response(x), 1), 0)
  else
    curve <- response
  ## area under the curve, normalised by x axis range:
  tryCatch((integrate(curve, log2(dose.range[1]), log2(dose.range[2]))$value) /
           log2(dose.range[2] / dose.range[1]), error=function(x) NA)
}


#' Return the IC50 of a dose-reponse curve fit
#'
#' If the estimated value is outside the dose range, `NA` is returned.
#'
#' @param dr4pl.fit Dose-response curve fit ([dr4pl()] result)
#' @return IC50 value (or `NA`)
get.IC50 <- function(dr4pl.fit) {
  ic50 <- unname(IC(dr4pl.fit, 50))
  if (ic50 > max(fit$data$Dose)) return(NA)
  ic50
}


#' Plot an annotated dose-response curve
#'
#' @param dr4pl.fit Dose-response curve fit ([dr4pl()] result)
#' @param auc AUC value (default: calculated from `dr4pl.fit`)
#' @param ic50 IC50 value (default: calculated from `dr4pl.fit`)
#' @return Plot with annotations
plot.curve.fit <- function(dr4pl.fit, auc=NULL, ic50=NULL) {
  if (is.null(auc)) auc <- compute.AUC(dr4pl.fit)
  if (is.null(ic50)) ic50 <- get.IC50(dr4pl.fit)

  plot <- plot(dr4pl.fit)
  if (!is.na(ic50))
    plot <- plot + geom_vline(xintercept=ic50, color="red", linetype = "dashed") +
      annotate("text", x=ic50, y=min(dr4pl.fit$data$Response), label="IC50", color="red", hjust=1.1)

  plot + annotate("label", x=max(dr4pl.fit$data$Dose) / 2, y=max(dr4pl.fit$data$Response),
                  vjust=1, size=5, label.padding=unit(0.3, "lines"),
                  label=paste0("IC50 = ", round(ic50, 1), "\nAUC = ", round(auc, 3)))
}


#' Plot multiple dose-response curves to PDF
#'
#' One plot (using [plot.curve.fit()]) per page.
#' The names of `dr4pl.fits` are used as titles for the plots.
#'
#' @param path Path to the output PDF file
#' @param dr4pl.fits Named list of dose-response curve fits ([dr4pl()] results)
#' @param aucs AUC values (default: calculated from `dr4pl.fits`)
#' @param ic50s IC50 values (default: calculated from `dr4pl.fits`)
plot.curve.fits <- function(path, dr4pl.fits, aucs=NULL, ic50s=NULL) {
  if (is.null(aucs)) aucs <- rep(NULL, length(dr4pl.fits))
  if (is.null(ic50s)) ic50s <- rep(NULL, length(dr4pl.fits))

  pdf(path)
  for (i in seq_along(dr4pl.fits)) {
    print(plot.curve.fit(dr4pl.fits[[i]], aucs[i], ic50s[i]) +
          labs(title=names(dr4pl.fits)[i]))
  }
  dev.off()
}
