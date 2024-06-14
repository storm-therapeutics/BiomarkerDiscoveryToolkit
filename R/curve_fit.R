## Functions for biomarker discovery (e.g. for drug sensitivity in cell lines)
## Here: utility functions for dose-response fits

#' Try to fit a dose-response curve using [dr4pl()], return `NULL` on failure
#'
#' The default parameters of the `dr4pl` call (`method.robust="absolute", method.init="logistic"`) were chosen to improve robustness.
#'
#' @param doses Vector of drug doses/concentrations
#' @param responses Vector or response values
#' @param ... Additional parameters passed to [dr4pl()]
#' @return `dr4pl` (curve fit) object or `NULL`
fit.robust <- function(doses, responses, ...) {
  tryCatch(dr4pl::dr4pl(doses, responses, method.robust="absolute", method.init="logistic", ...),
           error=function(e) NULL)
  ## TODO: raise a warning?
  ## TODO: return dose/response values on failure (for plotting)?
}


#' Read assay data from sheets in an Excel file and fit dose-response curves
#'
#' Sample names (e.g. cell lines) are taken from the names of the sheets.
#' Replicate measurements are expected in separate columns ("wide" format).
#'
#' @param excel.file Path to the Excel file
#' @param sheets Names of sheets to read (default: use all sheets)
#' @param dose.col Column (number) containing dose values
#' @param rep.cols Columns containing replicate response measurements (e.g. absorbance)
#' @param ref.row Reference row for relative responses (default: use row with dose 0)
#' @param ... Additional parameters passed to `dr4pl(method.robust="absolute")`
#' @return List of `dr4pl` curve fits (or `NULL` for failed fits)
#' @export
fit.data.sheets <- function(excel.file, sheets=NULL, dose.col=1, rep.cols=-1, ref.row=NULL, ...) {
  if (is.null(sheets)) sheets <- readxl::excel_sheets(excel.file) # use all sheets in the file

  data.all <- lapply(sheets, function(sheet) {
    ## suppress "New names:" messages from 'read_excel':
    data <- as.data.frame(suppressMessages(readxl::read_excel(excel.file, sheet)))
    ## keep only relevant columns:
    rep.cols <- (1:ncol(data))[rep.cols] # handle default value -1
    data <- data[, c(dose.col, rep.cols)]
    names(data) <- c("dose", paste0("rep", seq_along(rep.cols)))
    data
  })
  names(data.all) <- sheets

  ## transform to relative response (and long format):
  data.rel <- lapply(data.all, function(data) {
    if (is.null(ref.row)) {
      ref.row <- which(data$dose == 0)
    }
    if (length(ref.row) != 1) {
      warning("Could not determine reference row - using absolute response values")
      baseline <- 1
    } else {
      baseline <- median(as.numeric(data[ref.row, -1]), na.rm=TRUE)
      data <- data[-ref.row, ]
    }
    data.long <- as.data.frame(pivot_longer(data, -1, names_to="rep", names_prefix="[^0-9]+", values_to="response"))
    data.long$response <- data.long$response / baseline
    data.long
  })

  lapply(data.rel, function(data) fit.robust(data$dose, data$response, ...))
}


#' Read assay data in "long" format and fit dose-response curves
#'
#' @param data Data frame containing assay data
#' @param name.col Column containing sample (cell line) names
#' @param dose.col Column containing dose values
#' @param response.col Column containing response values
#' @param percent Are response values percentages?
#' @param ... Additional parameters passed to [dr4pl()]
#' @return List of `dr4pl` curve fits (or `NULL` for failed fits)
#' @export
fit.data.long <- function(data, name.col, dose.col, response.col, percent=FALSE, ...) {
  if (percent) data[[response.col]] <- data[[response.col]] / 100
  grouping <- if (name.col == 0) rownames(data) else data[[name.col]]
  parts <- split(data, grouping)
  lapply(parts, function(part) fit.robust(part[[dose.col]], part[[response.col]], ...))
}


#' Compute the area under the curve (AUC) of a dose-response fit
#'
#' Calculates the integral under the dose-response curve within the dose range and normalizes by the dose range (so the result is between 0 and 1).
#' Dose points are expected to be on a logarithmic scale.
#'
#' @param dr4pl.fit Dose-response curve fit ([dr4pl()] result)
#' @param dose.range Dose range (optional - taken from `dr4pl.fit` if missing)
#' @param bounded Whether response values have to be between 0 and 1
#' @return AUC value
#' @export compute.AUC
compute.AUC <- function(dr4pl.fit, dose.range, bounded=TRUE) {
  if (missing(dose.range))
    dose.range <- range(dr4pl.fit$data$Dose)
  response <- function(x) dr4pl::MeanResponse(coef(dr4pl.fit), 2^x) # undo the log
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
#' @export
get.IC50 <- function(dr4pl.fit) {
  ic50 <- unname(dr4pl::IC(dr4pl.fit, 50))
  if ((ic50 > max(dr4pl.fit$data$Dose)) || (ic50 < min(dr4pl.fit$data$Dose))) return(NA)
  ic50
}


#' Plot an annotated dose-response curve
#'
#' @param dr4pl.fit Dose-response curve fit ([dr4pl()] result)
#' @param auc AUC value (default: calculated from `dr4pl.fit`)
#' @param ic50 IC50 value (default: calculated from `dr4pl.fit`)
#' @return Plot with annotations
#' @export plot.curve.fit
plot.curve.fit <- function(dr4pl.fit, auc=NULL, ic50=NULL) {
  if (is.null(auc)) auc <- compute.AUC(dr4pl.fit)
  if (is.null(ic50)) ic50 <- get.IC50(dr4pl.fit)

  plot <- plot(dr4pl.fit) + labs(title="")
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
#' @export plot.curve.fits
plot.curve.fits <- function(path, dr4pl.fits, aucs=NULL, ic50s=NULL) {
  if (is.null(aucs)) aucs <- rep(NULL, length(dr4pl.fits))
  if (is.null(ic50s)) ic50s <- rep(NULL, length(dr4pl.fits))

  grDevices::pdf(path)
  for (i in seq_along(dr4pl.fits)) {
    print(plot.curve.fit(dr4pl.fits[[i]], aucs[i], ic50s[i]) +
          labs(title=names(dr4pl.fits)[i]))
  }
  grDevices::dev.off()
}
