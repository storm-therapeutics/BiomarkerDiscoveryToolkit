## Functions for biomarker discovery (e.g. for drug sensitivity in cell lines)
## Here: univariate analyses (considering e.g. one gene at a time)

#' Filter numeric features by applying a minimum threshold
#'
#' Features (i.e. columns in `data`) with values below the threshold in too many samples (rows) are removed.
#'
#' @param data Data matrix to be filtered
#' @param threshold Threshold value
#' @param fraction Fraction of samples that must pass the treshold (0: one sample, 1: all samples)
#' @return Data matrix after filtering
#' @export filter.threshold
filter.threshold <- function(data, threshold=1, fraction=1) {
  min.samples <- ceiling(nrow(data) * fraction)
  min.samples <- min(max(1, min.samples), nrow(data)) # must be between 1 and number of samples
  pass.count <- apply(data, 2, function(col) sum(col >= threshold, na.rm=TRUE))
  data[, pass.count >= min.samples]
}


#' Filter a response vector and data matrix to contain common samples
#'
#' Only samples that are common between `names(responses)` and `rownames(data)` are kept.
#' Features with missing values (`NA`) in all samples are removed as well.
#'
#' @param responses Named vector
#' @param data Named matrix
#' @return List containing responses and data after filtering
intersect.samples <- function(responses, data) {
  inter <- intersect(names(responses), rownames(data))
  if (length(inter) == 0) stop("No overlapping samples between `responses` and `data`")
  data <- data[inter, , drop=FALSE]
  all.na <- apply(data, 2, function(col) all(is.na(col)))
  if (any(all.na)) {
    warning("Removing ", sum(all.na), " feature(s) with missing values (NA) in all samples")
  }
  list(responses=responses[inter], data=data[, !all.na])
}


#' Test for univariate associations using correlations
#'
#' Every column of `data` is compared to `responses`.
#' Strength of association is measured by Spearman and Pearson correlation coefficients, and significance is estimated based on null distributions.
#' If `out.prefix` is not empty, two output files are created by appending ".csv" and ".pdf" to `out.prefix`.
#' The data frame containing results will be stored as a CSV file.
#' If `plot` is TRUE, a PDF file will be generated with the following plots:
#' Distributions of the Spearman/Pearson correlations and null distributions used to estimate p-values (via [plot.correlation.densities()]); scatter plots for the features with highest positive/negative correlations (`plot.cors` features in each case, via [plot.correlations()]).
#'
#' @param responses Named numeric vector of responses (e.g. drug effects, genetic dependencies)
#' @param data Named numeric matrix of feature data (e.g. gene expression)
#' @param out.prefix Path and filename prefix for output files (extensions will be appended)
#' @param sample.names Subset of sample names to use
#' @param n.null Number of repeats for null distribution
#' @param plot Create PDF with plots?
#' @param plot.cors Number of top hits (pos./neg. correlations) to plot
#' @param plot.rows Number of rows for arranging plots of top hits
#' @param plot.xlab X axis label for [plot.correlations()] plots
#' @param plot.ylab Y axis label for [plot.correlations()] plots
#' @return Numeric matrix of correlation coefficients and p-values
#' @export
correlation.analysis <- function(responses, data, out.prefix="", sample.names=NULL,
                                 n.null=10, plot=TRUE, plot.cors=15, plot.rows=3,
                                 plot.xlab="expression", plot.ylab="response") {
  ## match samples:
  ## TODO: move to separate function to avoid copies
  if (!is.null(sample.names)) {
    data <- data[rownames(data) %in% sample.names, ]
  }
  if ((length(responses) != nrow(data)) || !all(names(responses) == rownames(data))) {
    inter <- intersect.samples(responses, data)
    responses <- inter$responses
    data <- inter$data
  }

  cor.pvs <- correlations.with.pvalues(responses, data, methods=c("spearman", "pearson"),
                                       n.null=n.null, return.null=plot)
  ## beware: subsetting loses attributes! - stick them back on afterwards (but don't reset names!):
  attribs <- attributes(cor.pvs)[c("null.distributions", "null.sd")]
  cor.pvs <- cor.pvs[order(abs(cor.pvs[, "cor.spearman"]), abs(cor.pvs[, "cor.pearson"]), decreasing=TRUE), ]
  attributes(cor.pvs)[c("null.distributions", "null.sd")] <- attribs

  if (out.prefix != "") { # write output file(s)
    utils::write.csv(cor.pvs, paste0(out.prefix, ".csv"))

    if (plot) { # plot to PDF file
      grDevices::pdf(paste0(out.prefix, ".pdf"), paper="a4r", height=0, width=0)
      on.exit(grDevices::dev.off())
      plot.correlation.densities(cor.pvs, "spearman")
      plot.correlation.densities(cor.pvs, "pearson")
      if (plot.cors > 0) {
        is.negative <- cor.pvs[, "cor.spearman"] < 0
        plot.correlations(cor.pvs[!is.negative, ], responses, data, 1:plot.cors, plot.rows,
                          "Top positive correlations", xlab=plot.xlab, ylab=plot.ylab)
        plot.correlations(cor.pvs[is.negative, ], responses, data, 1:plot.cors, plot.rows,
                          "Top negative correlations", xlab=plot.xlab, ylab=plot.ylab)
      }
    }
  }

  invisible(cor.pvs)
}


#' Compare data for groups of responses using a statistical test
#'
#' For every column in `data`, compare the distributions of values based on the groups (levels) in `responses`.
#' Use the statistical test selected by parameter `test` to calculate p-values.
#' By default the test is [wilcox.test()] when `responses` has two levels (i.e. two groups), otherwise [kruskal.test()].
#' If a custom test function is supplied for `test`, it has to support the formula interface.
#'
#' @param responses Named factor of categorical responses (e.g. responders/non-responders)
#' @param data Named numeric matrix of feature data (e.g. gene expression)
#' @param out.prefix If set: Path and filename prefix for output files (extensions will be appended). Default: NULL.
#' @param sample.names Subset of sample names to use
#' @param test Statistical test to use (default: [wilcox.test()] for two groups, otherwise [kruskal.test()])
#' @param p.adj.method Method for multiple testing correction using [p.adjust()]. Default: BH
#' @param plot Create PDF with plots?
#' @param plot.hits Number of top hits to plot
#' @param plot.rows Number of rows for arranging plots of top hits
#' @param plot.xlab X axis label for plots
#' @param plot.ylab Y axis label for plots
#' @param ... Additional parameters passed to `test`
#' @return Data frame of results (p-values and basic statistics)
#' @export
#'
#' @examples
#'\dontrun{
#'responses=readRDS( "data/responses_group.rds")
#'data=readRDS("data/data_group.rds")
#'check=group.analysis(data=data,responses=responses)
#'}
group.analysis <- function(responses, data, out.prefix=NULL, sample.names=NULL, test=NULL, p.adj.method="BH",
                           plot=(!is.null(out.prefix)), plot.hits=15, plot.rows=3,
                           plot.xlab="", plot.ylab="expression", ...) {
  ## match samples:
  ## TODO: move to separate function to avoid copies
  if (!is.null(sample.names)) {
    data <- data[rownames(data) %in% sample.names, ]
  }
  if ((length(responses) != nrow(data)) || !all(names(responses) == rownames(data))) {
    inter <- intersect.samples(responses, data)
    responses <- inter$responses
    data <- inter$data
  }

  if (is.null(test)) test <- ifelse(length(levels(responses)) > 2, kruskal.test, wilcox.test)

  fullFrame <- data.frame(response=responses, data, check.names=FALSE)

  ## run statistical test for every column of `data`, comparing distributions according to grouping in `responses`
  results <- lapply(colnames(data), function(feature) {
    f <- formula(paste0("`", feature, "` ~ response"))
    stats <- test(f, data=fullFrame, ...)
    tab <- table(responses[!is.na(data[, feature])])
    names(tab) <- paste0("n.", names(tab))
    means <- tapply(data[, feature], responses, mean, na.rm=TRUE)
    names(means) <- paste0("mean.", names(means))
    data.frame(as.list(tab), as.list(means), statistic=stats$statistic, p.value=stats$p.value, check.names=FALSE)
  })

  thelper <- do.call(rbind, results)
  rownames(thelper) <- colnames(data)

  ## apply multiple testing correction
  thelper$p.adj <- p.adjust(thelper$p.value, method=p.adj.method)
  sortedFrame <- thelper[order(thelper$p.value), ]

  ## write output file if out.prefix is set (CSV)
  if (!is.null(out.prefix))
    utils::write.csv(sortedFrame, file=file.path(out.prefix, ".csv"))

  ## plot top results (depending on `plot` and other plot parameters)
  ## TODO: write plot to PDF
  if (plot)
  {
    plot.groups(stats, responses, data, 1:plot.hits, plot.rows, "Top hits by p-value", plot.xlab, plot.ylab)
    #library(gridExtra)

    #ggsave(filename = file.path("test.pdf"), plot = marrangeGrob(plots, nrow=plot.rows ,byrow=TRUE), width = 20, height = 20)
  }
  invisible(sortedFrame)
}


#' Compare groups of responses defined by mutational status using a statistical test
#'
#' For every column in `data`, compare the distributions of response values based on the groups in the column.
#' Use the statistical test selected by parameter `test` to calculate p-values.
#'
#' @param responses Named numeric vector of responses (e.g. drug effects, genetic dependencies)
#' @param data Named binary matrix of gene mutation data (e.g. damaging, likely LoF)
#' @param out.prefix Path and filename prefix for output files (extensions will be appended)
#' @param sample.names Subset of sample names to use
#' @param test Statistical test to use (default: [wilcox.test()])
#' @param min.samples Minimum number of samples in each group for testing
#' @param plot Create PDF with plots?
#' @param plot.hits Number of top hits to plot
#' @param plot.xlab X axis label for plots
#' @param plot.ylab Y axis label for plots
#' @return Data frame of results (p-values and basic statistics)
#' @export
mutation.analysis <- function(responses, data, out.prefix="", sample.names=NULL, test=wilcox.test, min.samples=5, plot=TRUE, plot.hits=10, plot.xlab="mutation", plot.ylab="response") {

  if (!is.null(sample.names)) {
    responses <- responses[sample.names]
  }

  # filter cell models not present in the response and select genes mutated and non-mutated in minimum 10 (5+5) cell models
  response.data <- merge(as.data.frame(responses), data, by=0)
  response.data <- response.data %>%
    drop_na(responses) %>%
    pivot_longer(c(3:ncol(response.data)), names_to = "Gene", values_to = "mutational_status") %>%
    group_by(Gene) %>%
    dplyr::filter(sum(mutational_status == TRUE) >= min.samples & sum(mutational_status == FALSE) >= min.samples)

  # perform selected test by Gene
  if (nrow(response.data) <= 2 * min.samples) stop("Sampling size too small for testing.")
  test_res <- response.data %>%
    do(t = test(responses ~ mutational_status, data=., paired=FALSE)) %>%
    summarise(across(Gene), p.value = t$p.value)
  # calculate medians and sampling sizes for each gene and mutational status
  other_stats <- response.data %>%
    group_by(Gene, mutational_status) %>%
    mutate(median=median(responses)) %>%
    mutate(n=n()) %>%
    select(Gene, mutational_status, n, median) %>%
    pivot_wider(names_from = mutational_status, values_from = c(n, median), values_fn = unique)
  # join results, order by p-value and export to csv
  response_vs_mutation_status <- as.data.frame(inner_join(other_stats, test_res, by="Gene") %>% arrange(p.value))
  colnames(response_vs_mutation_status) <- gsub("_", "\\.", colnames(response_vs_mutation_status))
  rownames(response_vs_mutation_status) <- response_vs_mutation_status$Gene
  response_vs_mutation_status$Gene <- NULL
  utils::write.csv(response_vs_mutation_status, file = paste0(out.prefix, ".csv"))

  if (plot) {
    top.results <- response_vs_mutation_status[c(1:plot.hits),] %>% mutate(Gene = rownames(response_vs_mutation_status)[c(1:plot.hits)])
    p <- response.data %>%
      dplyr::filter(Gene %in% rownames(top.results)) %>%
      left_join(top.results, by="Gene") %>%
      ggplot(aes(y = responses, x = mutational_status)) +
      geom_boxplot(outlier.shape = NA, width=0.5) +
      geom_jitter(width = 0.05, size = 0.1, alpha = 0.5) +
      xlab(plot.xlab) +
      ylab(plot.ylab) +
      facet_wrap(~ factor(Gene, levels = top.results$Gene) + paste("p =", format(p.value, digits=3)), scales = "free") +
      theme_bw()
    ggsave(plot = p, filename = paste0(out.prefix, ".pdf"), width = 210, height = 297, units = "mm")
  }

  invisible(response_vs_mutation_status)
}


#' Test for univariate associations with survival time using Cox regression
#'
#' For every column of `data`, a univariate Cox proportional hazard model ([coxph()]) is built to predict `responses`.
#' Models are ranked according to the significance of the coefficient (feature from `data`).
#'
#' @param responses Named vector of survival times (`survival::Surv()` object)
#' @param data Named numeric matrix of feature data (e.g. gene expression)
#' @param out.path Path to PDF output file (if `plot=TRUE`)
#' @param sample.names Subset of sample names to use
#' @param plot Create PDF with Kaplan-Meier plots?
#' @param plot.hits Number of top hits to plot
#' @param ... Further parameters passed to [plot.cox.pred()], e.g. `xlab`/`ylab`
#' @return List of Cox proportional hazard models for each feature, ordered by p-value
#' @export
survival.analysis <- function(responses, data, out.path="survival_analysis.pdf", sample.names=NULL,
                              plot=TRUE, plot.hits=10, ...) {
  ## match samples:
  ## TODO: move to separate function to avoid copies
  if (!is.null(sample.names)) {
    data <- data[rownames(data) %in% sample.names, ]
  }
  if ((length(responses) != nrow(data)) || !all(names(responses) == rownames(data))) {
    inter <- intersect.samples(responses, data)
    responses <- inter$responses
    data <- inter$data
  }

  merged <- data.frame(response=responses, data, check.names=FALSE)
  models <- lapply(colnames(data), function(feature)
    survival::coxph(as.formula(paste0("response ~ `", feature, "`")), merged))
  names(models) <- colnames(data)
  pvalues <- sapply(models, function(m) coef(summary(m))[1, 5])
  ord <- order(pvalues)

  if (plot) {
    base.curve <- survival::survfit(responses ~ 1)
    grDevices::pdf(out.path)
    for (i in 1:min(length(models), plot.hits)) {
      plot.cox.pred(models[[ord[i]]], merged, base.curve, ...)
    }
    grDevices::dev.off()
  }

  models[ord]
}
