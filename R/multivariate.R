## Functions for biomarker discovery (e.g. for drug sensitivity in cell lines)
## Here: multivariate analyses (considering features together)


#' Combined performance metrics for classification
#'
#' Package `caret` offers different model performance metrics for two-class classification.
#' Use all of them together.
#'
#' @inheritParams caret::defaultSummary
combined.summary <- function(data, lev=NULL, model=NULL) {
  def.summary <- caret::defaultSummary(data, lev, model)
  c2.summary <- caret::twoClassSummary(data, lev, model)
  pr.summary <- caret::prSummary(data, lev, model)
  c(def.summary, c2.summary, pr.summary)
}


#' Build multivariate models (random forest), including feature selection
#'
#' This function is based on [caret::rfe()].
#' It uses cross-validation to train and evaluate the performance of random forest models.
#' First using the full set of features, which generates estimates of feature importance.
#' Subsequently more models are trained and evaluated using only the top `N` ranked features, where `N` is given by the `n.features` parameter.
#' The best model/feature set is the one that achieves near-optimal performance using the fewest features.
#'
#' @param responses Named vector of responses (numeric or factor)
#' @param data Named data frame of feature data (e.g. gene expression)
#' @param n.features Numbers of features to consider (vector of integers)
#' @param n.trees Number of trees to include in random forest models
#' @param fit.func Custom model fitting function
#' @param rank.func Custom feature ranking function
#' @param metric Prediction performance metric (default: "RMSE" for regression, "Accuracy" for classification)
#' @param quick Use 5x5-fold cross-validation instead of 10x10-fold (default)
#' @param extended.summary Collect extended performance metrics (classification only)
#' @param ... Additional parameters passed to [caret::rfe()]
#' @return RFE object
#' @export
multivariate.analysis <- function(responses, data, n.features=c(5, 10, 20, 50, 100), n.trees=2000,
                                  fit.func=NULL, rank.func=NULL, metric=NULL,
                                  quick=FALSE, extended.summary=FALSE, ...) {
  ## match samples:
  ## TODO: move to separate function to avoid copies?
  if ((length(responses) != nrow(data)) || !all(names(responses) == rownames(data))) {
    inter <- intersect.samples(responses, data)
    responses <- inter$responses
    data <- inter$data
  }

  rf.funcs <- caret::rfFuncs # default RFE params/functions for random forest models
  ## adjust some of the defaults depending on use case:
  rf.funcs$selectSize <- caret::pickSizeTolerance
  if (!is.null(fit.func))
    rf.funcs$fit <- fit.func # custom model fitting function
  if (!is.null(rank.func))
    rf.funcs$rank <- rank.func # custom feature ranking function
  if (extended.summary) {
    rf.funcs$summary <- combined.summary # see above
  }
  ## 5x5-fold or 10x10-fold cross-validation?
  cv.num <- if (quick) 5 else 10
  rf.ctrl <- caret::rfeControl(rf.funcs, method="repeatedcv", number=cv.num, repeats=cv.num,
                               verbose=TRUE, returnResamp="all", saveDetails=TRUE)
  if (is.null(metric))
    metric <- ifelse(is.factor(responses), "Accuracy", "RMSE") # caret default
  ## call the caret RFE function - note the extra "ntree" parameter that gets
  ## passed to the "randomForest" function:
  caret::rfe(data, responses, n.features, metric=metric, rfeControl=rf.ctrl,
             ntree=n.trees, ...)
}


#' Plot predictions vs. observations for a multivariate regression result
#'
#' @param rfe RFE object (result from [multivariate.analysis()])
#' @param n.feat Number of features to use (default: optimum as selected by caret)
#' @param label Label for responses (type of data)
#' @param ... Additional parameters passed to [plot()], e.g. `asp=1`
#' @export
plot.predictions <- function(rfe, n.feat=NULL, label="", ...) {
  if (is.null(n.feat)) n.feat <- rfe$bestSubset # or 'optsize' (why do both exist?)
  pred <- rfe$pred[rfe$pred$Variables == n.feat, ]
  avg.pred <- tapply(pred$pred, pred$rowIndex, mean)
  row.inds <- sort(unique(pred$rowIndex))
  obs <- pred[match(row.inds, pred$rowIndex), "obs"]
  cor.reg <- cor(avg.pred, obs)
  plot(pred$obs, pred$pred, pch=20, col="grey", xlab=paste("Observed", label),
       ylab=paste("Predicted", label), main="Prediction performance", ...)
  points(obs, avg.pred, pch=21, bg="blue")
  mtext(paste0("r =", format(cor.reg, digits=2), ", n = ", length(rfe$fit$y)))
  reps <- paste0("CV repeats (x", rfe$control$repeats, ")")
  legend("topleft", c(reps, "CV average"), pch=c(20, 21), col=c("grey", "black"),
         pt.bg=c(NA, "blue"))
}


#' Extract feature importance from a multivariate analysis result
#'
#' @param rfe RFE object (result from [multivariate.analysis()])
#' @param return.means Return vector of mean importance over cross-validation folds (default) or full matrix?
#' @export
get.feature.ranking <- function(rfe, return.means=TRUE) {
  imp <- rfe$variables
  ## importance was only measured once (at the start) in each CV fold:
  imp <- imp[imp$Variables == max(imp$Variables), ]
  means <- tapply(imp$Overall, imp$var, mean)
  means <- sort(means, decreasing=TRUE)
  if (return.means) return(means)
  ## reformat into a matrix:
  wide <- pivot_wider(imp, id_cols=!Variables, names_from=Resample, values_from=Overall)
  imp <- as.matrix(wide[, -1])
  rownames(imp) <- wide[[1]]
  imp[names(means), ]
}


#' Plot feature ranking for a multivariate analysis result
#'
#' @param rfe RFE object (result from [multivariate.analysis()])
#' @param n.feat Number of features to show (default: optimum as selected by caret)
#' @export
plot.feature.ranking <- function(rfe, n.feat=NULL) {
  ranking <- get.feature.ranking(rfe, FALSE)
  if (is.null(n.feat)) {
    part <- ranking[rfe$optVariables, ]
  } else {
    part <- ranking[1:n.feat, ]
  }
  yrange <- c(min(0, min(part)), max(part)) # ensure y axis includes zero
  boxplot(t(part), las=3, ylab="Feature importance", main="Multivariate feature ranking", ylim=yrange)
  grid(nx=NA, ny=NULL)
  subtitle <- paste0("Top ", nrow(part), " of ", nrow(ranking), " features (",
                     rfe$control$repeats, "x", rfe$control$number, "-fold cross-validation)")
  mtext(subtitle)
}
