#===========================================#
#===========================================#
#' Pre-processing Step 
#'
#' Before applying GeneSurrounder, the observed and resampled 
#' differential expression of the genes should be calculated 
#'
#' @param expr A matrix (genes by samples) of expression values.
#' @param classLabels A factor with levels corresponding to class labels.
#' @param numResamples defaults to 1000. The number of resamples when calculating resampled differential expression. 
#' @importFrom limma treat lmFit
#' @export
# Calc gene level statistics & a null set of gene level stats (shuffle phenotype labels)
calcGeneTStats <- function(expr, classLabels, numResamples = 1000){
  # Calc gene level statistics
  # This code is being written using CurOvGradeKEGGnets[[2]]
  # 
  # Args:
  # expr: is a matrix of genes by samples
  # classLabels: is a vector of class labels (e.g. high vs low)
  # numrResamples: number of times the phenotype labels are shuffled
  # Returns:
  # observedStats: a vector of observed moderated t-statistics
  # permStats: a matrix of resampled moderated t statistics (resamplings are rows)

  # cf. d715_timecourse_contrasts, network_review_GSEAhat ? 
  # Should I save the fit so I have the gene p values etc...? 

  desMat <- model.matrix(~factor(classLabels))
  # treat is a limma function
  # Given a microarray linear model fit, compute moderated t-statistic, etc
  fit <- limma::treat(limma::lmFit(expr, desMat))
  observedStats <- fit$t[, 2]


  permStats <- sapply(1:numResamples, function(resampleLoopIndex) {
    # Shuffle the phenotype labels
    permLabels <- sample(classLabels, replace = FALSE)

    #Refit and recalculcate gene level statistics using permLabels
    permDesMat <- model.matrix(~factor(permLabels))
    permFit <- limma::treat(limma::lmFit(expr, permDesMat))
    #print(head(permFit$t[, 2]))
    return(permFit$t[, 2])
  })

  # Transpose permStats so the rows are resamplings
  permStats <- t(permStats)

  # List and return
  geneTStats <- list(observed = observedStats, resampled = permStats)

  return(geneTStats)	
}

#' Calculate sleuth observed and resampled statistics
#'
#' This method takes a \code{sleuth} object and calculates
#' observed as well as resampled statistics.
#'
#' @param obj the sleuth object. It must have 'whichModel' already fit
#'   using \code{sleuth_fit}, as well as the null model if the testType is
#'   'lrt'.
#' @param whichBeta this is required for both test types to indicate which column of the design
#'   matrix to permute. For "wt", this is the same as \code{whichTest}.
#' @param whichModel if "lrt" is the test type, this the alternative model; if "wt" is
#'   the test type, this is the model to use.
#' @param whichTest if "lrt", this is the test in the format of "[null model]:[alternative model]";
#'   if "wt", this is the same as whichBeta, the actual column in the design matrix to test
#' @param testType either "lrt" or "wt" to indicate which test, Likelihood Ratio Test or
#'   Wald Test, to use to calculate statistics
#' @param numResamples the number of times the samples should be resampled.
#'   Note that this method only uses unique resamplings, so if this is \code{NULL} or
#'   if the specified number of resamplings is equal to or more than the total unique
#'   combinations, this method will calculate statistics for all unique resamplings.
#' @param allPerms boolean for whether all permutations should be
#'   generated
#' @param numCores the number of cores to use for parallel computation of the resampled stats
#' @param ... these are additional parameters to pass to \code{link{sleuth_fit}}.
#'
#' @return a list containing two items:
#'   \itemize{
#'     \item{observed}{this is vector of the observed statistics; the likelihood ratio for "lrt", and
#'        wald statistic for "wt"}
#'     \item{resampled}{this is a matrix of number-of-features rows and numResamples columns, containing
#'        the resampled statistics}
#'   }
#' @importFrom methods is
#' @importFrom parallel mclapply
#' @importFrom sleuth sleuth_fit sleuth_lrt sleuth_wt
#' @export
calcGeneSleuthStats <- function(obj, whichBeta, whichModel = "full",
                                whichTest = "reduced:full", testType = "lrt",
                                numResamples = NULL, allPerms = TRUE,
                                numCores = 1L, ...) {
  stopifnot(is(obj, "sleuth"))

  if (missing(whichBeta)) {
    stop("'whichBeta' is missing. It must be specified to know which labels to permute")
  } else if (testType == "wt" && whichTest != whichBeta) {
    stop("for 'testType' 'wt', 'whichTest' and 'whichBeta' must match")
  }

  numCores <- check_cores(numCores)

  observedStats <- get_sleuth_stats(obj, whichTest, whichModel)

  design <- obj$fits[[whichModel]]$design_matrix
  classLabels <- design[, whichBeta]

  permMat <- returnResampleMat(classLabels, numResamples, allPerms = allPerms)

  permStats <- parallel::mclapply(1:nrow(permMat), function(i) {
    permLabels <- permMat[i, ]
    perm_design <- design
    perm_design[, whichBeta] <- permLabels
    perm_obj <- sleuth::sleuth_fit(obj, formula = perm_design, fit_name = whichModel, ...)
    if (testType == "lrt") {
      perm_obj <- sleuth::sleuth_lrt(perm_obj, null_model, whichModel)
      perm_stats <- perm_obj$tests$lrt[[whichTest]]$test_stat
    } else {
      perm_obj <- sleuth::sleuth_wt(perm_obj, whichBeta, whichModel)
      perm_stats <- perm_obj$tests$wt[[whichModel]][[whichTest]]$wald_stat
    }
    perm_stats
  }, mc.cores = numCores)
  names(permStats) <- paste0("resample_", 1:nrow(permMat))
  permStats <- as.matrix(as.data.frame(permStats))

  list(observed = observedStats, resampled = permStats)
}
