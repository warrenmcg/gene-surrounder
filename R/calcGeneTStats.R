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
#' @param obj the sleuth object
#' @param numResamples the number of times the samples should be resampled.
#'   Note that this method only uses unique resamplings, so if this is \code{NULL} or
#'   if the specified number of resamplings is equal to or more than the total unique
#'   combinations, this method will calculate statistics for all unique resamplings.
#' @param numCores the number of cores to use for parallel computation of the resampled stats
#' @param testType either "lrt" or "wt" to indicate which test, Likelihood Ratio Test or
#'   Wald Test, to use to calculate statistics
#' @param whichModel if "lrt" is the test type, this the alternative model; if "wt" is
#'   the test type, this is the model to use.
#' @param whichTest if "lrt", this is the test in the format of "[null model]:[alternative model]";
#'   if "wt", this is the same as whichBeta, the actual column in the design matrix to test
#' @param whichBeta this is required for both test types to indicate which column of the design
#'   matrix to permute. For "wt", this is the same as \code{whichTest}.
#' @param ... these are additional parameters to pass to \code{link{sleuth_fit}}.
#'
#' @return a list containing two items:
#'   \itemize{
#'     \item{observed}{this is vector of the observed statistics; the likelihood ratio for "lrt", and
#'        wald statistic for "wt"}
#'     \item{resampled}{this is a matrix of number-of-features rows and numResamples columns, containing
#'        the resampled statistics}
#'   }
#' @importFrom dplyr bind_rows
#' @importFrom parallel mclapply
#' @importFrom sleuth sleuth_fit sleuth_lrt sleuth_wt
#' @export
calcGeneSleuthStats <- function(obj, numResamples = NULL, numCores = 1L, testType = "lrt", whichModel = "full", whichTest = "reduced:full",
                                whichBeta, ...) {
  stopifnot(is(obj, "sleuth"))
  if(testType == "lrt") {
    models <- strsplit(whichTest, ":", fixed = TRUE)[[1]]
    null_model <- models[1]
    if (whichModel != models[2]) {
      stop("'whichModel' and the alternative model in 'whichTest' do not match")
    }
    observedStats <- obj$tests$lrt[[whichTest]]$test_stat
  } else {
    if (whichTest != whichBeta) {
      stop("for 'testType' 'wt', 'whichTest' and 'whichBeta' must match")
    }
    observedStats <- obj$tests$wt[[whichModel]][[whichTest]]$wald_stat
  }

  if (missing(whichBeta)) {
    stop("'whichBeta' is missing. It must be specified to know which labels to permute")
  }

  if ((!is(numCores, "integer") && !is(numCores, "numeric")) || numCores <= 0) {
    stop("numCores is an invalid value")
  } else {
    numCores <- as.integer(numCores)
  }

  design <- obj$fits[[whichModel]]$design_matrix
  classLabels <- design[, whichBeta]
  n_ctrl <- sum(classLabels == 0)
  n_samples <- length(classLabels)
  max_perms <- choose(n_samples, n_ctrl) - 1
  if (is.null(numResamples)) {
    message("'numResamples' was not specified. Generating all unique combinations.")
    permMat <- returnResampleMat(classLabels, numResamples, allPerms = TRUE)
    numResamples <- max_perms
  } else if (numResamples == max_perms) {
    message(paste0("'numResamples', ", numResamples, ", is equal to the maximum number of unique ",
                   "combinations: ", max_perms, ". Generating all unique combinations."))
    permMat <- returnResampleMat(classLabels, numResamples, allPerms = TRUE)
  } else if (numResamples > max_perms) {
    warning(paste0("'numResamples', ", numResamples, ", is greater than the maximum number of unique ",
                   "combinations: ", max_perms, ". Generating all unique combinations."))
    numResamples <- max_perms
    permMat <- returnResampleMat(classLabels, numResamples, allPerms = TRUE)
  } else {
    permMat <- returnResampleMat(classLabels, numResamples, allPerms = FALSE)
  }

  permStats <- parallel::mclapply(1:numResamples, function(i) {
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
  names(permStats) <- paste0("resample_", 1:numResamples)
  permStats <- as.matrix(as.data.frame(permStats))

  list(observed = observedStats, resampled = permStats)
}


#' Return Matrix of Unique Resamplings
#'
#' This method takes a named vector of class labels, and
#' returns a matrix of the resampled class labels matched to
#' the original sample names. This only generates unique
#' combinations (i.e. ordering does not matter).
#'
#' @param classLabels a named vector with the class labels to be permuted.
#'   at this time, only two-condition experiments are allowed.
#' @param numResamples the number of desired unique resamplings;
#'   if this value is equal to or greater than the maximum possible
#'   unique resamplings, all unique resamplings will be provided.
#'   The default is \code{NULL}, indicating that all combinations should be
#'   generated.
#' @param allPerms if \code{TRUE} (the default), it will override \code{numResamples}
#'   and generate all possible unique combinations.
#'
#' @return a matrix with n rows by m columns, with n equal to the number of
#'   resamplings, and m equal to the length of the class labels.
#' @importFrom utils combn
#' @export
returnResampleMat <- function(classLabels, numResamples = NULL, allPerms = TRUE) {
  if (length(unique(classLabels)) > 2) {
    stop("This method is only implemented for experiments with two conditions")
  }

  if (is.null(names(classLabels))) {
    stop("This method only works if the 'classLabels' vector is named")
  }

  if (is.null(numResamples) && !allPerms) {
    stop("If 'allPerms' is FALSE, 'numResamples' must be specified.")
  }

  if (is(classLabels, "factor")) {
    orig_levels <- levels(classLabels)
    orig_class <- "factor"
    control <- orig_levels[1]
  } else {
    orig_levels <- NULL
    orig_class <- class(classLabels)
    classLabels <- as.factor(classLabels)
    control <- levels(classLabels)[1]
  }
  n_ctrl <- sum(classLabels == control)
  orig_names <- names(classLabels)
  n_samples <- length(classLabels)

  max_perms <- choose(n_samples, n_ctrl) - 1
  if (is.null(numResamples)) numResamples <- max_perms

  if (allPerms && numResamples != max_perms) {
    warning(paste0("allPerms is TRUE, but the specified 'numResamples', ", numResamples,
                   ", doesn't match the total possible number of combinations, ",
                   max_perms, ". Ignoring 'numResamples'..."))
  } else if (!allPerms && numResamples >= max_perms) {
    warning(paste0("The specified 'numResamples', ", numResamples, ", is equal to or ",
                   "larger than the total possible number of combinations, ", max_perms,
                   ". Generating all possible combinations and will return only ", max_perms,
                   " combinations."))
    allPerms <- TRUE
  } else if ((!allPerms && numResamples > 10000) || (allPerms && max_perms > 10000)) {
    samples <- ifelse(allPerms, max_perms, numResamples)
    message(paste0("The specified number of combinations, ", samples, ", is more than 10,000.",
                   " This may take a little bit of time to generate unique permutations.",
                   " Please contact the developers if this is significantly delaying your analysis."))
  }

  if (allPerms) {
    # This is to make sure that the two groups are separate; this allows the combinations to work
    permLabels <- classLabels[order(classLabels)]
    sample_names <- names(permLabels)
    # This creates a combination of the full number of samples
    # Then matches the combination to the sample names to get the new order
    # This then gives the new class labels
    permMat <- combn(n_samples, n_ctrl, function(x) {
      rest <- setdiff(1:n_samples, x)
      fullPerm <- c(x, rest)
      sampleOrder <- sample_names[fullPerm]
      newLabels <- as.character(permLabels)
      names(newLabels) <- sampleOrder
      newLabels[sample_names]
    })
    # Transpose and then remove the first row, which is the original order
    permMat <- t(permMat)[-1,]
    # Rearrange the order of the columns to match the original sample order
    colnames(permMat) <- sample_names
    permMat <- permMat[, orig_names]
    mode(permMat) <- orig_class
  } else {
    # This mini-function guarantees unique combinations by keeping track of
    # which combinations have been generated already.
    # This makes use of a 'cache', and the '<<-' function to update the cache
    # from inside the function.
    # For more information about how '<<-' works, see '?`<<-`'
    # This function was taken from this Stats StackExchange answer:
    # https://stats.stackexchange.com/a/24713
    newperm <- function() {
      count <- 0
      repeat {
        p <- sample(classLabels)
        hash.p <- paste(p, collapse = "")
        if (is.null(cache[[hash.p]])) {
         break
        }
        count <- count+1
        if (count > 1000) {
          warning(paste0("Could not find a unique combination after 1000 tries. ",
                  "Returning repeated combination."))
          break
        }
      }
      cache[[hash.p]] <<- TRUE
      p
    }
    cache <- list()
    permMat <- replicate(numResamples, newperm())
    permMat <- t(permMat)
    colnames(permMat) <- orig_names
    mode(permMat) <- orig_class
  }
  permMat
}
