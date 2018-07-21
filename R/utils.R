#' Check number of cores
#'
#' This internal function checks if the number of cores
#' specified is reasonable.
#'
#' @param numCores the number of cores
#'
#' @return the corrected number of cores
check_cores <- function(numCores) {
  if (is.null(numCores) || is.na(suppressWarnings(as.integer(numCores))) ||
      numCores < 1 || numCores > parallel::detectCores()) {
    warning(paste0("'numCores' is an invalid value. It must be a number between ",
                "1 and the number of cores on your machine. 'numCores is being ",
                "set to 1 so your analysis can run."))
    numCores <- 1
  } else if (numCores > 1 && Sys.getenv("RSTUDIO") == "1") {
    warning(paste0("It appears that you are running Gene-Surrounder from within ",
                   "Rstudio.\nBecause of concerns with forking processes from a ",
                   "GUI, 'num_cores' is being set to 1.\nIf you wish to take ",
                   "advantage of multiple cores, please consider running ",
                   "Gene-Surrounder from the command line."))
    numCores <- 1
  } else {
    numCores <- as.integer(numCores)
  }

  numCores
}


#' Get sleuth Stats
#'
#' This internal function retrieves
#' the specified test statistics from
#' the sleuth object.
#'
#' @param obj the sleuth object. It must have 'whichModel' already fit
#'   using \code{sleuth_fit}, as well as the null model if the testType is
#'   'lrt'.
#' @param testType either "lrt" or "wt" to indicate which test, Likelihood Ratio Test or
#'   Wald Test, to use to calculate statistics
#' @param whichModel if "lrt" is the test type, this the alternative model; if "wt" is
#'   the test type, this is the model to use.
#' @param whichTest if "lrt", this is the test in the format of "[null model]:[alternative model]";
#'   if "wt", this is the same as whichBeta, the actual column in the design matrix to test
#'
#' @return a vector with the statistics from the specified test.
#' @importFrom sleuth sleuth_lrt sleuth_wt
get_sleuth_stats <- function(obj, testType, whichModel, whichTest) {
  # If the main model has not been fit yet, stop
  if (is.null(obj$fits[[whichModel]])) {
    stop(paste0("This sleuth object has not been fitted with the specified model, '", whichModel,
                "'. Please run 'sleuth_fit' with this model."))
  }

  if(testType == "lrt") {
    models <- strsplit(whichTest, ":", fixed = TRUE)[[1]]
    null_model <- models[1]
    # If the null model has not been fit yet, stop
    if (is.null(obj$fits[[null_model]])) {
      stop(paste0("This sleuth object has not been fitted with the null model, '", null_model,
                  "'. Please run 'sleuth_fit' with this model."))
    }
    if (whichModel != models[2]) {
      stop("'whichModel' and the alternative model in 'whichTest' do not match")
    }
    if (is.null(obj$tests$lrt[[whichTest]])) {
      message(paste0("Specified likelihood ratio test '",
                     whichTest, "' is missing. Generating it now..."))
      obj <- sleuth::sleuth_lrt(obj, null_model, whichModel)
    }
    observedStats <- obj$tests$lrt[[whichTest]]$test_stat
  } else {
    # if the testing has not been done yet, do it for them
    if (is.null(obj$tests$wt[[whichModel]][[whichTest]])) {
      message(paste0("Specified Wald test '",
                     whichTest, "' is missing. Generating it now..."))
      obj <- sleuth::sleuth_wt(obj, whichTest, whichModel)
    }
    observedStats <- obj$tests$wt[[whichModel]][[whichTest]]$wald_stat
  }
  observedStats
}

#' Check Resample Size
#'
#' This internal function checks to make sure
#' the specified number of resamples is correct.
#' If it was not specified, then it is set as the
#' maximum number of permutations.
#'
#' @param classLabels the vector of labels for the 
#'   experiment to determine maximum possible permutations
#' @param numResamples the number of resamples
#' @param allPerms boolean for whether all permutations should be
#'   generated
#'
#' @return a list with the corrected numResamples and allPerms
#' @importFrom methods is
check_resample_size <- function(classLabels, numResamples, allPerms) {
  ## Sanity checks for numResamples
  if (!is.null(numResamples)) {
    stopifnot(is(numResamples, 'numeric') || is(numResamples, 'integer'))
    stopifnot(length(numResamples) == 1)
    stopifnot(numResamples > 0)
  }
  ## Sanity checks for allPerms
  stopifnot(is.logical(allPerms))
  stopifnot(length(allPerms) == 1)

  ## coerce numResamples to an integer
  numResamples <- as.integer(numResamples)

  ## determine the maximum number of permutations
  ## the minus 2 is to remove the original arrangement
  ## as well as its direct inverse
  classLabels <- as.factor(classLabels)
  n_ctrl <- sum(classLabels == levels(classLabels)[1])
  n_samples <- length(classLabels)
  max_perms <- choose(n_samples, n_ctrl) - 2

  if (is.null(numResamples)) {
    message("'numResamples' was not specified. Generating all unique combinations.")
    numResamples <- max_perms
    allPerms <- TRUE
  } else if (numResamples == max_perms) {
    message(paste0("'numResamples', ", numResamples, ", is equal to the maximum number of unique ",
                   "combinations: ", max_perms, ". Generating all unique combinations."))
    allPerms <- TRUE
  } else if (numResamples > max_perms) {
    warning(paste0("'numResamples', ", numResamples, ", is greater than the maximum number of unique ",
                   "combinations: ", max_perms, ". Generating all unique combinations."))
    numResamples <- max_perms
    allPerms <- TRUE
  } else if (allPerms && numResamples < max_perms) {
    warning(paste0("allPerms is TRUE, but the specified 'numResamples', ", numResamples,
                   ", doesn't match the total possible number of combinations, ",
                   max_perms, ". Ignoring 'numResamples'..."))
  }

  if (numResamples > 10000) {
    message(paste0("The specified number of combinations, ", numResamples, ", is more than 10,000.",
                   " This may take a little bit of time to generate unique permutations.",
                   " Please contact the developers if this is significantly delaying your analysis."))
  }

  list(numResamples, allPerms)
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
  } else {
    orig_levels <- NULL
    orig_class <- class(classLabels)
    classLabels <- as.factor(classLabels)
  }

  check_res <- check_resample_size(classLabels, numResamples, allPerms)
  numResamples <- check_res$numResamples
  allPerms <- check_res$allPerms

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
    # remove also the last row, which is the direct opposite of the original order
    permMat <- t(permMat)[-1, ]
    permMat <- permMat[-nrow(permMat), ]
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
