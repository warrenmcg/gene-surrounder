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
  ## the minus 1 is to remove the original arrangement
  classLabels <- as.factor(classLabels)
  n_ctrl <- sum(classLabels == levels(classLabels)[1])
  n_samples <- length(classLabels)
  max_perms <- choose(n_samples, n_ctrl) - 1

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
