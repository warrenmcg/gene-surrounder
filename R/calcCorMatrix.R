#' Calculate large gene expression correlation matrix in parallel
#'
#' This method uses the ff package to calculate a large
#' correlation matrix in parallel.
#'
#' @param exprMatrix a N x M matrix or data.frame, with N columns corresponding to
#'   each samples in the experiment observed, and M rows corresponding to the
#'   the genes observed.
#' @param exprName A string with the name of the expression matrix.
#' @param size the number of elements per block to split the data. The
#'   default is 2000, which results in decent performance for typical
#'   usage (order of 10-20K features). If the number of features is
#'   not evenly divisible by the size parameter, a remainder block is used.
#' @param numCores the number of cores to use for parallel calculations
#' @param corMethod string with the method of calculating correlation. Supported
#'   methods are 'pearson' (default), 'spearman', or 'kendall'. See
#'   ?stats::cor for more details.
#' @param useMethod string with what observations should be used. Supported
#'   options are "everything" (default), "all.obs", "complete.obs",
#'   "na.or.complete", or "pairwise.complete.obs". See ?stats::cor for more
#'   details on what each option entails.
#'
#' @return an 'ff' object pointing to a file on disk containing the
#'   gene correlation matrix
#' @importFrom ff ff
#' @importFrom parallel mclapply
#' @importFrom stats cor
#' @export
calcCorMatrix <- function(exprMatrix, exprName, size = 2000, numCores = 1L,
                        corMethod = "pearson", useMethod = "everything") {
  ## this code is adapted from two sources:
  ## 1) bigcor.R: http://www.dr-spiess.de/scripts/bigcor.R
  ## 2) bigcorPar.R: https://gist.github.com/bobthecat/5024079
  corMethod <- match.arg(corMethod, c('pearson', 'spearman', 'kendall'))
  useMethod <- match.arg(useMethod, c('everything', 'all.obs', 'complete.obs',
                                      'na.or.complete',
                                      'pairwise.complete.obs'))
  if (ncol(data) > nrow(data)) {
    stop(paste0("There are more columns than row. Are you ",
                "sure that each sample is a column, and each ",
                "transcript is a row?"))
  } else {
    data <- as.data.frame(t(data))
  }

  g_names <- colnames(data)
  if (is.null(t_names)) {
    warning(paste0('There are no transcript IDs as column names. The ',
                   'correlation matrix will have no IDs for row names or ',
                   'column names.'))
  }

  NCOL <- ncol(data)
  rest <- NCOL %% size
  large <- NCOL - rest
  nblocks <- NCOL %/% size

  group <- rep(1:nblocks, each = size)
  if (rest > 0) group <- c(group, rep(nblocks + 1, rest))

  SPLIT <- split(1:NCOL, group)

  COMBS <- expand.grid(1:length(SPLIT), 1:length(SPLIT))
  COMBS <- t(apply(COMBS, 1, sort))
  COMBS <- unique(COMBS)

  corMAT <- ff::ff(vmode = 'single', dim = c(NCOL, NCOL),
                   dimnames = list(g_names, g_names))
  attr(corMAT, "expr") <- exprName
  attr(corMAT, "method") <- corMethod

  result <- parallel::mclapply(1:nrow(COMBS), function(i) {
    comb <- COMBS[i, ]
    g1 <- SPLIT[[comb[1]]]
    g2 <- SPLIT[[comb[2]]]
    COR <- stats::cor(data[, g1], data[, g2], method = corMethod, use = useMethod)
    corMAT[g1, g2] <- COR
    corMAT[g2, g1] <- t(COR)
  }, mc.cores = numCores)
  suppressMessages(gc())

  return(corMAT)
}
