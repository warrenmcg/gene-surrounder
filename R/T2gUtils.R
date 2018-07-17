#' Aggregate transcript-level p-values to gene-level
#'
#' This takes a vector or a matrix of transcript-level
#' p-values and aggregates them to the gene-level.
#'
#' @param pvals vector or matrix of transcript-level p-values.
#'   if a vector, names should be transcript IDs; if a matrix,
#'   rownames should be transcript IDs.
#' @param mean_obs named vector of the mean of the observed normalized
#'   expression for each transcript
#' @param t2g a data frame with two columns; the first contains
#'   the transcript IDs, and should match the names of pvals; the second
#'   contains the gene IDs. All transcripts should have one
#'   and only one corresponding gene ID.
#' @param weight_func which function should be applied to mean_obs
#'   to weight the p-values for aggregation using the lancaster method?
#'   the default is the identity function. Note that whatever
#'   function is applied should yield only positive values.
#' @param numCores the number of cores to use. If a p-value matrix is supplied
#'   and numCores is set to 1, this function will use data.table to calculate
#'   the matrix gene-level matrix. If it's more than one, it will use
#'   parallel's mclapply, which is faster based on benchmarks.
#'
#' @return a vector or matrix of gene-level p-values. The names of
#'   the vector, or the rownames of the matrix, will be the provided
#'   gene IDs.
#'
#' @importFrom aggregation lancaster
#' @importFrom data.table as.data.table
#' @importFrom parallel mclapply detectCores
#' @export
transToGenePvals <- function(pvals, mean_obs, t2g, weight_func = identity, numCores = 1) {
  if (is(pvals, "numeric")) {
    target_ids <- names(pvals)
  } else if (is(pvals, "matrix") && is(pvals[,1], "numeric")) {
    target_ids <- rownames(pvals)
  } else {
    stop("pvals is not a numeric vector or a matrix")
  }

  stopifnot(is(numCores, "integer") || is(numCores, "numeric"))
  stopifnot(numCores >= 1 && numCores <= parallel::detectCores())
  numCores <- as.integer(numCores)

  if (!is(t2g, "data.frame")) {
    stop("'t2g' must be a data.frame")
  } else if (ncol(t2g) != 2) {
    stop(paste0("'t2g' must have two columns: the first for the transcript IDs, ",
                "and the second for the corresponding gene IDs"))
    colnames(t2g) <- c("target_id", "gene_id")
  }

  if (is.null(target_ids)) {
    stop("'pvals' must have names if it is a vector, or rownames if it is matrix")
  } else if (!all(target_ids %in% t2g$target_id)) {
    stop(paste0("At least one of the transcript IDs in the provided ",
                "'pvals' cannot be found in the 't2g' data frame.\n",
                "Make sure that the first column of 't2g' contains ",
                "transcript IDs and that all of the pvals transcript IDs ",
                "have an entry in 't2g'."))
  }

  if (is.null(names(mean_obs))) {
    stop("'mean_obs' must be a named vector using transcript IDs for names")
  } else if (!all(target_ids %in% names(mean_obs))) {
    stop(paste0("At least one of the transcript IDs in the provided ",
                "'pvals' cannot be found among the names of 'mean_obs'"))
  } else {
    mean_obs <- mean_obs[target_ids]
  }

  ## If more than one core is used, and the supplied p-values are a matrix,
  ## the faster method using parallel's mclapply is used
  if (numCores > 1 && !is.null(dim(pvals))) {
    new_pval_list <- parallel::mclapply(1:ncol(pvals), function(i) {
      df <- data.frame(target_id = target_ids, mean_obs = mean_obs, pval = pvals[,i])
      dt <- merge(t2g, data.table::as.data.table(df), by = "target_id")
      pval_dt <- dt[, list(new_pval = as.numeric(aggregation::lancaster(pval, weight_func(mean_obs)))), by = "gene_id"]
      pvals <- pval_dt[["new_pval"]]
      names(pvals) <- pval_dt[["gene_id"]]
      pvals
    }, mc.cores = numCores)
    new_pvals <- as.matrix(as.data.frame(new_pvals))
  } else {
    ## If only one core is used, or if the supplied p-values are a vector,
    ## the faster method using data.table is used
    pvals_df <- data.frame(target_id = target_ids, mean_obs = mean_obs, pvals)
    
    t2g <- data.table::as.data.table(t2g)
    full_dt <- merge(t2g, data.table::as.data.table(pvals_df), by = "target_id")
    pval_cols <- colnames(full_dt)[-1:-3]
    ## the code below takes the pval_cols, and uses the
    ## p-value aggregation method from the aggregation package
    ## on each column, producing a new aggregated data table
    ## with the gene_id as the first column, and the aggregated
    ## p-values as the other columns.
    ## see ?aggregation::lancaster for details on that method
    ## see stackoverflow answer here for more details on how lapply
    ## works with a data table:
    ## https://stackoverflow.com/a/16513949
    new_dt <- full_dt[,
      lapply(.SD, function(col) {
        as.numeric(aggregation::lancaster(col, weight_func(mean_obs)))
      }),
      by = "gene_id",
      .SDcols = pval_cols
    ]
    new_pvals <- new_dt[,..pval_cols]
    if (is.null(dim(new_pvals))) {
      names(new_pvals) <- new_dt[["gene_id"]]
    } else {
      rownames(test) <- new_dt[["gene_id"]]
    }
  }
  new_pvals
}

