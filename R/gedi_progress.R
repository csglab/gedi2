# ==============================================================================
# GEDI Progress and Logging Utilities
# Unified interface for progress bars and status messages
# ==============================================================================

#' Create Progress Bar (Method 4: Base R txtProgressBar)
#'
#' @description
#' Creates a simple text-based progress bar for iterative operations.
#' Use this for loops over samples or other countable iterations.
#'
#' @param total Integer, total number of iterations
#' @param label Character, optional label to display (default: "Progress")
#' @param verbose Integer, verbosity level (0 = silent, 1+ = show progress)
#'
#' @return Progress bar object or NULL if verbose = 0
#'
#' @keywords internal
#' @noRd
create_progress_bar <- function(total, label = "Progress", verbose = 1) {
  if (verbose < 1 || total <= 1) {
    return(NULL)
  }

  if (nchar(label) > 0) {
    message(label)
  }

  pb <- txtProgressBar(min = 0, max = total, style = 3, width = 50)
  return(pb)
}


#' Update Progress Bar
#'
#' @description
#' Updates the progress bar to the current iteration.
#'
#' @param pb Progress bar object (from create_progress_bar)
#' @param current Integer, current iteration number
#'
#' @return NULL (invisible)
#'
#' @keywords internal
#' @noRd
update_progress_bar <- function(pb, current) {
  if (!is.null(pb)) {
    setTxtProgressBar(pb, current)
  }
  invisible(NULL)
}


#' Close Progress Bar
#'
#' @description
#' Closes the progress bar and prints a newline.
#'
#' @param pb Progress bar object (from create_progress_bar)
#'
#' @return NULL (invisible)
#'
#' @keywords internal
#' @noRd
close_progress_bar <- function(pb) {
  if (!is.null(pb)) {
    close(pb)
  }
  invisible(NULL)
}


#' Log Stage Message (Method 7: Progress with custom messages)
#'
#' @description
#' Prints a stage message for multi-stage operations.
#' Format: "  [Stage X/Y] Message"
#'
#' @param stage Integer, current stage number
#' @param total Integer, total number of stages
#' @param message Character, description of this stage
#' @param verbose Integer, verbosity level (0 = silent, 1+ = show)
#'
#' @return NULL (invisible)
#'
#' @keywords internal
#' @noRd
log_stage <- function(stage, total, message, verbose = 1) {
  if (verbose >= 1) {
    message(sprintf("  [%d/%d] %s", stage, total, message))
  }
  invisible(NULL)
}


#' Log Start Message
#'
#' @description
#' Prints a start message for an operation.
#' Format: "Computing [operation]..."
#'
#' @param operation Character, name of the operation
#' @param verbose Integer, verbosity level (0 = silent, 1+ = show)
#'
#' @return NULL (invisible)
#'
#' @keywords internal
#' @noRd
log_start <- function(operation, verbose = 1) {
  if (verbose >= 1) {
    message("Computing ", operation, "...")
  }
  invisible(NULL)
}


#' Log Complete Message
#'
#' @description
#' Prints a completion message for an operation.
#' Format: "[OK] [operation] computed: [details]"
#'
#' @param operation Character, name of the operation
#' @param details Character, optional details about dimensions/size
#' @param verbose Integer, verbosity level (0 = silent, 1+ = show)
#'
#' @return NULL (invisible)
#'
#' @keywords internal
#' @noRd
log_complete <- function(operation, details = NULL, verbose = 1) {
  if (verbose >= 1) {
    if (is.null(details)) {
      message("[DONE] ", operation, " complete")
    } else {
      message("[DONE] ", operation, " computed: ", details)
    }
  }
  invisible(NULL)
}


#' Log Info Message
#'
#' @description
#' Prints an informational message.
#' Format: "  [info]"
#'
#' @param message Character, the message to display
#' @param verbose Integer, verbosity level (0 = silent, 1+ = show)
#' @param indent Logical, whether to indent the message (default: TRUE)
#'
#' @return NULL (invisible)
#'
#' @keywords internal
#' @noRd
log_info <- function(message, verbose = 1, indent = TRUE) {
  if (verbose >= 1) {
    if (indent) {
      message("  ", message)
    } else {
      message(message)
    }
  }
  invisible(NULL)
}


#' Log Cached Message
#'
#' @description
#' Prints a message indicating cached data is being used.
#' Format: "[CACHED] Using cached [item]"
#'
#' @param item Character, name of the cached item
#' @param verbose Integer, verbosity level (0 = silent, 1+ = show)
#'
#' @return NULL (invisible)
#'
#' @keywords internal
#' @noRd
log_cached <- function(item, verbose = 1) {
  if (verbose >= 1) {
    message("[CACHED] Using cached ", item)
  }
  invisible(NULL)
}


#' Format Dimensions Message
#'
#' @description
#' Creates a formatted string for matrix dimensions.
#' Format: "X rows x Y cols" or "X genes x Y cells"
#'
#' @param nrow Integer, number of rows
#' @param ncol Integer, number of columns
#' @param row_label Character, label for rows (default: "rows")
#' @param col_label Character, label for columns (default: "cols")
#'
#' @return Character string with formatted dimensions
#'
#' @keywords internal
#' @noRd
format_dims <- function(nrow, ncol, row_label = "rows", col_label = "cols") {
  sprintf("%d %s x %d %s", nrow, row_label, ncol, col_label)
}
