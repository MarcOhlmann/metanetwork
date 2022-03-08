#' Pipe 
#'
#' Like dplyr, metanetwork also uses the pipe function, \code{\%>\%} to turn
#' function composition into a series of imperative statements.
#'
#' @importFrom magrittr %>%
#' @name %>%
#' @rdname pipe
#' @export
#' @examples
#' 
#' library(metanetwork)
#' meta_angola %>% attach_layout() %>% ggmetanet()
NULL
