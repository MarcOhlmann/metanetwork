#' Pipe 
#'
#' Like dplyr, metanetwork also uses the pipe function, \code{\%>\%} to turn
#' function composition into a series of imperative statements.
#'
#' @importFrom magrittr %>%
#' @name %>%
#' @rdname pipe
#' @return an object of the class of the output of the last called method/function
#' @export
#' @examples
#' 
#' library(metanetwork)
#' data("meta_angola")
#' meta_angola %>% attach_layout() %>% ggmetanet()
NULL
