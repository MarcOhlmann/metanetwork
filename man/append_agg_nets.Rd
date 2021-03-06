% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/append_agg_nets.R
\name{append_agg_nets}
\alias{append_agg_nets}
\alias{append_agg_nets.metanetwork}
\title{append aggregated networks}
\usage{
append_agg_nets(metanetwork)

\method{append_agg_nets}{metanetwork}(metanetwork)
}
\arguments{
\item{metanetwork}{object of class 'metanetwork'}
}
\value{
object of class 'metanetwork'

\code{NULL}
}
\description{
Method to append aggregated metawebs and aggregated local networks using
the hierarchy describes in \code{trophicTable}
}
\examples{
library(metanetwork)
data(meta_angola)
append_agg_nets(meta_angola)

}
