# This file is part of metanetwork

# metanetwork is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# metanetwork is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with metanetwork.  If not, see <http://www.gnu.org/licenses/>


#' Norway soil metanetwork
#' metanetwork built from:
#' Calderon-Sanou et al. 2021, Data from: 
#' Calderon-Sanou, I., Munkemuller, T., Zinger, L., Schimann, H., Yoccoz, N. G., Gielly, L., ... & Thuiller, W. (2021).
#'  Cascading effects of moth outbreaks on subarctic soil food webs. Scientific reports, 11(1), 1-12.
#' @docType data
#'
#' @format
#' \describe{
#' A object of class 'metanetwork'
#' \item{metaweb}{The metaweb from Calderon-Sanou et al. 2021, containing 40 groups and 204 interactions, a \code{igraph} object} 
#' \item{abTable}{Abundance table built from eDNA data in disturbed (moth outbreaks) and non-disturbed sites, a \code{matrix}}
#' \item{trophicTable}{Trophic table, a three column \code{data.frame} 
#' with three different taxonomic levels (trophic_group, trophic_class and taxa)}
#' }
#'
#' @usage data(meta_norway)
#' @source \url{https://www.nature.com/articles/s41598-021-94227-z}
"meta_norway"

