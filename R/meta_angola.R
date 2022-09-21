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


#' Angola fishery metanetwork
#' metanetwork built from:
#' Angelini & Velho 2011, Data from: 
#' Angelini, R., Velho, VF. (2011)
#'  Ecosystem structure and trophic analysis of Angolan fishery landings. Scientia Marina 75(2)
#' @docType data
#'
#' @format
#' \describe{
#' A object of class 'metanetwork'
#' \item{metaweb}{The metaweb from Angelini & Velho 2011, containing 28 groups and 127 interactions, a \code{igraph} object} 
#' \item{abTable}{Abundance table built from biomass at two dates: 1986 and 2003, a \code{matrix}}
#' \item{trophicTable}{Taxonomic table, a three column \code{data.frame} 
#' with three different taxonomic levels (species (or group), phylum and kingdom)}
#' }
#' @usage data(meta_angola)
#' @source \url{https://scientiamarina.revistas.csic.es/index.php/scientiamarina/article/view/1254}
"meta_angola"

