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


#' European vertebrates metanetwork
#' metanetwork built using data from:
#' O'Connor, L. M., Pollock, L. J., Braga, J., Ficetola, G. F., Maiorano, L., Martinez-Almoyna, C., ...  & Thuiller, W. (2020). 
#' Unveiling the food webs of tetrapods across Europe through the prism of the Eltonian niche. Journal of Biogeography, 47(1), 181-192.
#' and
#' Maiorano, L., Montemaggiori, A., Ficetola, G. F., O'connor, L., & Thuiller, W. (2020). 
#' TETRA-EU 1.0: a species-level trophic metaweb of European tetrapods. Global Ecology and Biogeography, 29(9), 1452-1457.
#' @docType data
#'
#' @format
#' \describe{
#' A object of class 'metanetwork'
#' \item{metaweb}{The metaweb from Maiorano et al. 2020, O'Connor et al 2020, containing 1101 species and 49013 interactions, a \code{igraph} object} 
#' \item{trophicTable}{Trophic table, a two columns \code{data.frame} 
#' with a column containing species name and a column containing Stochastic Block Model groups inferred in O'Connor et al 2020}
#' }
#' @usage data(meta_vrtb)
#' @source \url{https://onlinelibrary.wiley.com/doi/abs/10.1111/geb.13138}, \url{https://onlinelibrary.wiley.com/doi/abs/10.1111/jbi.13773}
"meta_vrtb"

