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


# extended Moran index
get_ext_Moran_index <- function(y_coord,g,beta = NULL){
    A = as.matrix(get.adjacency(g,attr = "weight"))
    X = y_coord
    Xmean = mean(X)
    n = vcount(g)
    W = sum(A)
    I = (n/W) * (t(A%*%(X-Xmean))%*%X) / (t(X)%*%X)
    #penalty with number of distance inferior to label size
    x_y_coords = cbind(V(g)$TL,y_coord)
    sum(dist(x_y_coords) < ggnet.default$label.size)/length(dist(x_y_coords))
    I = I - sum(dist(x_y_coords) < ggnet.default$label.size)/length(dist(x_y_coords))
  return(I)
}

get_coords_list <- function(g,beta,nrep){
  K = compute_diffusion_kernel(g,beta)
  coords_list = lapply(1:nrep,function(k)
    return(
      get_nodes_position_TL_tsne(
        g,V(g)$TL,beta = beta,
        TL_tsne.config = TL_tsne.default)
    )
  )
}




    
