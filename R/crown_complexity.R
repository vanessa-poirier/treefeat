# This script includes 2 functions to calculate crown structural complexity measurement from a TLS or MLS point cloud of a single tree
# crown_pc and tree_pc are lidR objects
# Dependencies include VoxR and rTwig

#' @title Crown fullness
#' 
#' @description Proportion of points located in the crown compared to the entire tree, 
#' scaled by the height of the crown. ALternatively, a measure of crown density.
#' 
#' @param tree_pc lidR object. A point cloud file of an individual tree
#' Designed for single-spectral TLS (terrestrial laser scanning) or 
#' high-resolution MLS (mobile laser scanning) data.
#' @param crown_pc lidR object. A point cloud file of an individual tree crown.
#' Designed for single-spectral TLS (terrestrial laser scanning) or 
#' high-resolution MLS (mobile laser scanning) data.
#' Tree crown point cloud can be obtained from the classify_crown function.
#' @return Numeric. A proportion (values range from 0 to 1)
#' @export

crown_fullness <- function(crown_pc, tree_pc) {
    # crown_pc = lidR object. Point cloud of the tree's crown. Can be obtained by running classify_crown
    # tree_pc = lidR object. Point cloud of a single tree (TLS or high-resolution MLS is preferable). Must be the same tree as crown_pc

    crown_pts <- nrow(crown_pc)
    total_pts <- nrow(tree_pc)

    crown_height <- max(crown_pc[complete.cases(crown_pc$Z),]$Z) - min(crown_pc[complete.cases(crown_pc$Z),]$Z)

    crown_fullness <- (crown_pts/total_pts)/crown_height

    return(crown_fullness)
}


#' @title Gap fraction
#' 
#' @description The proportion of the 2D projected crown voxels that contain less than 20 points. 
#' Calculated by voxelating the crown (with a default resolution of 10 cm), projecting it into 2D space, 
#' and calculating the proportion of ‘gaps’, whereby gaps are characterized as all voxels 
#' containing less than 20 points.
#' 
#' @param crown_pc lidR object. A point cloud file of an individual tree crown.
#' Designed for single-spectral TLS (terrestrial laser scanning) or 
#' high-resolution MLS (mobile laser scanning) data.
#' Tree crown point cloud can be obtained from the classify_crown function.
#' @param res Numeric. Resolution of the voxels used to calculate crown projected area and gap fraction. 
#' Default is 0.1 meters.
#' @param npts_empty Numeric. Minimum number of points required for a voxel to be considered not empty (i.e. containing vegetation).
#' Default is 20 points. Value used in crown projected area and gap fraction calculations.
#' @return Numeric. A proportion (values range from 0 to 1)
#' @export

gap_fraction <- function(crown_pc, res = 0.1, npts_empty = 20) {
    # crown_pc = lidR object. Point cloud of the tree's crown. Can be obtained by running classify_crown
    # res = numeric. Resolution of the voxels. Default is 0.1 meters

    vox_pc <- VoxR::vox(crown_pc, res=res, full.grid=FALSE) # res is resolution of voxels

    vox_pc_xy <- VoxR::project_voxels(vox_pc, "xy") # create 2D grid from the voxels
    vox_pc_xy_gap <- subset(vox_pc_xy, npts<npts_empty) # create 2D grid from voxels that contain no more than 20 points

    #~Gap fraction (2D measurement) = proportion of 2D grid spaces (as seen from above) that contain voxels (res=0.1) with no more than 20 points
    gap_fraction <- nrow(vox_pc_xy_gap)/nrow(vox_pc_xy)

    return(gap_fraction)
}


#' @title Per-point curvatures
#' 
#' @description The coherence of each point as considered by its 20 nearest neighbors. 
#' Based on the eigenstructure method of computing 3D structural coherence. see Gersztenkorn and Marfurt 1999.
#' This function is meant to be used in the crown_curvature function.
#' 
#' @param points lidR object. A point cloud file of an object.
#' Designed for single-spectral TLS (terrestrial laser scanning) or 
#' high-resolution MLS (mobile laser scanning) data.
#' @param k Numeric. Number of nearest neighbors to be considered when computing curvatures. 
#' Default is 20.
#' @return Numeric. A vector of numbers.
#' @export

compute_curvature <- function(points, k = 20) {
    knn <- FNN::get.knn(points, k = k)
    curvatures <- as.numeric(nrow(points))
  
    for (i in 1:nrow(points)) {
        neighbor_idx <- knn$nn.index[i, ]
        neighbors <- points[neighbor_idx, ]
        
        cov_matrix <- stats::cov(neighbors)
        eigen_vals <- eigen(cov_matrix)$values
        
        # sort eigenvalues ascending: λ0 (smallest), λ1, λ2
        eigen_vals <- sort(eigen_vals)
        
        # Surface variation (curvature proxy)
        curvature <- eigen_vals[1] / sum(eigen_vals)
        curvatures[i] <- curvature
    }
    return(curvatures)
}


#' @title Crown curvature
#' 
#' @description The mean, standard deviation, and 90th percentile crown curvature 
#' as computed by compute_curvature when applied to crown_pc.
#' 
#' @param crown_pc lidR object. A point cloud file of an individual tree crown.
#' Designed for single-spectral TLS (terrestrial laser scanning) or 
#' high-resolution MLS (mobile laser scanning) data.
#' Tree crown point cloud can be obtained from the classify_crown function.
#' @param k Numeric. Number of nearest neighbors to be considered when computing curvatures. 
#' Default is 20.
#' @return List of 3 numeric values corresponding to the mean, standard deviation, 
#' and 90th percentile curvature values, respectively.
#' @export

crown_curvature <- function(crown_pc, k = 20) {
    # crown_pc = lidR object. Point cloud of the tree's crown. Can be obtained by running classify_crown
    # k = numeric. Number of nearest neighbors to be considered when computing curvatures

    curvatures <- compute_curvature(crown_pc, k = k) # compute individual curvatures (per point) with nearest neighbors k = 20
    mean_curvature <- mean(curvatures)
    sd_curvature <- sd(curvatures)
    p90_curvature <- quantile(curvatures, 0.9) # the threshold value below which 90% of the curvature values fall

    return(list(mean_curvature = mean_curvature, sd_curvature = sd_curvature, p90_curvature = p90_curvature))
}


## OTHER FUNCTIONS ALREADY AUTHORED IN OTHER PACKAGES BUT USED IN THIS PACKAGE ##

# Box dimension == rTwig function
# box_dimension <- function(crown_pc) {
#     # crown_pc = lidR object. Point cloud of the tree's crown. Can be obtained by running classify_crown

#     fractal <- (rTwig::box_dimension(crown_pc))[[2]]$slope
#     return(fractal)
# }


# Rumple_index == lidR function
# rumple_index <- function(crown_pc) {
#     # crown_pc = lidR object. Point cloud of the tree's crown. Can be obtained by running classify_crown

#     rumple <- lidR::rumple_index(crown_pc$X, crown_pc$Y, crown_pc$Z)

#     return(rumple)
# }
