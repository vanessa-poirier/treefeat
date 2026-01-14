# This script includes 1 function which can be used to extract ALL tree, trunk, and crown features 
# described in this package from a laz or las file of an individual tree

#' @title Extraction of tree architectural features
#' 
#' @description This function calculates tree archiectural features pertaining to  . In effect, it combines 
#' all functions available in the entireity of the package and adds a few. Particularly, some
#' ratios of width and height of different parts of the tree (e.g. ratio of crown height to
#' tree height).
#' 
#' @details The entire list of features rendered (totaling 23) is as follows:
#' tree height, crown height, crown width (at the upper limit, lower limit, and middle of the 
#' crown), trunk width, crown projected area, ratio of crown projected area to tree height, 
#' ratio of crown height to tree height, ratio of crown width to tree height, 
#' ratio of crown width to trunk width, ratio of upper crown width to lower crown width,
#' ratio of lower crown width to middle crown width, ratio of upper crown width to middle 
#' crown width, crown fullness, gap fraction, box dimension, curvature mean, curvature standard
#' deviation, curvature p90, rumple index.
#' 
#' @param file las/laz file or lidR object. A point cloud file of an individual tree
#' Designed for single-spectral TLS (terrestrial laser scanning) or 
#' high-resolution MLS (mobile laser scanning) data.
#' @param res Numeric. Resolution/intervals (i.e. height in meters) at which to slice the point cloud (for classifying crown points).
#' @param denoise Logical. Determines whether point cloud noise should be removed (determined by lidR::classify_noise).
#' @param max_crown_height Numeric. The maximum height at which the crown can start (presented as proportional height). 
#' Default is 0.4 (i.e. the crown must start =< 40% up the tree).
#' @param res_vox Numeric. Resolution of the voxels used to calculate crown projected area and gap fraction. 
#' Default is 0.1 meters.
#' @param npts_empty Numeric. Minimum number of points required for a voxel to be considered not empty (i.e. containing vegetation).
#' Default is 20 points. Value used in crown projected area and gap fraction calculations.
#' @param slice_radius Numeric. Half the height in meters of the slices used to calculate crown width and trunk widths. 
#' Default is 0.1 meters.
#' @param top_crown Numeric. Proportional crown height (values range from 0 to 1) at which the width of the top of the crown 
#' will be calculated. Default is 0.9.
#' @param middle_crown Numeric. Proportional crown height (values range from 0 to 1) at which the width of the middle of the crown
#' will be calculated. Default is 0.5.
#' @param bottom_crown Numeric. Proportional crown height (values range from 0 to 1) at which the width of the bottom of the crown
#' will be calculated. Default is 0.1.
#' @param top_trunk Numeric. Proportional crown height (values range from 0 to 1) at which the width of the top of the trunk
#' will be calculated. Default is 0.6.
#' @param middle_trunk Numeric. Proportional crown height (values range from 0 to 1) at which the width of the middle of the trunk
#' will be calculated. Default is 0.5.
#' @param bottom_trunk Numeric. Proportional crown height (values range from 0 to 1) at which the width of the bottom of the trunk
#' will be calculated. Default is 0.4.
#' @param k Numeric. Number of nearest neighbors to be considered when computing curvatures. Default is 20
#' @return A dataframe containing all tree architecural features.
#' @examples
#' tree <- lidR::readTLSLAS(file)
#' features <- extract_features(tree)
#' @export

extract_features <- function(file, res = 0.1, denoise = TRUE, max_crown_height = 0.4, res_vox = 0.1, npts_empty = 20,
                            slice_radius = 0.1, top_crown = 0.9, middle_crown = 0.5, bottom_crown = 0.1,
                            top_trunk = 0.6, middle_trunk = 0.5, bottom_trunk = 0.4, k = 20) {
    if (endsWith(file, ".las") | endsWith(file, ".laz")) {
        tree_pc <- lidR::readTLSLAS(file) # load point cloud
    } else if (class(file)[1] == "LAS") {
        tree_pc <- file
    } else {
        stop("Invalid file_type specified. Use las/laz file or lidR object.")
    }

    # Calculate tree height
    total_height <- tree_height(tree_pc)

    # Classify crown
    tree_pcs <- classify_crown(tree_pc, denoise = denoise)
    crown_pc <- tree_pcs$crown_pc
    trunk_pc <- tree_pcs$trunk_pc

    # Calculate crown height
    c_height <- crown_height(crown_pc)

    # Calculate crown projected area
    cpa <- crown_projected_area(crown_pc, res = res_vox, npts_empty = npts_empty)

    c_width <- crown_width(crown_pc, slice_radius = slice_radius, top_height = top_crown, middle_height = middle_crown, bottom_crown = bottom_height)

    t_width <- trunk_width(trunk_pc, slice_radius = slice_radius, top_height = top_trunk, middle_height = middle_trunk, bottom_crown = bottom_trunk)

    # Ratio of tree dimensions
    crown2tree_ratio_height <- c_height/total_height
    crown2tree_area_height <- cpa/total_height
    crown2trunk_ratio_diameter <- c_width[[2]]/t_width
    crowndiameter_treeheight_ratio <- c_width[[2]]/total_height
    crowndiameter_crownheight_ratio <- c_width[[2]]/c_height
    top_middle_ratio <- c_width[[1]]/c_width[[2]]
    bottom_middle_ratio <- c_width[[3]]/c_width[[2]]
    top_bottom_ratio <- c_width[[1]]/c_width[[3]]

    # Crown structural complexity
    crown_full <- crown_fullness(crown_pc, tree_pc)
    gap_fract <- gap_fraction(crown_pc, res = res_vox, npts_empty = npts_empty)
    fractal <- (rTwig::box_dimension(crown_pc))[[2]]$slope
    curvatures <- crown_curvature(crown_pc, k = k)
    rumple <- lidR::rumple_index(crown_pc$X, crown_pc$Y, crown_pc$Z)

    ## Put all metrics into one list ##
    return(list(tree_height = total_height, crown_height = c_height, crown_projected_area = cpa,
                trunk_diameter = t_width,
                crown2tree_ratio_height = crown2tree_ratio_height, 
                crown_area_tree_height_ratio = crown2tree_area_height,
                crown2trunk_ratio_diameter = crown2trunk_ratio_diameter, 
                crowndiameter_treeheight_ratio = crowndiameter_treeheight_ratio, 
                crowndiameter_crownheight_ratio = crowndiameter_crownheight_ratio, 
                top_diameter = c_width[[1]], middle_diameter = c_width[[2]], bottom_diameter = c_width[[3]], 
                top_middle_ratio = top_middle_ratio, bottom_middle_ratio = bottom_middle_ratio, 
                top_bottom_ratio = top_bottom_ratio, crown_fullness = crown_full, gap_fraction = gap_fract,
                box_dimension = fractal, curvature_mean = curvatures[[1]], curvature_sd = curvatures[[2]], 
                curvature_p90 = curvatures[[3]], rumple_index = rumple))
}
