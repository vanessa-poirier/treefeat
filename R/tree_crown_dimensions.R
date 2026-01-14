# This script includes 6 functions to calculate absolute tree and crown dimensions from a TLS or MLS point cloud of a single tree
# Output units will match the units of the X, Y, Z coordinates of the point cloud. 
# N.B. If X, Y, Z are not in meters, slice_radius must be modified to reflect the same units as X, Y, Z.
# dependencies: VoxR

#' @title Tree height
#' 
#' @description Calculates tree height as the difference between the minimum and maximum Z values of the point cloud. 
#' 
#' @param tree_pc lidR object. A point cloud file of an individual tree
#' Designed for single-spectral TLS (terrestrial laser scanning) or 
#' high-resolution MLS (mobile laser scanning) data.
#' @return Numeric. Tree height (units match units of the point cloud X, Y, Z coordinates)
#' @export

tree_height <- function(tree_pc) {
    # tree_pc = lidR object. Point cloud of a single tree (TLS or high-resolution MLS is preferable)

    tree_height <- max(tree_pc[complete.cases(tree_pc$Z),]$Z) - min(tree_pc[complete.cases(tree_pc$Z),]$Z)
    return(tree_height)
}


#' @title Crown height
#' 
#' @description Calculates crown height as the difference between the minimum and maximum Z values of the 
#' tree crown's point cloud. 
#' 
#' @param crown_pc lidR object. A point cloud file of an individual tree crown.
#' Designed for single-spectral TLS (terrestrial laser scanning) or 
#' high-resolution MLS (mobile laser scanning) data.
#' Tree crown point cloud can be obtained from the classify_crown function.
#' @return Numeric. Crown height (units match units of the point cloud X, Y, Z coordinates)
#' @export

crown_height <- function(crown_pc) {
    # crown_pc = lidR object. Point cloud of the tree's crown. Can be obtained by running classify_crown

    crown_height <- max(crown_pc[complete.cases(crown_pc$Z),]$Z) - min(crown_pc[complete.cases(crown_pc$Z),]$Z)
    return(crown_height)
}


#' @title Crown projected areaList containing 3 numeric values corresponding to the crown width at the top, middle, and bottom of the crown, respectively.
#' 
#' @description Calculates crown projected area as the area the crown would occupy if lain flat on the ground. 
#' In other words, the 2D area of the crown as seen from above. Calculated by voxelating the crown, 
#' projecting it into 2D space, and multiplying the number of pixels by the resolution of the voxels.
#' 
#' @param crown_pc lidR object. A point cloud file of an individual tree crown.
#' Designed for single-spectral TLS (terrestrial laser scanning) or 
#' high-resolution MLS (mobile laser scanning) data.
#' Tree crown point cloud can be obtained from the classify_crown function.
#' @param res Numeric. Resolution of the voxels used to calculate crown projected area and gap fraction. 
#' Default is 0.1 meters.
#' @param npts_empty Numeric. Minimum number of points required for a voxel to be considered not empty (i.e. containing vegetation).
#' Default is 20 points. Value used in crown projected area and gap fraction calculations.
#' @return Numeric. A proportion (values range from 0 to 1).
#' @export

crown_projected_area <- function(crown_pc, res = 0.1, npts_empty = 20) {
    # crown_pc = lidR object. Point cloud of the tree's crown. Can be obtained by running classify_crown
    # res = numeric. Resolution of the voxels. Default is 0.1 meters
    # npts_empty = numeric. Minimum number of points required for a voxel to be considered not empty (i.e. containing vegetation)

    vox_pc <- VoxR::vox(crown_pc, res=res, full.grid=FALSE) # res is resolution of the voxels

    vox_pc_xy <- VoxR::project_voxels(vox_pc, "xy") # create 2D grid from the voxels
    vox_pc_xy_gap <- subset(vox_pc_xy, npts<npts_empty) # create 2D grid from voxels that contain no more than 20 points

    # Crown projected area
    crown_proj_area <- nrow(vox_pc_xy)*(res^2) # multiply the number of pixels making up the crown by the area of the voxels

    return(crown_proj_area)
}


#' @title Crown width/diameter
#' 
#' @description Crown width is calculated by calculating crown radius at 5° intervals, 
#' computing the mean, and multiplying this mean by 2.
#' 
#' @details Details of the calculation are as follows:
#' (1) Isolated the middle of the crown by slicing the crown at its midpoint (at 50% of its maximum height), 
#' keeping 10 cm above and below for a total slice thickness of 20 cm. (2) Measured the radial distance 
#' from the tree’s central axis to the furthest point in the sliced crown. (3) Repeated this measure in 
#' 5° intervals along the center axis and calculated the average. (4) Multiplied this average by two to 
#' estimate crown diameter. 
#' 
#' @param crown_pc lidR object. A point cloud file of an individual tree crown.
#' Designed for single-spectral TLS (terrestrial laser scanning) or 
#' high-resolution MLS (mobile laser scanning) data.
#' Tree crown point cloud can be obtained from the classify_crown function.
#' @param slice_radius Numeric. Half the height in meters of the slices used to calculate crown widths. 
#' Default is 0.1 meters.
#' @param top_height Numeric. Proportional crown height (values range from 0 to 1) at which the width of the top of the crown 
#' will be calculated. Default is 0.9.
#' @param middle_height Numeric. Proportional crown height (values range from 0 to 1) at which the width of the middle of the crown
#' will be calculated. Default is 0.5.
#' @param bottom_height Numeric. Proportional crown height (values range from 0 to 1) at which the width of the bottom of the crown
#' will be calculated. Default is 0.1.
#' @return List containing 3 numeric values corresponding to the crown width at the top, middle, and bottom of the crown, respectively.
#' @export

crown_width <- function(crown_pc, slice_radius = 0.1, top_height = 0.9, middle_height = 0.5, bottom_height = 0.1) {
    # crown_pc = lidR object. Point cloud of the tree's crown. Can be obtained by running classify_crown. X, Y, Z coordinates should be in meters.
    # slice_radius = Height in meters to be included below and above the slicing point. Default is 0.1 meters.
    # top_height = numeric. Proportional crown height (values range from 0 to 1) at which the width of the top of the crown will be calculated
    # middle_height = numeric. Proportional crown height (values range from 0 to 1) at which the width of the middle of the crown will be calculated
    # bottom_height = numeric. Proportional crown height (values range from 0 to 1) at which the width of the bottom of the crown will be calculated

    crown_pc <- data.frame(X = crown_pc$X, Y = crown_pc$Y, Z = crown_pc$Z) # make into xyz dataframe

    # find heights at which corresponds to the low, high, and middle part of the crown (as determined by top_height, middle_height, and bottom_height, respectively)
    crown_high <- ((max(crown_pc$Z) - min(crown_pc$Z))*top_height) + min(crown_pc$Z) # height top_height% up the crown
    crown_middle <- ((max(crown_pc$Z) - min(crown_pc$Z))*middle_height) + min(crown_pc$Z) # height (Z coordinate) middle_height% up the tree's crown
    crown_low <- ((max(crown_pc$Z) - min(crown_pc$Z))*bottom_height) + min(crown_pc$Z) # height bottom_height% up the crown

    crown_heights <- list(crown_high, crown_middle, crown_low)

    ## calculate widths at these crown heights
        # Create slices at each crown height
        # get points with Z +/- slice thickness for crown middle, high, and low
            crown_high_slice <- subset(crown_pc, Z<(crown_high+slice_radius) & Z>(crown_high-slice_radius)) # note that the xyz coords are in meters
            crown_middle_slice <- subset(crown_pc, Z<(crown_middle+slice_radius) & Z>(crown_middle-slice_radius))
            crown_low_slice <- subset(crown_pc, Z<(crown_low+slice_radius) & Z>(crown_low-slice_radius))

            crown_slices <- list(crown_high_slice, crown_middle_slice, crown_low_slice)

            # For each crown slice:
            crown_diameters <- list()

                for (c in 1:length(crown_slices)) {
                    
                    # Shift xyz coords with centroid_2d (x,y,z) and crown_middle (or high or low) being th origin
                    xyz <- crown_slices[[c]] %>% 
                        mutate(X = centroid_2d[1] - X,
                                Y = centroid_2d[2] - Y,
                                Z = Z)
                
                ## Transform coordinates into spherical coordinates ##
                    ptphi <- xyz %>% 
                                mutate(#p = sqrt((X^2)+(Y^2)+(Z^2)),
                                        theta = atan(Y/X), # atan returns value in radians
                                        #phi = acos(Z/(sqrt((X^2)+(Y^2)+(Z^2)))), # acos returns value in radians. And phi ranges from 0 to pi
                                        r = sqrt((X^2)+(Y^2))) # r is the distance (in xy plane) from the point to the Z axis of origin

                    # bin angles every 5 degrees (i.e. 5*(pi/180))
                    angles <- c(0,5,10,15,20,25,30,35,40,45,50,55,60,65,70,75,80,85,90)
                    max_r <- c()

                    # get max 'r' (xy-plane distance from point to Z-axis of origin) for each angle category
                    for (a in 1:(length(angles)-1)) {
                        df <- subset(ptphi, theta > (angles[a]*(pi/180)) & theta <= (angles[a+1]*(pi/180)))
                        rs <- max(df$r, na.rm=TRUE)
                        max_r <- append(max_r, rs)
                        }

                    # calculate the mean of all the max radii from each phi category
                    if (any(is.infinite(max_r)) == TRUE) {
                            max_r[!is.finite(max_r)] <- NA # turn -Inf into NA
                            diameter <- 2*(mean(max_r, na.rm=TRUE)) # multiply by 2 to get the diameter
                        } else {
                            diameter <- 2*(mean(max_r, na.rm=TRUE))
                        }

                    # save the diameter of that slice
                    crown_diameters[[c]] <- diameter
                }
    
        names(crown_diameters) <- c("top_crown_diameter", "middle_crown_diameter", "bottom_crown_diameter" )

    return(crown_diameters)
}


#' @title Trunk width/diameter
#' 
#' @description Trunk width is calculated by calculating crown radius at 5° intervals at different heights
#' along the trunk, computing the mean, and multiplying this mean by 2.
#' 
#' @details Equivalent to the crown width calculation but measured 10 times per trunk at 10% height intervals. 
#' The mean trunk width across all height intervals is taken. 
#' 
#' @param trunk_pc lidR object. A point cloud file of an individual tree trunk.
#' Designed for single-spectral TLS (terrestrial laser scanning) or 
#' high-resolution MLS (mobile laser scanning) data.
#' Tree trunk point cloud can be obtained from the classify_crown function.
#' @param slice_radius Numeric. Half the height in meters of the slices used to calculate the trunk width. 
#' Default is 0.1 meters.
#' @param top_height Numeric. Proportional crown height (values range from 0 to 1) at which the width of the top of the trunk
#' will be calculated. Default is 0.6.
#' @param middle_height Numeric. Proportional crown height (values range from 0 to 1) at which the width of the middle of the trunk
#' will be calculated. Default is 0.5.
#' @param bottom_height Numeric. Proportional crown height (values range from 0 to 1) at which the width of the bottom of the trunk
#' will be calculated. Default is 0.4.
#' @return Numeric. Trunk width (units match units of the point cloud X, Y, Z coordinates)
#' @export

trunk_width <- function(trunk_pc, slice_radius = 0.1, top_height = 0.6, middle_height = 0.5, bottom_height = 0.4) {
    # trunk_pc = lidR object. Point cloud of the tree's trunk. Can be obtained by running classify_crown. X, Y, Z coordinates hsould be in meters.
    # slice_radius = Height in meters to be included below and above the slicing point. Default is 0.1 meters.
    # top_height = numeric. Proportional trunk height (values range from 0 to 1) at which the uppermost trunk width will be calculated
    # middle_height = numeric. Proportional trunk height (values range from 0 to 1) at which the middle-most trunk width will be calculated
    # bottom_height = numeric. Proportional crown height (values range from 0 to 1) at which the bottommost trunk width will be calculated

    # Find heights at middle of trunk
    trunk_pc <- data.frame(X = trunk_pc$X, Y = trunk_pc$Y, Z = trunk_pc$Z) # make lidR object trunk_pc into xyz dataframe

        trunk_high <- ((max(trunk_pc$Z) - min(trunk_pc$Z))*top_height) + min(trunk_pc$Z)
        trunk_middle <- ((max(trunk_pc$Z) - min(trunk_pc$Z))*middle_height) + min(trunk_pc$Z)
        trunk_low <- ((max(trunk_pc$Z) - min(trunk_pc$Z))*bottom_height) + min(trunk_pc$Z)

        trunk_heights <- list(trunk_high, trunk_middle, trunk_low)

        # determine trunk center axis (x,y coordinate)
        knn <- FNN::get.knn(trunk_pc, k = 10) # calculate knn distance for each point
        knn_means <- rowMeans(knn$nn.dist, na.rm = TRUE) # average the knn distance for each point (which I will use as weights to find the centroid)
        weights <- 1 - (knn_means - min(knn_means, na.rm = TRUE)) / (max(knn_means, na.rm = TRUE)-min(knn_means, na.rm = TRUE)) # inverse these weights (so that less distance carries a higher weight) and proportionalize them

        # find centroid of trunk in xy plane
        centroid_X <- sum(trunk_pc$X * weights, na.rm = TRUE) / sum(weights, na.rm = TRUE)
        centroid_Y <- sum(trunk_pc$Y * weights, na.rm = TRUE) / sum(weights, na.rm = TRUE)
        #centroid_2d <- c(mean(trunk_pc$X), mean(trunk_pc$Y)) # NON weighted centroid
        centroid_2d <- c(centroid_X, centroid_Y)

        # Create slices at each trunk height
        # get points with Z +/- slice thickness for middle, high, and low
        trunk_high_slice <- subset(trunk_pc, Z<(trunk_high+slice_radius) & Z>(trunk_high-slice_radius)) # note that the xyz coords are in meters
        trunk_middle_slice <- subset(trunk_pc, Z<(trunk_middle+slice_radius) & Z>(trunk_middle-slice_radius))
        trunk_low_slice <- subset(trunk_pc, Z<(trunk_low+slice_radius) & Z>(trunk_low-slice_radius))

        trunk_slices <- list(trunk_high_slice, trunk_middle_slice, trunk_low_slice)

        # For each crown slice:
        trunk_diameters <- list()

        for (c in 1:length(trunk_slices)) {
                    
            # Shift xyz coords with centroid_2d (x,y,z) and crown_middle (or high or low) being th origin
            xyz <- trunk_slices[[c]] %>% 
                    dplyr::mutate(X = centroid_2d[1] - X,
                            Y = centroid_2d[2] - Y,
                            Z = Z)

            ## Spherical coords ##
            ptphi <- xyz %>% 
                        dplyr::mutate(#p = sqrt((X^2)+(Y^2)+(Z^2)), # for spherical coords
                        theta = atan(Y/X), # atan returns value in radians
                        #phi = acos(Z/(sqrt((X^2)+(Y^2)+(Z^2)))), # acos returns value in radians. And phi ranges from 0 to pi
                        r = sqrt((X^2)+(Y^2)))

            # bin angles every 5 degrees (i.e. 5*(pi/180))
            angles <- c(0,5,10,15,20,25,30,35,40,45,50,55,60,65,70,75,80,85,90)
            max_r <- c()

            # get max 'r' (xy-plane distance from point to Z-axis of origin) for each angle category
            for (a in 1:(length(angles)-1)) {
                    df <- subset(ptphi, theta > (angles[a]*(pi/180)) & theta <= (angles[a+1]*(pi/180)))
                    rs <- max(df$r, na.rm=TRcrown_width <- function(crown_pc, slice_radius = 0.1, top_height = 0.9, middle_height = 0.5, bottom_height = 0.1)UE)
                    max_r <- append(max_r, rs)
                }

            # calculate the mean of all the max radii from each phi category
            if (any(is.infinite(max_r)) == TRUE) {
                    max_r[!is.finite(max_r)] <- NA # turn -Inf into NA
                    diameter <- 2*(mean(max_r, na.rm=TRUE)) # multiply by 2 to get the diameter
                } else {
                    diameter <- 2*(mean(max_r, na.rm=TRUE))
                }

            # save the diameter of that slice
            trunk_diameters[[c]] <- diameter
                }

        # trunk diameter average
        trunk_diameter <- mean(unlist(trunk_diameters))
    
    return(trunk_diameter)
}


#' @title Crown shape
#' 
#' @description Calculates crown shape based on the crown widths at the top, middle, and bottom of the crown 
#' (as computed by the function crown_width). Returns crown shape as one of four options: conal (where bottom 
#' crown width is the smallest, followed by the middle and the top), vase (where bottom crown width is the 
#' largest, followed by the middle and the top), column (where all widths are of the same length, allowing 
#' for a buffer zone), and oval (where the middle width is the largest).
#' Since there are only four shapes, this is a rough estimate (i.e. oval and round shapes are essentially
#' the same here).
#' The function therefore first computes crown_widths then determines crown shape from these values.
#' 
#' @param crown_pc lidR object. A point cloud file of an individual tree crown.
#' Designed for single-spectral TLS (terrestrial laser scanning) or 
#' high-resolution MLS (mobile laser scanning) data.
#' Tree crown point cloud can be obtained from the classify_crown function.
#' @param buffer Numeric. Proportion of the crown widths (value ranges from 0 to 1). The permissable buffer zone before
#' the function determines whether one width is smaller of larger than another.
#' @param slice_radius Numeric. Half the height in meters of the slices used to calculate crown widths. 
#' Default is 0.1 meters.
#' @param top_height Numeric. Proportional crown height (values range from 0 to 1) at which the width of the top of the crown 
#' will be calculated. Default is 0.9.
#' @param middle_height Numeric. Proportional crown height (values range from 0 to 1) at which the width of the middle of the crown
#' will be calculated. Default is 0.5.
#' @param bottom_height Numeric. Proportional crown height (values range from 0 to 1) at which the width of the bottom of the crown
#' will be calculated. Default is 0.1.
#' @return Character.
#' @export 

crown_shape <- function(crown_pc, buffer = 0.1, slice_radius = 0.1, top_height = 0.9, middle_height = 0.5, bottom_height = 0.1) {
    widths <- crown_width(crown_pc = crown_pc, slice_radius = slice_radius, top_height = top_height, middle_height = middle_height,
                            bottom_height = bottom_height)

    if ((widths[[3]]-buffer > widths[[2]]) # top crown width is bigger than middle crown
        & (widths[[3]]+buffer > widths[[2]]) # &
        & (widths[[2]]-buffer < widths[[1]]) # middle crown width is bigger than bottom crown
        & (widths[[2]]+buffer < widths[[1]])) {
        shape <- "vase"
    } else if ((widths[[1]]-buffer > widths[[2]]) # bottom crown width is bigger than middle crown
        & (widths[[1]]+buffer > widths[[2]]) # &
        & (widths[[2]]-buffer < widths[[3]]) # middle crown width is bigger than top crown
        & (widths[[2]]+buffer < widths[[3]])) {
        shape <- "cone"
    } else if ((widths[[3]]-buffer <= widths[[2]]) # top crown width is same size as middle crown
        & (widths[[3]]+buffer >= widths[[2]]) # &
        & (widths[[2]]-buffer <= widths[[1]]) # middle crown width is same size as bottom crown
        & (widths[[2]]+buffer >= widths[[1]])) {
        shape <- "column"
    } else ((widths[[3]]-buffer < widths[[2]]) # top crown width is smaller than middle crown
        & (widths[[3]]+buffer < widths[[2]]) # &
        & (widths[[2]]-buffer > widths[[1]]) # middle crown width is larger than bottom crown
        & (widths[[2]]+buffer > widths[[1]])) 
        {
        shape <- "oval"
    }

    return(shape)
}