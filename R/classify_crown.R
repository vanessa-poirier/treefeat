# This script is to classify a tree point cloud into crown and non-crown points based on a KNN and time-series algorithm
# returns a list of two point clouds (lidR objects): the first being the crown points, the second being the trunk points (or non-crown points)
# dependent libraries: lidR, dplyr, FNN

#' @title Classify tree crown and trunk points in a point cloud
#' 
#' @description Separates a tree point cloud into the crown and trunk.
#' 
#' @param file las/laz file or lidR object. A point cloud file of an individual tree
#' Designed for single-spectral TLS (terrestrial laser scanning) or 
#' high-resolution MLS (mobile laser scanning) data.
#' @param res Numeric. Resolution/intervals (i.e. height in meters) at which to slice the point cloud.
#' @param denoise Logical. Determines whether point cloud noise should be removed (determined by lidR::classify_noise).
#' @param max_crown_height Numeric. The maximum height at which the crown can start (presented as proportional height). 
#' Default is 0.4 (i.e. the crown must start =< 40% up the tree).
#' @return list of two point clouds, the first being the tree crown and the second being the trunk
#' @export

classify_crown <- function(file, res = 0.1, denoise = TRUE, max_crown_height = 0.4) {
    # res = resolution/intervals (i.e. height in meters) at which to slice the point cloud
    # denoise = logical. Whether point cloud noise should be removed (determined by lidR::classify_noise)
    # max_crown_height = the maximum height at which the crown can start (presented as proportional height). 
    #        Default is 0.4 (i.e. the crown must start =< 40% up the tree)

    ## Load and prepare point cloud file
    if (endsWith(file, ".las") | endsWith(file, ".laz")) {
        tree_pc <- lidR::readTLSLAS(file) # load point cloud
    } else if (class(file)[1] == "LAS") {
        tree_pc <- file
    } else {
        stop("Invalid file_type specified. Use las/laz file or lidR object.")
    }

    if (denoise == TRUE) {
        tree_pc <- lidR::classify_noise(tree_pc, lidR::sor(k = 10, m = 1, quantile = FALSE)) # remove noise from point cloud
        tree_pc <- lidR::filter_poi(tree_pc, Classification==18) # filter to keep only non-noise points (classed as '18')
    } else {
        tree_pc <- tree_pc
    }

    # Create an xyz dataframe from the point cloud
    tree_df <- data.frame(X = tree_pc$X, Y = tree_pc$Y, Z = tree_pc$Z)

    # Slice point cloud every 'res' meters (default is 10 cm)
    slices <- seq(min(tree_df$Z), max(tree_df$Z), by = res) # Create a vector of slice break points
    tree_df$slice_id <- findInterval(tree_df$Z, slices) # Assign a slice ID to each point


    ## Calculate KNN dispersion of points in each slice ##
        # Compute KNN-Based Dispersion Per Slice
    get_k <- function(n) { # n is the number of points in a given slice
        if (n <= 100) {
            k = n - 1 # k must be number of points - 1
        } else {
            k = 100
        }
        return(k)
    }  # Number of nearest neighbors on a sliding scale --> min value is 'number of points - 1' and max is '100'

    # create mean_knn_distance function to calculate mean knn distance per slice
    calculate_mean_knn_distance = function(x, k) {
        coords <- data.frame(X=x$X, Y=x$Y, Z=x$Z)
        knn <- FNN::get.knn(coords, k = k)
        mean_knn_dist <- mean(rowMeans(knn$nn.dist))
        return(mean_knn_dist)
    }

    dispersion_by_slice <- tree_df %>%
        dplyr::group_by(slice_id) %>%
        dplyr::filter(n() > 10) %>%  # Ensure that each slice has at least 10 values (enough to run a meaningful KNN)
        dplyr::summarise(
            slice_id = unique(slice_id),
            k = get_k(n()),
            mean_knn_distance = calculate_mean_knn_distance(.data, k) # .data calls to the data subsetted by the beforehand grouping and filtering
            )

        # Find Sudden Changes in Dispersion
    dispersion_by_slice <- dispersion_by_slice %>%
    dplyr::mutate(
        dispersion_change = c(NA, abs(diff(mean_knn_distance)))
    )

        # at which height is there a maximum difference in dispersion from the previous height?
        # based on average knn dispersion per slice
        # ignore lower and upper portion of tree (otherwise the algorithm has a hard time finding the crown start)
    #lower_portion <- 1 # option to ignore lower part of tree and only start considering crown start at > 10cm high
    upper_portion <- ceiling((max(dispersion_by_slice$slice_id) - 
                        min(dispersion_by_slice$slice_id))*max_crown_height + min(dispersion_by_slice$slice_id)) 
                        # finds slice ID that is 40% up the tree
    
    crown_start_slice <- dispersion_by_slice$slice_id[which(dispersion_by_slice$dispersion_change == 
                                                                max(dispersion_by_slice$dispersion_change[-c(1,which(dispersion_by_slice$slice_id > upper_portion))]))] # ignore the upper portion of tree

    crown_start <- min(subset(tree_df, slice_id == crown_start_slice)$Z)

    ## Classify crown points ##
    # crown points will be all points with Z > the crown_start
    crown_pc <- lidR::filter_poi(tree_pc, Z >= crown_start)
    trunk_pc <- lidR::filter_poi(tree_pc, Z < crown_start)

    return(list(crown = crown_pc, trunk = trunk_pc))
}