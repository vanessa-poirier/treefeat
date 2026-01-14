# Treefeat
An R package for computing tree architectural features from LiDAR point clouds

_Author: Vanessa Poirier_

_Release Date: January 2026_

# Description
Compute individual tree architectural/structural features from a LiDAR point cloud.
 Functions in this package can accept a LiDAR point cloud (las/laz file or a loaded lidR object)
 of an individual tree and compute summary features (e.g. tree height/width, crown height/width,
 trunk height/width, and crown projected area), architecural features related to overall shape
 (e.g. ratios between tree/trunk and crown heights/widths, and crown shapes), and crown complexity
 features (e.g. gap fraction, crown fullness, and crown curvature).

While each function can be used individually to compute 1 or more features, the function
 extract_features calculates all available features at once to provide a summary of the tree's
 architecture.
 
Another notable function is classify_crown which separates the point cloud into crown and trunk,
 this step is required for most other functions in the package.
 This package is robust against lower resolution TLS (terrestrial laser scanning) or MLS (mobile
 laser scanning) point clouds. While other available packages are adapted for TLS (often necessitating
 quantitative structural models (QSMs) or alpha shapes which are not always possible to compute for
 lower resolution point clouds) or ALS (aerial LiDAR), this package is adapted for MLS data.
 Of course, higher resolution point clouds will give more reliable calculations.
