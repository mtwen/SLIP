#' fMRI data used in our paper
#'
#' A dataset contain the regions and their region-averaged hemodynamic response
#'
#' @format A list with 6 compoenents
#' \describe{
#'   \item{dat}{a 360 x 264 data matrix (360 scans, 264 regions)}
#'   \item{regions}{names of regions. For example, "2001xlylxl", "2001" indicates
#'   the Tzourio-Mazoyer code in the Anatomical Automatic Labeling (AAL) partition
#'   and "xlylxl" indicates the sub-region of the certain ROI.}
#'   \item{coord}{a data.frame with four columns: the first three columns indicate
#'   the (x, y, z) location in 3D fMRI image and the last column indicates the
#'   corresponding region.}
#'   \item{dimx}{dimension of x axis of 3D fMRI image.}
#'   \item{dimy}{dimension of y axis of 3D fMRI image.}
#'   \item{dimz}{dimension of z axis of 3D fMRI image.}
#' }
#' @source original data from the subject P2 on
#' \url{http://www.cs.cmu.edu/~fmri/science2008/data.html}
#' @references "Predicting Human Brain Activity Associated with the Meanings of Nouns,"
#' T. M. Mitchell, S. V. Shinkareva, A. Carlson, K.M. Chang, V. L. Malave, R. A. Mason,
#' and M. A. Just, Science, 320, 1191, May 30, 2008. DOI: 10.1126/science.1152876.

"fmri.data"
