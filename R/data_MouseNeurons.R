#' @title A subset of mouse brain cells data
#' @name MouseNeurons
#' @docType data
#'
#' @description A subset of mouse brain cells data used in To el al. (2022). This is used to evaluate the ability of Lamp5 gene to discriminate three types of glutamatergic neurons. A full data is publicly available at \url{https://portal.brain-map.org/atlases-and-data/rnaseq/mouse-v1-and-alm-smart-seq}.
#'
#' @usage data(MouseNeurons)
#'
#' @format A data frame with 860 observations from 23 clusters and 7 variables:
#' \describe{
#'   \item{\code{sample_name}}{name of each observation.}
#'   \item{\code{subclass_label}}{a factor with 3 levels (types) of glutamatergic neurons, i.e., L2/3 IT (Layer 2/3 Intratelencephalic), L4 (Layer 4) and L5 PT (Layer 5 Pyramidal Tract) neurons.}
#'   \item{\code{genotype_id}}{the mouse genotype (which yield 23 clusters).}
#'   \item{\code{sex}}{the gender of mouse.}
#'   \item{\code{age_days}}{the age of mouse, in days.}
#'   \item{\code{Slc17a7_cpm}}{count per million of Slc17a7 (Solute Carrier Family 17 Member 7) gene expression.}
#'   \item{\code{Lamp5_cpm}}{count per million of Lamp5 (Lysosomal Associated Membrane Protein Family Member 5) gene expression.}
#'}
#'
#' @references
#' To, D-K., Adimari, G., Chiogna, M. and Risso, D. (2022)
#' ``Receiver operating characteristic estimation and threshold selection criteria in three-class classification problems for clustered data''. \emph{Statistical Methods in Medical Research}, \bold{7}, 31, 1325-1341.
#'
#' @keywords data
#'
NULL
