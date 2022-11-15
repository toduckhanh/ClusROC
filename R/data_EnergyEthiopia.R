#' @title A subset of energy choice data in 4 cities of Ethiopia
#' @name EnergyEthiopia
#' @docType data
#'
#' @description A subset of energy choice data used in Alem et al. (2016). The authors are used the full dataset to investigate the determinants of household cooking fuel choice and energy transition in urban Ethiopia. A full data is publicly available at \url{https://doi.org/10.1016/j.eneco.2016.06.025}.
#'
#' @usage data(EnergyEthiopia)
#'
#' @format A data frame with 2088 observations from 1123 households (or clusters) in the capital Addis Ababa and 9 variables:
#' \describe{
#'   \item{\code{uqid}}{the id of household (which yield 1123 clusters).}
#'   \item{\code{energy2}}{a factor with 3 levels (types) of cooking energy state at each time (2000, 2004, 2009), i.e., 1 (clean fuel only - electricity, gas and kerosene), 2 (a mix of clean and biomass), 3 (biomass fuel only - firewood, charcoal, dung and crop residues).}
#'   \item{\code{hhs}}{household size.}
#'   \item{\code{hhs_ft}}{a factor with 4 levels of household size: small (1 \eqn{\le} hhs \eqn{\le} 4); medium (5 \eqn{\le} hhs \eqn{\le} 8); large ((9 \eqn{\le} hhs \eqn{\le} 12)); very large (hhs \eqn{\ge} 13).}
#'   \item{\code{lrconsaeu}}{log of real consumption per adult equivalent units.}
#'   \item{\code{lfirewood_pr}}{Firewood log price.}
#'   \item{\code{lcharcol_pr}}{Charcoal log price.}
#'   \item{\code{lkerosene_pr}}{Kernosene log price.}
#'   \item{\code{lelectric_pr}}{Electricity log price.}
#'}
#'
#' @references
#' Alem, Y., Beyene, A. D., Köhlin, G., & Mekonnen, A. (2016). "Modeling household cooking fuel choice: A panel multinomial logit approach". \emph{Energy Economics}, \bold{59}, 129-137.
#'
#' @keywords data
#'
NULL
