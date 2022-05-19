#' The Error Function (or Gaussian Error function).
#'
#' For any \eqn{x \geq 0}, this returns \eqn{\operatorname{Erf}(x)}. Here the
#' function \eqn{\operatorname{Erf}(x)} is defined to be
#' \deqn{\operatorname{Erf}(x) := \frac{2}{\sqrt{\pi}} \int_{0}^{x} e^{-t^{2}} d t}
#'
#' @param x : A non-negative real number.
#'
#' @return : A non-negative real value corresponding to the error function.
#' @export
#'
#' @examples
#' \dontrun{
#' erf(0) # should return 0
#' }
erf <- function(x) {
  base::return(2 * stats::pnorm(x * base::sqrt(2)) - 1)
}

# TODO: Check if this is the mean of the Folded normal, not the inverse folded normal mean!

#' The Inverse Folded Normal mean value, for a fixed standard deviation.
#'
#' @param theta : A non-negative real number.
#' @param sigma : A positive real number.
#'
#' @return
#' @export
mean_theta_sigma <- function(theta, sigma) {
  # TODO: Assertion check theta non-negative
  # TODO: Assertion check sigma positive. Can we have sigma = 0?
  x <- ((theta / sigma)**2 / 2)
  base::return((sqrt(pi) * sqrt(x) * fAsianOptions::erf(sqrt(x)) + exp(-x)) * sqrt(2 / pi) * sigma)
}


