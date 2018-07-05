
#' Circular ggplot scales
#'
#' Extensions of \code{scale_x_continuous} and \code{scale_y_continuous} that
#' are useful for circular axes.
#'
#' Behind the screens, radians will be used, but the plot will display
#' transformed values.
#'
#' The continuous transformations are \code{"radians"}, \code{"degrees"} and
#' \code{"hours"}. For these, \code{nticks} and \code{digits} can be set.
#'
#' The categorical transformations are \code{"cardinal"}, which prints character
#' labels of the directions (left, right, up, down), and \code{"compass"}, which
#' prints compass directions (North, South, East, West). Both of these assume
#' the circular data starts on the right and moves counterclockwise. If not,
#' appropriate transformations must first be taken, or the labels will be
#' non-sensical.
#'
#' Two additional labels are \code{"texpi"} and \code{"texnegpi"}, which print
#' nice labels for TeX to interpret starting at 0 and \code{-pi} respectively.
#' This is useful if TikZ is used.
#'
#' @param units The units to display on the axis. See `Details` for the options.
#' @param nticks The number of ticks to display. Only relevant for continous
#'   scales.
#' @param digits Numeric; The number of digits for continuous scales.
#' @param limits Numeric vector; The limits of the plot. Must two numbers that
#'   are  \code{2*pi} apart.
#' @param scale_function The scale function from ggplot to use. Either
#'   \code{scale_x_continuous} or \code{scale_y_continuous}. For convenience
#'   this option can be circumvented by directly using \code{scale_x_circular}
#'   or \code{scale_y_circular}.
#' @param ... Additional arguments to be passed to  \code{scale_x_continuous} or
#'   \code{scale_y_continuous}.
#'
#' @return An object of type \code{ScaleContinuousPosition} that can be added to
#'   any existing \code{ggplot}.
#' @export
#'
#' @examples
#' p <- plot_batmixfit(params = cbind(mu = c(-pi/2, 0, pi/2, pi), kp = 4,
#'                               lam = c(-.9, .2, .8, 0), alph = .25))
#'
#' p + scale_x_circular(units = "compass")
#' p + scale_x_circular(units = "compass", limits = c(-1, 2*pi - 1))
#'
#' p + scale_x_circular(units = "cardinal")
#'
#' p + scale_x_circular(units = "texpi")
#' p + scale_x_circular(units = "texnegpi")
#'
#' p + scale_x_circular(units = "hours", nticks = 24)
#'
scale_circular <- function(units = "degrees", nticks = 4,
                           digits = 0, limits = c(0, 2 * pi),
                           scale_function = ggplot2::scale_x_continuous, ...) {

  if (round(abs(limits[1] - limits[2]) - 2*pi, 3) != 0) {
    stop("Limits must be a range of length 2*pi.")
  }

  units_categ <- units %in% c("texpi", "texnegpi", "cardinal", "compass")

  if (units_categ & nticks != 4) {
    warning(paste0("For `units = ", units, "`, setting `nticks` to 4."))
    nticks <- 4
  }

  if (units_categ) {
    if (units == "texpi") {
      brks   <- seq(0, 2 * pi, length.out = 5)
      limits <- c(0, 2 * pi)
    } else if (units == "texnegpi") {
      brks   <- seq(-pi, pi, length.out = 5)
      limits <- c(-pi, pi)
    } else {
      # These are the options for cardinal or compass directions.
      brks   <- seq(-pi, 3*pi, length.out = 9)
    }
  } else {
    brks <- seq(limits[1], limits[2], length.out = nticks + 1)
  }

  conv_brks <- switch(
    units,
    radians  = round(brks, digits),
    degrees  = round(brks * 180 / pi, digits),
    hours    = round(brks * 12 / pi, digits),
    texpi    = c("$0$", "$\\frac{\\pi}{2}$", "$\\pi$",
                 "$\\frac{3\\pi}{2}$", "$2\\pi$"),
    texnegpi = c("$-\\pi$", "$-\\frac{\\pi}{2}$", "$0$",
                 "$\\frac{\\pi}{2}$", "$\\pi$"),
    cardinal = c("Right", "Up", "Left", "Down")[c(3:4, 1:4, 1:3)],
    compass  = c("East", "North", "West", "South")[c(3:4, 1:4, 1:3)])

  scale_function(breaks = brks,
                 labels = conv_brks,
                 limits = limits,
                 ...)
}

#' @export
#' @describeIn scale_circular
scale_x_circular <- function(...) {
  scale_circular(scale_function = ggplot2::scale_x_continuous, ...)
}

#' @export
#' @describeIn scale_circular
scale_y_circular <- function(...) {
  scale_circular(scale_function = ggplot2::scale_y_continuous, ...)
}



