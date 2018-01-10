#' Compute the earthquake FMD
#'
#' @param m Earthquake magnitude vector
#' @param mbin Magnitude binning number
#' @return The earthquake frequency magnitude distribution (FMD)
#' @examples
#' beta <- log(10); mc <- 2; mbin <- 0.1
#' m <- mc - mbin / 2 + rexp(100, beta)
#' fmd(m, mbin)
get_fmd <- function(m, mbin) {
  mmin <- min(m, na.rm = T)
  mmax <- max(m, na.rm = T)
  mi <- seq(floor(mmin / mbin) * mbin, ceiling(mmax / mbin) * mbin, mbin)
  ni <- hist(m, breaks = c(mi, mmax + mbin) - mbin / 2)$counts
  Ni <- sapply(1:length(ni), function(i) sum(ni[i:length(ni)]))
  return(data.frame(m=mi, n=ni, N=Ni))
}

#' Add together two numbers.
#'
#' @param x A number.
#' @param y A number.
#' @return The sum of \code{x} and \code{y}.
#' @examples
#' add(1, 1)
#' add(10, 1)
get_mc <- function(m, shape, mbin) {
  if(shape == "angular") mc_angular(m, mbin)
  if(shape == "curved") mc_curved(m, mbin)
}

#' Compute the FMD beta-value
#'
#' @param m Earthquake magnitude vector
#' @param mc Completeness magnitude number
#' @param mbin Magnitude binning number
#' @return The value of beta for a given m vector
#' @examples
#' beta <- log(10); mc <- 2; mbin <- 0.1
#' m <- round(mc - mbin / 2 + rexp(100, beta), digits = log10(1/mbin))
#' get_beta(m, mc, mbin)
get_beta <- function(m, mc, mbin) {
  1 / (mean(m[which(m > mc - mbin/2)]) - (mc - mbin / 2))
}


