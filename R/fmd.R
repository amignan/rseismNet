#' Earthquake FMD Computation
#'
#' Compute the earthquake frequency magnitude distribution (FMD) for a vector of
#' earthquake magnitudes \code{m} with binning \code{mbin}.
#'
#' @param m a numeric vector of earthquake magnitudes
#' @param mbin the magnitude binning value
#' @return The earthquake FMD data frame made of the 3 columns:
#'    \item{mi}{the magnitude bins}
#'    \item{ni}{the non-cumulative number of earthquakes of magnitude \code{mi}}
#'    \item{Ni}{the cumulative number of earthquakes of magnitude \eqn{\ge} \code{mi}}
#' @examples
#' beta <- log(10); mc <- 2
#' m <- mc - mbin / 2 + rexp(100, beta)
#' mdistr <- fmd(m)
#' View(mdistr)
fmd <- function(m, mbin = 0.1) {
  mmin <- min(m, na.rm = T)
  mmax <- max(m, na.rm = T)
  mi <- seq(floor(mmin / mbin) * mbin, ceiling(mmax / mbin) * mbin, mbin)
  ni <- hist(m, breaks = c(mi, mmax + mbin) - mbin / 2, plot = F)$counts
  Ni <- sapply(1:length(ni), function(i) sum(ni[i:length(ni)]))
  return(data.frame(mi = mi, ni = ni, Ni = Ni))
}

#' Estimate the completeness magnitude
#'
#' Estimate the completeness magnitude mc from the distribution of the magnitude m
#' using two different methods, the mode or a KS test of the Gutenberg-Richter law.
#'
#' @param m a numeric vector of earthquake magnitudes
#' @param mbin Magnitude binning number
#' @param method Method following which the completeness magnitude mc is computed.
#' Options include "mode" and "KStest"
#' @return The completeness magnitude mc
#' @examples
#' theta <- list(kappa = 3 * log(10), beta = log(10), mc = 2)
#' m.angular <- efmd.sim(1e3, theta)
#' mc.val(m.angular, "mode")          #correctly estimate mc
#' mc.val(m.angular, "KStest")        #can significantly over-estimate mc
#'
#' theta <-
#' m.curved <- bfmd(1e3, theta)
#' mc.val(m.angular, "mode")          #can significantly under-estimate mc
#' mc.val(m.curved, "KStest")         #tendency to over-estimate mc
mc.val <- function(m, method, mbin = 0.1) {
  if (method == "mode") {
    dens <- density(m, from = min(m) - 1, to = max(m) + 1)
    return(round(dens$x[which(dens$y == max(dens$y))], digits = log10(1 / mbin))[1])
  }
  if (method == "KStest") {
    mdistr <- fmd(m, mbin)
    mmin <- mc.val(m, "mode", mbin)
    mscan <- seq(mmin, max(m) - mbin, mbin)
    nscan <- length(mscan)
    betai <- sapply(1:nscan, function(i) beta.mle(m[ which(m > mscan[i] - mbin / 2)],
                                                  mscan[i], mbin))
    pred <- sapply(1:nscan, function(i) exp(-betai[i]*(mscan[i:nscan]-mscan[i])))
    obs <- sapply(1:nscan, function(i) mdistr$n[which(mdistr$mi >= mscan[i])])
    D <- sapply(1:nscan, function(i) max(abs(cumsum(-diff(unlist(pred[i]))) -
                                             cumsum(-diff(unlist(obs[i]))/
                                                      max(unlist(obs[i]))) )))
    return(mscan[which(D == min(D[is.finite(D)]))][1])
  }
}

#' Compute the FMD beta-value
#'
#' @param m a numeric vector of earthquake magnitudes
#' @param mc Completeness magnitude number
#' @param mbin Magnitude binning number
#' @return The value of beta for a given m vector
#' @examples
#' beta <- log(10); mc <- 2; mbin <- 0.1
#' m <- round(mc - mbin / 2 + rexp(100, beta), digits = log10(1/mbin))
#' get_beta(m, mc, mbin)
beta.mle <- function(m, mc, mbin) {
  1 / (mean(m[which(m > mc - mbin / 2)]) - (mc - mbin / 2))
}

#' Compute the pdf of the elemental FMD
#'
#' @param m a numeric vector of earthquake magnitudes
#' @param mc the ompleteness magnitude number
#' @param mbin Magnitude binning number
#' @return The value of beta for a given m vector
#' @examples
#' ...
efmd.pdf <- function(m, theta) {
  indcomplete <- which(m >= theta$mc)
  n.incomplete <- exp((theta$kappa - theta$beta) * (m[-indcomplete] - theta$mc))
  n.complete <- exp(-theta$beta * (m[indcomplete] - theta$mc))
  f_angular <- c(n.incomplete, n.complete)
  c <- 1 / (1 / (theta$kappa - theta$beta) + 1 / theta$beta)
  return(c * f_angular)
}

#' Elemental FMD Simulation
#'
#' Simulate the elemental earthquake frequency magnitude distribution (eFMD) by
#' applying the Inverse Method to the angular FMD model (Mignan, 2012)
#' with probability density function \code{efmd.pdf}.
#'
#' If \code{mbin} is not specified it assumes the default value of \code{0.1}.
#'
#' The angular FMD model has density
#'
#' @param N number of earthquakes
#' @param mc Completeness magnitude number
#' @param mbin Magnitude binning number
#' @return A numeric vector of \code{N} earthquake magnitudes.
#' @references Mignan, A. (2012), Functional shape of the earthquake frequency-magnitude
#' distribution and completeness magnitude, J. Geophys. Res., 117, B08302,
#' \href{http://onlinelibrary.wiley.com/doi/10.1029/2012JB009347/full}{doi: 10.1029/2012JB009347}
#' @seealso \code{efmd.pdf} for the probability density function of the eFMD.
#' @examples
#' theta <- list(kappa = 3 * log(10), beta = log(10), mc = 2)
#' m.sim <- efmd.sim(1e3, theta)
#' mdistr <- fmd(m.sim)
#' plot(mdistr$mi, mdistr$Ni, log = "y")
#' points(mdistr$mi, mdistr$ni)
efmd.sim <- function(N, theta, mbin = 0.1) {
  N.complete <- round(N * (theta$kappa - theta$beta) / theta$kappa)
  m.complete <- theta$mc - mbin / 2 + rexp(N.complete, theta$beta)
  m.incomplete <- theta$mc - mbin/2 -rexp(N - N.complete, (theta$kappa - theta$beta))
  return(round(c(m.incomplete, m.complete), digits = log10(1 / mbin)))
}

