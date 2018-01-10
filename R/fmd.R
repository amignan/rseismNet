#' Compute the earthquake FMD
#'
#' @param m Earthquake magnitude vector
#' @param mbin Magnitude binning number
#' @return The earthquake frequency magnitude distribution (FMD)
#' @examples
#' beta <- log(10); mc <- 2; mbin <- 0.1
#' m <- mc - mbin / 2 + rexp(100, beta)
#' mdistr <- fmd(m, mbin)
#' plot(mdistr$mi, mdistr$Ni, log = "y")
#' points(mdistr$mi, mdistr$ni)
fmd <- function(m, mbin) {
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
#' @param m Earthquake magnitude vector
#' @param mbin Magnitude binning number
#' @param method Method following which the completeness magnitude mc is computed.
#' Options include "mode" and "KStest"
#' @return The completeness magnitude mc
#' @examples
#' m <- 
#' mc <- mc.val(m, 0.1, "mode")
mc.val <- function(m, mbin, method) {
  if (method == "mode") {
    dens <- density(m, from = min(m) - 1, to = max(m) + 1)
    return(round(dens$x[which(dens$y == max(dens$y))], digits = log10(1 / mbin))[1])
  }
  if (method == "KStest") {
    mdistr <- fmd(m, mbin)
    mmin <- mc.val(m, mbin, "mode")
    mscan <- seq(mmin, max(m)-mbin, mbin)
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
#' @param m Earthquake magnitude vector
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


