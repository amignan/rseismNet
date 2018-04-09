#' The Earthquake FMD
#'
#' Compute the earthquake frequency magnitude distribution (FMD) for a vector of
#' earthquake magnitudes \code{m} with binning \code{mbin}.
#'
#' Magnitude intervals have the form
#' \out{(<i>m<sub>i</sub></i> - <i>m<sub>bin</sub></i>/2, <i>m<sub>i</sub></i> + <i>m<sub>bin</sub></i>/2]},
#' as defined in function \code{graphics::hist} (for its default \code{right = TRUE}).
#'
#' @param m a numeric vector of earthquake magnitudes
#' @param mbin the magnitude binning value (if not provided, \code{mbin = 0.1})
#' @return The earthquake FMD data frame of 3 parameters:
#'    \item{mi}{the magnitude bins}
#'    \item{ni}{the non-cumulative number of earthquakes of magnitude \code{mi}}
#'    \item{Ni}{the cumulative number of earthquakes of magnitude \eqn{\ge} \code{mi}}
#' @examples
#' mbin <- 0.1
#' beta <- log(10); mc <- 2
#' m <- mc - mbin / 2 + rexp(1e3, beta)
#' mdistr <- fmd(m)
#' View(mdistr)
#' @export
fmd <- function(m, mbin = 0.1) {
  mmin <- min(m, na.rm = T)
  mmax <- max(m, na.rm = T)
  mi <- seq(floor(mmin / mbin) * mbin, ceiling(mmax / mbin) * mbin, mbin)
  ni <- hist(m, breaks = c(mi, mmax + mbin) - mbin / 2, plot = F)$counts
  Ni <- sapply(1:length(ni), function(i) sum(ni[i:length(ni)]))
  return(data.frame(mi = mi, ni = ni, Ni = Ni))
}

#' Completeness Magnitude FMD-based Estimation
#'
#' Estimate the completeness magnitude \out{m<sub>c</sub>} from the
#' earthquake frequency magnitude distribution (FMD) using different published methods.
#'
#' \code{method = "mode"} calculates the mode of the vector of magnitudes
#' \code{m}. Applies to angular FMDs (Mignan, 2012), otherwise systematically
#' underestimates \out{m<sub>c</sub>}.
#'
#' \code{method = "mbass"} ("median-based analysis of the segment slope")
#' determines the main breakpoints of the earthquake FMD. \out{m<sub>c</sub>}
#' is defined as the change point that corresponds to the smallest probability of
#' making an error when rejecting the null-hypothesis in a Wilcoxon-Mann-Whitney test
#' (Amorese, 2007).
#'
#' \code{method = "gft"} estimates the goodness-of-fit between the cumulative
#' number of earthquakes observed and predicted by the Gutenberg-Richter model.
#' \out{m<sub>c</sub>} is defined as the lowest magnitude bin at which a fixed
#' threshold \emph{R} is first met. \emph{R} is defined as a normalized absolute
#' difference, fixed to 0.95. If the threshold is not reached, 0.90 is used. If again
#' the threshold is not reached, the \code{method = "mode"} is used instead
#' (Wiemer and Wyss, 2000).
#'
#' Both \code{"mode"} and \code{"mbass"} methods are non-parametric while \code{"gft"}
#' depends on the fitting of the Gutenberg-Richter model (see the function
#' \code{beta.mle}). For a general review of FMD-based \out{m<sub>c</sub>}
#' estimation methods, see Mignan and Woessner (2012). For further comparisons of
#' \code{"mbass"} and \code{"gft"}, see Mignan and Chouliaras (2014).
#'
#' @param m a numeric vector of earthquake magnitudes
#' @param method the method to be used: \code{"mode"}, \code{"mbass"}, or \code{"gft"}
#' (read Details)
#' @param mbin the magnitude binning value (if not provided, \code{mbin = 0.1})
#' @return The numeric value of \out{m<sub>c</sub>}.
#' @references Amorese, D. (2007), Applying a Change-Point Detecion Method on
#' Frequency-Magnitude Distributions, Bull. Seismol. Soc. Am., 97, 1742-1749,
#' \href{https://pubs.geoscienceworld.org/ssa/bssa/article-abstract/97/5/1742/146470/applying-a-change-point-detection-method-on}{doi: 10.1785/0120060181}
#' @references Mignan, A. (2012), Functional shape of the earthquake frequency-magnitude
#' distribution and completeness magnitude, J. Geophys. Res., 117, B08302,
#' \href{http://onlinelibrary.wiley.com/doi/10.1029/2012JB009347/full}{doi: 10.1029/2012JB009347}
#' @references Mignan, A., Woessner, J. (2012), Estimating the magnitude of completeness
#' for earthquake catalogs, Community Online Resource for Statistical Seismicity Analysis,
#' \href{http://www.corssa.org/export/sites/corssa/.galleries/articles-pdf/Mignan-Woessner-2012-CORSSA-Magnitude-of-completeness.pdf}{doi: 10.5078/corssa-00180805}
#' @references Mignan, A., Chouliaras, G. (2014), Fifty Years of Seismic Network
#' Performance in Greece (1964-2013): Spatiotemporal Evolution of the Completeness Magnitude,
#' Seismol. Res. Lett., 85, 657-667
#' \href{https://pubs.geoscienceworld.org/ssa/srl/article-abstract/85/3/657/315375/fifty-years-of-seismic-network-performance-in}{doi: 10.1785/0220130209}
#' @references Wiemer, S., Wyss, M. (2000), Minimum Magnitude of Completeness in
#' Earthquake Catalogs: Examples from Alaska, the Western United States, and Japan,
#' Bull. Seismol. Soc. Am.,
#' \href{https://pubs.geoscienceworld.org/ssa/bssa/article-abstract/90/4/859/120531/minimum-magnitude-of-completeness-in-earthquake}{90, 859-869}
#' @seealso \code{beta.mle}; \code{bfmd.sim}; \code{efmd.sim}; \code{fmd}
#' @examples
#' # Estimate mc for an angular FMD
#' theta <- list(kappa = 3 * log(10), beta = log(10), mc = 2)
#' m.angular <- efmd.sim(1e3, theta)
#' mc.val(m.angular, "mode")
#' mc.val(m.angular, "mbass")
#' mc.val(m.angular, "gft")
#'
#' # Estimate mc for a curved FMD
#' theta <- list(beta = log(10), mu = 2, sigma = 0.5)
#' m.curved <- bfmd.sim(1e4, theta)
#' mc.mode <- mc.val(m.curved, "mode")
#' mc.mbass <- mc.val(m.curved, "mbass")
#' mc.gft <- mc.val(m.curved, "gft")
#' mdistr <- fmd(m.curved)
#' plot(mdistr$mi, mdistr$Ni, log = "y")
#' points(mdistr$mi, mdistr$ni)
#' abline(v=c(mc.mode, mc.mbass, mc.gft), lty=c("dotted","solid","dashed"))
#'
#' # download the Southern California relocated catalogue of Hauksson et al. (2012)
#' url <- "http://service.scedc.caltech.edu/ftp/catalogs/"
#' dat <- "hauksson/Socal_DD/hs_1981_2016_comb_K4_A.cat_so_SCSN_v2q"
#' seism <- scan(paste(url, dat, sep = ""), what = "character", sep = "\n")
#' mbin <- 0.1
#' m <- round(as.numeric(substr(seism, start=63, stop=67)), digits = log10(1/mbin))
#' mc.mode <- mc.val(m, "mode")
#' mc.mbass <- mc.val(m, "mbass")
#' mc.gft <- mc.val(m, "gft")
#' mdistr <- fmd(m)
#' plot(mdistr$mi, mdistr$Ni, log = "y")
#' points(mdistr$mi, mdistr$ni)
#' abline(v=c(mc.mode, mc.mbass, mc.gft), lty=c("dotted","solid","dashed"))
#' @export
mc.val <- function(m, method, mbin = 0.1) {
  if (method == "mode") {
    dens <- density(m, from = min(m) - 1, to = max(m) + 1)
    return(round(dens$x[which(dens$y == max(dens$y))], digits = log10(1 / mbin))[1])
  }
  if (method == "mbass") {
    mdistr <- fmd(m, mbin)

    N <- length(mdistr$mi) - 1
    sl <- numeric(N)
    for(i in 1:N) sl[i] <- (log10(mdistr$ni[i+1]) -
      log10(mdistr$ni[i])) / (mdistr$mi[i+1] - mdistr$mi[i])
    mincr_corr <- mdistr$mi[2:length(mdistr$mi)]

    sl[which(is.infinite(sl) == T)] <- NA
    indnoNA <- which(is.na(sl) == F)
    N <- length(indnoNA)
    sl <- sl[indnoNA]; mincr_corr <- mincr_corr[indnoNA]

    niter <- 3
    j <- 0   #iterations
    k <- 0   #discontinuities
    SA <- numeric(N); pva <- numeric(); tau <- numeric()
    while(j < niter){
      for(i in 1:N) SA[i] <- abs(2 * sum(rank(sl)[1:i]) - i * (N + 1))
      indmax <- which.max(SA)[1]
      xn1 <- sl[1:indmax]
      xn2 <- sl[-(1:indmax)]
      if((indmax[1] > 2) && (indmax[1] <= (N-2)) &&
         (wilcox.test(xn1, xn2, exact = F, correct = T)[3] < 0.05)){
        k <- k + 1
        pva[k] <- wilcox.test(xn1, xn2, exact = F, correct = T)[3]
        tau[k] <- indmax[1]
        if(k > 1){
          medsl1 <- median(sl[1:n0])
          medsl2 <- median(sl[-(1:n0)])
          for(i in seq(1, n0, 1)) sl[i] <- sl[i] + medsl1
          for(i in seq(n0 + 1, length(sl), 1)) sl[i] <- sl[i] + medsl2
        }
        medsl1 <- median(sl[1:indmax[1]])
        medsl2 <- median(sl[-(1:indmax[1])])
        for(i in seq(1, indmax[1], 1)) sl[i] <- sl[i] - medsl1
        for(i in seq(indmax[1] + 1, length(sl), 1)) sl[i] <- sl[i] - medsl2
        n0 <- indmax[1]
      }
      j <- j + 1
    }
    Vpva <- as.vector(pva, mode = "numeric")
    ip <- order(Vpva)
    Break1 <- c(signif(mincr_corr[tau[ip[1]]]))
    Break2 <- c(signif(mincr_corr[tau[ip[2]]]))
    return(Break1)
  }
  if (method == "gft") {
    mdistr <- fmd(m, mbin)
    mmin <- mc.val(m, "mode", mbin)
    mscan <- seq(mmin, max(m) - mbin, mbin)
    nscan <- length(mscan)

    bi <- sapply(1:nscan, function(i) beta.mle(m[which(m > mscan[i] - mbin / 2)],
                                               mscan[i], mbin)/log(10))
    ai <- sapply(1:nscan, function(i) log10(length(which(m > mscan[i] - mbin / 2)))+
                   bi[i]*mscan[i])
    pred <- sapply(1:nscan, function(i) 10^(ai[i]-bi[i]*mdistr$mi))
    R <- sapply(1:nscan, function(i) sum(abs(mdistr$N[which(mdistr$mi > mscan[i]-mbin / 2)]-
                  unlist(pred[,i])[which(mdistr$mi > mscan[i]-mbin / 2)]))/
                  sum(mdistr$N[which(mdistr$mi > mscan[i]- mbin / 2)]))

    indGFT <- which(R <= 0.05)    #95% confidence
    if(length(indGFT) != 0) {
      mc <- mscan[indGFT[1]]
    } else {
      indGFT <- which(R <= 0.10)  #90% confidence
      if(length(indGFT) != 0) {
        mc <- mscan[indGFT[1]]
      } else {
        mc <- mmin
      }
    }
    return(mc)
  }
}

#' Gutenberg-Richter \eqn{\beta}-value
#'
#' Estimate the \eqn{\beta}-value (i.e. slope) of the Gutenberg-Richter model
#' (Gutenberg and Richter, 1944) by using the maximum likelihood estimation method (Aki, 1965).
#'
#' Note that \eqn{\beta} = \emph{b} log(10).
#'
#' @param m a numeric vector of earthquake magnitudes
#' @param mc the completeness magnitude value
#' @param mbin the magnitude binning value (if not provided, \code{mbin = 0.1})
#' @return The numeric value of \eqn{\beta}.
#' @references Aki, K. (1965), Maximum likelihood estimate of b in the formula log N =
#' a - bM and its confidence limits, Bull. Earthquake Res. Inst. Univ. Tokyo, 43, 237-239
#' @references Gutenberg, B., Richter, C.F. (1944), Frequency of earthquakes in California,
#' Bull. Seismol. Soc. Am.,
#' \href{https://authors.library.caltech.edu/47734/1/185.full.pdf}{34, 184-188}
#' @seealso \code{mc.val}
#' @examples
#' beta <- log(10); mc <- 2; mbin <- 0.1
#' m <- round(mc - mbin / 2 + rexp(1e3, beta), digits = log10(1/mbin))
#' beta.mle(m, mc, mbin)
#'
#' theta <- list(kappa = 3 * log(10), beta = 1.2*log(10), mc = 1.5)
#' m.angular <- efmd.sim(1e3, theta)
#' beta.mle(m.angular, theta$mc, mbin)
#' @export
beta.mle <- function(m, mc, mbin = 0.1) {
  1 / (mean(m[which(m > mc - mbin / 2)]) - (mc - mbin / 2))
}

#' \eqn{\chi}-value
#'
#' Estimate the \eqn{\chi}-value (i.e. slope) of the incomplete part of the
#' elemental frequency-magnitude distribution with \eqn{\chi} = \eqn{\kappa} - \eqn{\beta},
#' \eqn{\kappa} representing the earthquake detection parameter (Mignan, 2012).
#'
#' \eqn{\chi} is estimated similarly to \eqn{\beta}, by using the maximum likelihood
#' estimation method (Aki, 1965).
#'
#' @param m a numeric vector of earthquake magnitudes
#' @param mc the completeness magnitude value
#' @param mbin the magnitude binning value (if not provided, \code{mbin = 0.1})
#' @return The numeric value of \eqn{\chi}.
#' @references Aki, K. (1965), Maximum likelihood estimate of b in the formula log N =
#' a - bM and its confidence limits, Bull. Earthquake Res. Inst. Univ. Tokyo, 43, 237-239
#' @references Mignan, A. (2012), Functional shape of the earthquake frequency-magnitude
#' distribution and completeness magnitude, J. Geophys. Res., 117, B08302,
#' \href{http://onlinelibrary.wiley.com/doi/10.1029/2012JB009347/full}{doi: 10.1029/2012JB009347}
#' @seealso \code{beta.mle}; \code{efmd.sim}; \code{mc.val}
#' @examples
#' theta <- list(kappa = 2 * log(10), beta = log(10), mc = 2)
#' m.angular <- efmd.sim(1e4, theta)
#' mdistr <- fmd(m.angular)
#' plot(mdistr$mi, mdistr$Ni, log = "y")
#' points(mdistr$mi, mdistr$ni)
#' chi <- chi.mle(m.angular, theta$mc)
#' chi + theta$beta   # = kappa
#' beta <- beta.mle(m.sim, theta$mc)
#' abline(v = theta$mc, lty = "dotted", col = "red")
#' abline(a = log10(mdistr$ni[which(mdistr$mi >= theta$mc)[1]]) +
#'    beta / log(10) * theta$mc, b = -beta / log(10), col = "red")
#' abline(a = log10(mdistr$ni[which(mdistr$mi <= theta$mc)[length(which(mdistr$mi <= theta$mc))]]) -
#'    chi / log(10) * theta$mc, b = chi / log(10), col = "red")
#' @export
chi.mle <- function(m, mc, mbin = 0.1) {
  1 / ((mc - mbin / 2) - mean(m[which(m <= mc - mbin / 2)]))
}

#' Simulation of the Elemental FMD
#'
#' Simulate the elemental earthquake frequency magnitude distribution (eFMD) by
#' applying the Inversion Method (Devroye, 1986) to the angular FMD model of
#' Mignan (2012).
#'
#' The angular FMD model is defined as an Asymmetric Laplace distribution. It has an
#' angular shape in the log-lin space and corresponds to the case where the completeness
#' magnitude \out{m<sub>c</sub>} is constant (read more in Mignan, 2012; Mignan and
#' Chen, 2016).
#'
#' @param N the number of earthquakes to simulate
#' @param theta a list of 3 parameters:
#' * \code{kappa} the detection parameter value
#' * \code{beta}  the Gutenberg-Richter parameter value
#' * \code{mc}    the completeness magnitude value
#' @param mbin the magnitude binning value (if not provided, \code{mbin = 0.1})
#' @return A numeric vector of \code{N} earthquake magnitudes.
#' @references Devroye, L. (1986), Non-Uniform Random Variate Generation,
#' Springer-Verlag New York Inc., New York, 843 pp.
#' @references Mignan, A. (2012), Functional shape of the earthquake frequency-magnitude
#' distribution and completeness magnitude, J. Geophys. Res., 117, B08302,
#' \href{http://onlinelibrary.wiley.com/doi/10.1029/2012JB009347/full}{doi: 10.1029/2012JB009347}
#' @references Mignan, A., Chen, C.-C. (2016), The Spatial Scale of Detected Seismicity,
#' Pure Appl. Geophys., 173, 117-124,
#' \href{https://link.springer.com/article/10.1007/s00024-015-1133-7}{doi: 10.1007/s00024-015-1133-7}
#' @seealso \code{beta.mle}; \code{mc.val}
#' @examples
#' theta <- list(kappa = 3 * log(10), beta = log(10), mc = 2)
#' m.sim <- efmd.sim(1e4, theta)
#' mdistr <- fmd(m.sim)
#' plot(mdistr$mi, mdistr$Ni, log = "y")
#' points(mdistr$mi, mdistr$ni)
#' @export
efmd.sim <- function(N, theta, mbin = 0.1) {
  N.complete <- round(N * (theta$kappa - theta$beta) / theta$kappa)
  m.complete <- theta$mc - mbin / 2 + rexp(N.complete, theta$beta)
  m.incomplete <- theta$mc - mbin/2 -rexp(N - N.complete, (theta$kappa - theta$beta))
  return(round(c(m.incomplete, m.complete), digits = log10(1 / mbin)))
}

#' PDF of the Bulk FMD
#'
#' Compute the probability density function (PDF) of the bulk
#' frenquency-magnitude distribution (FMD), as defined by Ogata and Katsura (2006)
#' (see also Ringdal, 1975; Ogata and Katsura, 1993).
#'
#' The bulk FMD model is the product of the Gutenberg-Richter model and a detection function
#' defined as the cumulative Normal distribution. Its FMD has a curved shape in the log-lin
#' space and corresponds to the case where the completeness magnitude
#' \out{m<sub>c</sub>} is variable (the FMD curvature representing the
#' \out{m<sub>c</sub>} distribution; read more in Mignan, 2012; Mignan and Chen, 2016).
#'
#' @param m a numeric vector of earthquake magnitudes
#' @param theta a list of 3 parameters:
#' * \code{beta}  the Gutenberg-Richter parameter value
#' * \code{mu}    the mean value of the cumulative normal distribution
#' * \code{sigma} the standard deviation of the cumulative normal distribution
#' @return A numeric vector of densities.
#' @references Mignan, A. (2012), Functional shape of the earthquake frequency-magnitude
#' distribution and completeness magnitude, J. Geophys. Res., 117, B08302,
#' \href{http://onlinelibrary.wiley.com/doi/10.1029/2012JB009347/full}{doi: 10.1029/2012JB009347}
#' @references Mignan, A., Chen, C.-C. (2016), The Spatial Scale of Detected Seismicity,
#' Pure Appl. Geophys., 173, 117-124,
#' \href{https://link.springer.com/article/10.1007/s00024-015-1133-7}{doi: 10.1007/s00024-015-1133-7}
#' @references Ogata, Y., Katsura, K. (1993), Analysis of temporal and spatial
#' heterogeneity of magnitude frequency distribution inferred from earthquake
#' catalogues, Geophys. J. Int.,
#' \href{http://onlinelibrary.wiley.com/doi/10.1111/j.1365-246X.1993.tb04663.x/abstract}{113, 727-738}
#' @references Ogata, Y., Katsura, K. (2006), Immediate and updated forecasting of
#' aftershock hazard, Geophys. Res. Lett., 33, L10305,
#' \href{http://onlinelibrary.wiley.com/doi/10.1029/2006GL025888/full}{doi: 10.1029/2006GL025888}
#' @references Ringdal, F. (1975), On the estimation of seismic detection thresholds,
#' Bull. Seismol. Soc. Am.,
#' \href{https://pubs.geoscienceworld.org/ssa/bssa/article-abstract/65/6/1631/117495/on-the-estimation-of-seismic-detection-thresholds}{65, 1631-1642}
#' @seealso \code{beta.mle}; \code{bfmd.sim}
#' @export
bfmd.pdf <- function(m, theta) {
  theta$beta * exp(-theta$beta * (m - theta$mu) -
    0.5 * theta$beta^2 * theta$sigma^2) * pnorm(m, mean = theta$mu, sd = theta$sigma)
}

#' Simulation of the Bulk FMD
#'
#' Simulate the bulk earthquake frequency magnitude distribution (bFMD) by
#' applying the Thinning Method (Lewis and Shedler, 1979) to the curved FMD model of Ringdal (1975);
#' Ogata and Katsura (1993; 2006).
#'
#' The bulk FMD model is the product of the Gutenberg-Richter model and a detection function
#' defined as the cumulative Normal distribution. Its FMD has a curved shape in the log-lin
#' space and corresponds to the case where the completeness magnitude
#' \emph{\out{m<sub>c</sub>}} is variable (the FMD curvature representing the
#' \emph{\out{m<sub>c</sub>}} distribution; read more in Mignan, 2012; Mignan and Chen, 2016).
#'
#' @param N the approximate number of earthquakes to simulate
#' @param theta a list of 3 parameters:
#' * \code{beta}  the Gutenberg-Richter parameter value
#' * \code{mu}    the mean value of the cumulative normal distribution
#' * \code{sigma} the standard deviation of the cumulative normal distribution
#' @param mbin the magnitude binning value (if not provided, \code{mbin = 0.1})
#' @param mmin the minimum magnitude value (if not provided, \code{mmin = 0})
#' @param mmax the maximum magnitude value (if not provided, \code{mmax = 9})
#' @return A numeric vector of approximatively \code{N} earthquake magnitudes.
#' @references Lewis, P.A.W., Shedler, G.S. (1979), Simulation of nonhomogeneous poisson
#' processes by thinning, Naval Res. Logistics, 26, 403-413
#' \href{http://onlinelibrary.wiley.com/doi/10.1002/nav.3800260304/abstract}{doi: 10.1002/nav.3800260304}
#' @references Mignan, A. (2012), Functional shape of the earthquake frequency-magnitude
#' distribution and completeness magnitude, J. Geophys. Res., 117, B08302,
#' \href{http://onlinelibrary.wiley.com/doi/10.1029/2012JB009347/full}{doi: 10.1029/2012JB009347}
#' @references Mignan, A., Chen, C.-C. (2016), The Spatial Scale of Detected Seismicity,
#' Pure Appl. Geophys., 173, 117-124,
#' \href{https://link.springer.com/article/10.1007/s00024-015-1133-7}{doi: 10.1007/s00024-015-1133-7}
#' @references Ogata, Y., Katsura, K. (1993), Analysis of temporal and spatial
#' heterogeneity of magnitude frequency distribution inferred from earthquake
#' catalogues, Geophys. J. Int.,
#' \href{http://onlinelibrary.wiley.com/doi/10.1111/j.1365-246X.1993.tb04663.x/abstract}{113, 727-738}
#' @references Ogata, Y., Katsura, K. (2006), Immediate and updated forecasting of
#' aftershock hazard, Geophys. Res. Lett., 33, L10305,
#' \href{http://onlinelibrary.wiley.com/doi/10.1029/2006GL025888/full}{doi: 10.1029/2006GL025888}
#' @references Ringdal, F. (1975), On the estimation of seismic detection thresholds,
#' Bull. Seismol. Soc. Am.,
#' \href{https://pubs.geoscienceworld.org/ssa/bssa/article-abstract/65/6/1631/117495/on-the-estimation-of-seismic-detection-thresholds}{65, 1631-1642}
#' @seealso \code{beta.mle}; \code{bfmd.pdf}
#' @examples
#' theta <- list(beta = log(10), mu = 2, sigma = 0.3)
#' m.sim <- bfmd.sim(1e3, theta)
#' mdistr <- fmd(m.sim)
#' plot(mdistr$mi, mdistr$Ni, log = "y")
#' points(mdistr$mi, mdistr$ni)
#' @export
bfmd.sim <- function(N, theta, mbin = 0.1, mmin = 0, mmax = 9) {
  lambda.max <- max(bfmd.pdf(seq(mmin, mmax, mbin), theta))
  lambda.sum <- sum(bfmd.pdf(seq(mmin, mmax, mbin), theta))
  lambda.star <- lambda.sum
  N.star <- round(N * (mmax - mmin) / mbin)
  m.star <- round(runif(N.star, min = mmin, max = mmax), digits = log10(1 / mbin))
  m.curved <- sapply(1:length(m.star), function(i)
    if(runif(1) <= bfmd.pdf(m.star[i], theta)/lambda.star) m.star[i])
  return(unlist(m.curved))
}

