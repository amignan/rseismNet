#' rseismNet: Earthquake Frequency-Magnitude Distribution & Network Statistics
#'
#' A list of statistical functions to describe the earthquake frequency
#' magnitude distribution (FMD). The functions follow the two following rules:
#' (1) the complete FMD side (m>=mc) is governed by beta and depends on the fault
#' network properties; (2) the incomplete FMD part (m<mc) is governed by kappa and
#' depends on the seismic network properties (with m the earthquake magnitude and
#' mc the completeness magnitude). Both beta and kappa are computed using the Aki
#' approach. Any FMD shape is function of the mc spatial distribution, computed
#' using the Bayesian Magnitude of Completeness (BMC) method.
#'
#' @section rseismNet basic functions:
#' beta.mle, bfmd.pdf, bfmd.sim, efmd.sim, fmd, mc.val
#' @section rseismNet mapping functions:
#' bmc, bmc.bayes, bmc.prior, bmc.prior.default, d.geogr2km, mc.geogr
#' @references Aki, K. (1965), Maximum likelihood estimate of b in the formula log N =
#' a - bM and its confidence limits, Bull. Earthquake Res. Inst. Univ. Tokyo, 43, 237-239
#' @references Amorese, D. (2007), Applying a Change-Point Detecion Method on
#' Frequency-Magnitude Distributions, Bull. Seismol. Soc. Am., 97, 1742-1749,
#' \href{https://pubs.geoscienceworld.org/ssa/bssa/article-abstract/97/5/1742/146470/applying-a-change-point-detection-method-on}{doi: 10.1785/0120060181}
#' @references Devroye, L. (1986), Non-Uniform Random Variate Generation,
#' Springer-Verlag New York Inc., New York, 843 pp.
#' @references Gutenberg, B., Richter, C.F. (1944), Frequency of earthquakes in California,
#' Bull. Seismol. Soc. Am.,
#' \href{https://authors.library.caltech.edu/47734/1/185.full.pdf}{34, 184-188}
#' @references Kraft, T., Mignan, A., Giardini, D. (2013), Optimization of a large-scale
#' microseismic monitoring network in northern Switzerland, Geophys. J. Int., 195, 474-490,
#' \href{https://academic.oup.com/gji/article/195/1/474/601101}{doi: 10.1093/gji/ggt225}
#' @references Lewis, P.A.W., Shedler, G.S. (1979), Simulation of nonhomogeneous poisson
#' processes by thinning, Naval Res. Logistics, 26, 403-413
#' \href{http://onlinelibrary.wiley.com/doi/10.1002/nav.3800260304/abstract}{doi: 10.1002/nav.3800260304}
#' @references Mignan, A., Werner, M.J., Wiemer, S., Chen, C.-C., Wu, Y.-M. (2011),
#' Bayesian Estimation of the Spatially Varying Completeness Magnitude of Earthquake
#' Catalogs, Bull. Seismol. Soc. Am., 101, 1371-1385,
#' \href{https://pubs.geoscienceworld.org/bssa/article-lookup/101/3/1371}{doi: 10.1785/0120100223}
#' @references Mignan, A. (2012), Functional shape of the earthquake frequency-magnitude
#' distribution and completeness magnitude, J. Geophys. Res., 117, B08302,
#' \href{http://onlinelibrary.wiley.com/doi/10.1029/2012JB009347/full}{doi: 10.1029/2012JB009347}
#' @references Mignan, A., Woessner, J. (2012), Estimating the magnitude of completeness
#' for earthquake catalogs, Community Online Resource for Statistical Seismicity Analysis,
#' \href{http://www.corssa.org/export/sites/corssa/.galleries/articles-pdf/Mignan-Woessner-2012-CORSSA-Magnitude-of-completeness.pdf}{doi: 10.5078/corssa-00180805}
#' @references Mignan, A., Jiang, C., Zechar, J.D., Wiemer, S., Wu, Z., Huang, Z. (2013),
#' Completeness of the Mainland China Earthquake Catalog and Implications for the Setup of
#' the China Earthquake Forecast Texting Center, Bull. Seismol. Soc. Am., 103, 845-859,
#' \href{https://pubs.geoscienceworld.org/ssa/bssa/article/103/2A/845/331723/completeness-of-the-mainland-china-earthquake}{doi: 10.1785/0120120052}
#' @references Mignan, A., Chouliaras, G. (2014), Fifty Years of Seismic Network
#' Performance in Greece (1964-2013): Spatiotemporal Evolution of the Completeness Magnitude,
#' Seismol. Res. Lett., 85, 657-667
#' \href{https://pubs.geoscienceworld.org/ssa/srl/article-abstract/85/3/657/315375/fifty-years-of-seismic-network-performance-in}{doi: 10.1785/0220130209}
#' @references Mignan, A., Chen, C.-C. (2016), The Spatial Scale of Detected Seismicity,
#' Pure Appl. Geophys., 173, 117-124,
#' \href{https://link.springer.com/article/10.1007/s00024-015-1133-7}{doi: 10.1007/s00024-015-1133-7}
#' @references Panzera, F., Mignan, A., Vogfjord, K.S. (2017), Spatiotemporal evolution of
#' the completeness magnitude of the Icelandic earthquake catalogue from 1991 to 2013, J.
#' Seismol., 21, 615-630,
#' \href{https://link.springer.com/article/10.1007/s10950-016-9623-3}{doi: 10.1007/s10950-016-9623-3}
#' @references Tormann, T., Wiemer, S., Mignan, A. (2014), Systematic survey of
#' high-resolution b value imaging along Californian faults: inference on asperities, J.
#' Geophys. Res. Solid Earth, 119, 2029-2054,
#' \href{http://onlinelibrary.wiley.com/doi/10.1002/2013JB010867/full}{doi: 10.1002/2013JB010867}
#' @references Wiemer, S., Wyss, M. (2000), Minimum Magnitude of Completeness in
#' Earthquake Catalogs: Examples from Alaska, the Western United States, and Japan,
#' Bull. Seismol. Soc. Am.,
#' \href{https://pubs.geoscienceworld.org/ssa/bssa/article-abstract/90/4/859/120531/minimum-magnitude-of-completeness-in-earthquake}{90, 859-869}
#' @docType package
#' @name rseismNet

NULL
