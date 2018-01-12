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
#' @section rseismNet functions:
#' beta.mle, efmd.pdf, efmd.sim, fmd, mc.val
#' @references Mignan, A. (2012), Functional shape of the earthquake frequency-magnitude
#' distribution and completeness magnitude, J. Geophys. Res., 117, B08302,
#' \href{http://onlinelibrary.wiley.com/doi/10.1029/2012JB009347/full}{doi: 10.1029/2012JB009347}
#' @references Mignan, A., Woessner, J. (2012), Estimating the magnitude of completeness
#' for earthquake catalogs, Community Online Resource for Statistical Seismicity Analysis,
#' \href{http://www.corssa.org/export/sites/corssa/.galleries/articles-pdf/Mignan-Woessner-2012-CORSSA-Magnitude-of-completeness.pdf}{doi: 10.5078/corssa-00180805}
#' @docType package
#' @name rseismNet

NULL
