#' Completeness Magnitude Mapping
#'
#' Map the completeness magnitude \emph{\out{m<sub>c</sub>}} for a seismicity data
#' frame \code{seism}, following one of three possible geographical mappings (see Details).
#'
#' \code{"mapping = grid"} computes \emph{\out{m<sub>c</sub>}} for earthquakes located in
#' each cell of bin size \code{dbin}.
#'
#' \code{"mapping = circle.cst"} computes \emph{\out{m<sub>c</sub>}} for earthquakes located
#' in cylinders centered on each cell and with radius \code{R} (see the smoothing impact of
#' increasing R in Mignan et al., 2011).
#'
#' \code{"mapping = circle.opt"} computes \emph{\out{m<sub>c</sub>}} for earthquakes located
#' in cylinders centered on each cell and with variable radius estimated from the seismic
#' network spatial density (Mignan et al., 2011). This method aims at minimizing
#' \emph{\out{m<sub>c</sub>}} spatial heterogeneities by using a small radius in the dense
#' parts of the network where \emph{\out{m<sub>c</sub>}} changes faster. The method depends
#' on the default BMC prior model parameters defined in function \code{bmc.prior.default}.
#'
#' @param seism an earthquake catalog data frame of parameters:
#' * \code{lon} the earthquake longitude
#' * \code{lat} the earthquake latitude
#' * \code{m}   the earthquake magnitude
#' * \code{...} other earthquake parameters
#' @param method the method to be used to evaluate \emph{\out{m<sub>c</sub>}}: \code{"mode"}, \code{"mbass"}, or \code{"gft"}
#' (see Details of function \code{mc.val})
#' @param mapping the mapping to be used: \code{"grid"}, \code{"circle.cst"}, or
#' \code{"circle.opt"} (see Details)
#' @param mbin the magnitude binning value (if not provided, \code{mbin = 0.1})
#' @param box a vector of the minimum longitude, maximum longitude, minimum latitude and
#' maximum latitude, in this order (if not provided, \code{box} is calculated from
#' the geographical limits of \code{seism})
#' @param dbin the spatial binning value (if not provided, \code{dbin} is calculated such that
#' the map is made of 10 longitudinal cells based on \code{box})
#' @param nmin the minimum number of earthquakes required to evaluate
#' \emph{\out{m<sub>c</sub>}} (if not provided, \code{nmin = 4} for \code{method = "mode"},
#' otherwise \code{nmin = 50}; see explanation in Mignan et al., 2011)
#' @param R the cylinder radius in kilometers for \code{"mapping = circle.cst"} (if not provided,
#' \code{R = 30})
#' @param stations the seismic network data frame for \code{"mapping = circle.opt"},
#' of parameters:
#' * \code{lon} the seismic station longitude
#' * \code{lat} the seismic station latitude
#' * \code{...} other station attributes
#' @param kth the \out{<i>k</i><sup>th</sup>} seismic station used for distance calculation
#' for \code{"mapping = circle.opt"} (if not provided, \code{kth = 4})
#' @return The data frame of 3 parameters:
#' * \code{lon} the longitude of the cell center
#' * \code{lat} the latitude of the cell center
#' * \code{mc.obs}  the completeness magnitude \emph{\out{m<sub>c</sub>}} in the cell
#' @references Mignan, A., Werner, M.J., Wiemer, S., Chen, C.-C., Wu, Y.-M. (2011),
#' Bayesian Estimation of the Spatially Varying Completeness Magnitude of Earthquake
#' Catalogs, Bull. Seismol. Soc. Am., 101, 1371-1385,
#' \href{https://pubs.geoscienceworld.org/bssa/article-lookup/101/3/1371}{doi: 10.1785/0120100223}
#' @seealso \code{bmc.prior.default}; \code{mc.val}
#' @examples
#' # download the Southern California relocated catalogue of Hauksson et al. (2012)
#' url <- "http://service.scedc.caltech.edu/ftp/catalogs/"
#' cat <- "hauksson/Socal_DD/hs_1981_2011_06_comb_K2_A.cat_so_SCSN_v01"
#' dat <- scan(paste(url, cat, sep = ""), what = "character", sep = "\n")
#' yr <- as.numeric(substr(dat, start=1, stop=4))
#' lat <- as.numeric(substr(dat, start=35, stop=42))
#' lon <- as.numeric(substr(dat, start=44, stop=53))
#' m <- as.numeric(substr(dat, start=63, stop=67))
#' seism <- data.frame(yr = yr, lon = lon, lat = lat, m = m)
#'
#' # reduce catalogue size for faster computation
#' seism <- subset(seism, yr >= 2008)
#'
#' # map mc in a grid
#' mc.grid <- mc.geogr(seism, "mode", "grid", dbin = 0.1)
#' image(matrix(mc.grid$mc.obs, nrow=length(unique(mc.grid$lon)), ncol=length(unique(mc.grid$lat))))
#'
#' # map mc with cylinder smoothing (this takes a few minutes!)
#' mc.cst <- mc.geogr(seism, "mbass", "circle.cst", dbin = 0.1)
#' image(matrix(mc.cst$mc.obs, nrow=length(unique(mc.cst$lon)), ncol=length(unique(mc.cst$lat))))
#'
#' # download the Southern California seismic network data
#' url <- "http://service.scedc.caltech.edu/station/weblist.php"
#' dat <- scan(url, what = "character", sep = "\n", skip = 7)
#' network <- substr(dat, start = 1, stop = 2)
#' sta.name <- substr(dat, start = 5, stop = 9)
#' sta.lat <- as.numeric(substr(dat, start = 52, stop = 59))
#' sta.lon <- as.numeric(substr(dat, start = 61, stop = 70))
#' sta.on <- as.numeric(substr(dat, start = 78, stop = 81))
#' sta.off <- as.numeric(substr(dat, start = 89, stop = 92))
#' stations <- data.frame(lon = sta.lon, lat = sta.lat, name = sta.name)
#' stations <- subset(stations, (network == "CI" & sta.off > min(seism$yr) & sta.on < max(seism$yr)))
#' stations <- subset(stations, (duplicated(name) == F))
#'
#' # map mc with optimal sample size (this takes a few minutes!)
#' # (minimizes mc heterogeneities while maximizing sample size)
#' mc.opt <- mc.geogr(seism, "mode", "circle.opt", dbin = 0.1, stations = stations)
#' image(matrix(mc.opt$mc.obs, nrow=length(unique(mc.opt$lon)), ncol=length(unique(mc.opt$lat))))
mc.geogr <- function(seism, method, mapping, mbin = 0.1, box = NULL, dbin = NULL, nmin = NULL,
                     R = 30, stations = NULL, kth = 4, n.sample = 0) {
  if(is.null(box)) box <- c(floor(min(seism$lon)), ceiling(max(seism$lon)),
                            floor(min(seism$lat)), ceiling(max(seism$lat)))
  if(is.null(dbin)) dbin <- (box[2]-box[1])/10
  if(is.null(nmin)) {
    if(method == "mode") nmin <- 4 else nmin <- 50
  }

  grid <- expand.grid(lon = seq(box[1], box[2], dbin), lat = seq(box[3], box[4], dbin))
  grid.n <- nrow(grid)
  if (mapping == "grid") {
    ind.cell <- sapply(1:grid.n, function(i) which(seism$lon >= grid$lon[i] - dbin / 2 &
                                                   seism$lon < grid$lon[i] + dbin / 2 &
                                                   seism$lat >= grid$lat[i] - dbin / 2 &
                                                   seism$lat < grid$lat[i] + dbin / 2))
  }
  if (mapping == "circle.cst") {
    pt.grid <- matrix(c(grid$lon, grid$lat), nrow = grid.n, ncol = 2)
    pt.seism <- matrix(c(seism$lon, seism$lat), nrow = nrow(seism), ncol = 2)
    r_m <- sapply(1:grid.n, function(i) geosphere::distHaversine(pt.grid[i, ], pt.seism))
    ind.cell <- sapply(1:grid.n, function(i) which(r_m[,i] * 1e-3 <= R))
  }
  if (mapping == "circle.opt") {
    sta.n <- nrow(stations)
    pt.grid <- matrix(c(grid$lon, grid$lat), nrow = grid.n, ncol = 2)
    pt.sta <- matrix(c(stations$lon, stations$lat), nrow = sta.n, ncol = 2)
    pt.seism <- matrix(c(seism$lon, seism$lat), nrow = nrow(seism), ncol = 2)
    r_m <- sapply(1:grid.n, function(i) geosphere::distHaversine(pt.grid[i, ], pt.seism))
    d_m <- sapply(1:grid.n, function(i) geosphere::distHaversine(pt.grid[i, ], pt.sta))
    d.kth <- sapply(1:grid.n, function(i) sort(d_m[,i])[kth] * 1e-3)  # in km
    params <- bmc.prior.default(kth)
    mc.pred <- (params$c1 * d.kth ^ params$c2 + params$c3)
    dlow <- (((mc.pred - params$sigma) - params$c3) / params$c1) ^ (1 / params$c2)
    dhigh <- (((mc.pred + params$sigma) - params$c3) / params$c1) ^ (1 / params$c2)
    R <- (dhigh-dlow)/2
    Rmin <- sqrt(2) * dbin / 2 * 111
    R[which(R <= Rmin)] <- Rmin
    ind.cell <- sapply(1:grid.n, function(i) which(r_m[,i] * 1e-3 <= R))
  }

  if(n.sample > 0) {
    mc.unc <- sapply(1:grid.n, function(i) if(length(unlist(ind.cell[i])) >= nmin) {
      sapply(1:n.sample, function(j) mc.val(
          sample(seism$m[unlist(ind.cell[i])], replace = T), method, mbin
        ))
    } else {NA})
    mc.cell <- sapply(1:grid.n, function(i) round(mean(mc.unc[[i]], na.rm = T),
                                                  digits = log10(1 / mbin)))
    sigma.cell <- sapply(1:grid.n, function(i) sd(mc.unc[[i]], na.rm = T))
  } else {
    mc.cell <- sapply(1:grid.n, function(i) if(length(unlist(ind.cell[i])) >= nmin)
    {mc.val(seism$m[unlist(ind.cell[i])], method, mbin)} else {NA})
    sigma.cell <- rep(NA, grid.n)
  }
  return(data.frame(grid, mc.obs = unlist(mc.cell), sigma.obs = sigma.cell))
}

#' Default BMC Prior
#'
#' List the parameters of the default prior model of the Bayesian Magnitude of
#' Completeness (BMC) method, as defined in Mignan et al. (2011) for different
#' \out{<i>k</i><sup>th</sup>} values (only \code{kth = 3}, \code{4} and \code{5}
#' allowed, otherwise returns \code{NULL}).
#'
#' The generic model is defined as the prior model evaluated for the Taiwan earthquake data
#' (Mignan et al., 2011), as it represents the best constrained data set so far (read more
#' on this in Mignan and Chouliaras, 2014). It often provides a better prior once calibrated
#' to other data than a new fit to the data (see example given for the function
#' \code{bmc.prior}). It is also used for rapid mapping of the optimal
#' \emph{\out{m<sub>c</sub>}} map (see function \code{mc.geogr}).
#'
#' @param kth the \out{<i>k</i><sup>th</sup>} seismic station used for distance calculation (if not provided,
#' \code{kth = 4})
#' @return the BMC prior parameter list
#' * \code{c1}, \code{c2}, \code{c3} the empirical parameters
#' * \code{sigma} the standard error
#' * \code{kth}   the kth seismic station used for distance calculation
#' * \code{support}   the information supporting the prior model (here
#' \code{support = "default"})
#' @references Mignan, A., Werner, M.J., Wiemer, S., Chen, C.-C., Wu, Y.-M. (2011),
#' Bayesian Estimation of the Spatially Varying Completeness Magnitude of Earthquake
#' Catalogs, Bull. Seismol. Soc. Am., 101, 1371-1385,
#' \href{https://pubs.geoscienceworld.org/bssa/article-lookup/101/3/1371}{doi: 10.1785/0120100223}
#' @references Mignan, A., Chouliaras, G. (2014), Fifty Years of Seismic Network
#' Performance in Greece (1964-2013): Spatiotemporal Evolution of the Completeness Magnitude,
#' Seismol. Res. Lett., 85, 657-667
#' \href{https://pubs.geoscienceworld.org/ssa/srl/article-abstract/85/3/657/315375/fifty-years-of-seismic-network-performance-in}{doi: 10.1785/0220130209}
#' @seealso \code{bmc}; \code{bmc.prior}; \code{mc.geogr}
#' @examples
#' # map the predicted mc for a set of simulated stations
#' box <- c(-5, 5, -5, 5); dbin <- 0.1  #degrees
#' sta.n <- 30
#' stations <- data.frame(lon = rnorm(sta.n), lat=rnorm(sta.n))
#' grid <- expand.grid(lon = seq(box[1], box[2], dbin), lat = seq(box[3], box[4], dbin))
#' grid.n <- nrow(grid)
#'
#' kth <- 4
#' params <- bmc.prior.default(kth)
#' pt.grid <- matrix(c(grid$lon, grid$lat), nrow = grid.n, ncol = 2)
#' pt.sta <- matrix(c(stations$lon, stations$lat), nrow = sta.n, ncol = 2)
#' d_m <- sapply(1:grid.n, function(i) geosphere::distHaversine(pt.grid[i, ], pt.sta))
#' d.kth <- sapply(1:grid.n, function(i) sort(d_m[,i])[kth] * 1e-3)  # in km
#' mc.pred <- (params$c1 * d.kth ^ params$c2 + params$c3)
#' image(unique(grid$lon), unique(grid$lat),
#'   matrix(mc.pred, nrow=length(unique(grid$lon)), ncol=length(unique(grid$lat))))
#' points(stations, pch = 2)
bmc.prior.default <- function(kth) {
  params <- NULL
  if(kth == 3) params <- list(c1 = 4.81, c2 = 0.0883, c3 = -4.36,
                              sigma = 0.19, kth = 3, support = "default")
  if(kth == 4) params <- list(c1 = 5.96, c2 = 0.0803, c3 = -5.80,
                              sigma = 0.18, kth = 4, support = "default")
  if(kth == 5) params <- list(c1 = 9.42, c2 = 0.0598, c3 = -9.60,
                              sigma = 0.18, kth = 5, support = "default")
  return(params)
}

#' BMC Prior
#'
#' Define the prior model of the Bayesian Magnitude of Completeness (BMC) method (Mignan et
#' al., 2011) by fitting a function of the form xx with \emph{\out{m<sub>c</sub>}} the observed completeness magnitude
#' in any given cell and d the distance to the kth nearest seismic station
#'
#' \code{support = "calibrated"} uses the generic BMC prior defined by function
#' \code{bmc.prior.generic} and calibrates it to the \code{mc.obs} data by shifting the
#' residual average to 0, substracting it from c3 (see e.g., Mignan et al., 2013; Mignan
#' and Chouliaras, 2014).
#'
#' \code{support = "data"} directly fits the prior function to the \code{mc.obs} data
#' by using the Nonlinear Least Squares \code{stats::nls} function. If the estimation fails,
#' \code{support = "calibrated"} is automatically used instead.
#'
#' @param mc.obs a data frame of the observed mc map, defined by the function \code{mc.geogr}
#' @param stations the seismic network data frame of parameters:
#' * \code{lon} the seismic station longitude
#' * \code{lat} the seismic station latitude
#' * \code{...} other station attributes
#' @param kth the \code{kth = 3} seismic station used for distance calculation (if not provided,
#' \code{kth = 4})
#' @param support the information supporting the prior model: \code{"calibrated"} or
#' \code{"data"} (if not provided, \code{support = "calibrated"} - read Details)
#' @return A list of
#' @return the BMC prior parameter list
#' * \code{c1}, \code{c2}, \code{c3} the empirical parameters
#' * \code{sigma} the standard error
#' * \code{kth}   the kth seismic station used for distance calculation
#' * \code{support}   the information supporting the prior model
#' @return the input data frame
#' * \code{mc}   the completeness magnitude value per cell
#' * \code{d.kth}   the distance to the kth nearest seismic station per cell
#' @references Mignan, A., Werner, M.J., Wiemer, S., Chen, C.-C., Wu, Y.-M. (2011),
#' Bayesian Estimation of the Spatially Varying Completeness Magnitude of Earthquake
#' Catalogs, Bull. Seismol. Soc. Am., 101, 1371-1385,
#' \href{https://pubs.geoscienceworld.org/bssa/article-lookup/101/3/1371}{doi: 10.1785/0120100223}
#' @references Mignan, A., Jiang, C., Zechar, J.D., Wiemer, S., Wu, Z., Huang, Z. (2013),
#' Completeness of the Mainland China Earthquake Catalog and Implications for the Setup of
#' the China Earthquake Forecast Texting Center, Bull. Seismol. Soc. Am., 103, 845-859,
#' \href{https://pubs.geoscienceworld.org/ssa/bssa/article/103/2A/845/331723/completeness-of-the-mainland-china-earthquake}{doi: 10.1785/0120120052}
#' @references Mignan, A., Chouliaras, G. (2014), Fifty Years of Seismic Network
#' Performance in Greece (1964-2013): Spatiotemporal Evolution of the Completeness Magnitude,
#' Seismol. Res. Lett., 85, 657-667
#' \href{https://pubs.geoscienceworld.org/ssa/srl/article-abstract/85/3/657/315375/fifty-years-of-seismic-network-performance-in}{doi: 10.1785/0220130209}
#' @seealso \code{bmc}; \code{bmc.prior.default}; \code{mc.geogr}
#' @examples
#' # download the Southern California relocated catalogue of Hauksson et al. (2012)
#' url <- "http://service.scedc.caltech.edu/ftp/catalogs/"
#' cat <- "hauksson/Socal_DD/hs_1981_2011_06_comb_K2_A.cat_so_SCSN_v01"
#' dat <- scan(paste(url, cat, sep = ""), what = "character", sep = "\n")
#' yr <- as.numeric(substr(dat, start=1, stop=4))
#' lat <- as.numeric(substr(dat, start=35, stop=42))
#' lon <- as.numeric(substr(dat, start=44, stop=53))
#' m <- as.numeric(substr(dat, start=63, stop=67))
#' seism <- data.frame(yr = yr, lon = lon,lat = lat, m = m)
#' mc.obs <- mc.geogr(seism, "mode", "grid", dbin = 0.1)
#'
#' # download the Southern California seismic network data
#' url <- "http://service.scedc.caltech.edu/station/weblist.php"
#' dat <- scan(url, what = "character", sep = "\n", skip = 7)
#' network <- substr(dat, start = 1, stop = 2)
#' sta.name <- substr(dat, start = 5, stop = 9)
#' sta.lat <- as.numeric(substr(dat, start = 52, stop = 59))
#' sta.lon <- as.numeric(substr(dat, start = 61, stop = 70))
#' sta.on <- as.numeric(substr(dat, start = 78, stop = 81))
#' sta.off <- as.numeric(substr(dat, start = 89, stop = 92))
#' stations <- data.frame(lon = sta.lon, lat = sta.lat, name = sta.name)
#' stations <- subset(stations, (network == "CI" & sta.off > min(seism$yr) & sta.on < max(seism$yr)))
#' stations <- subset(stations, (duplicated(name) == F))
#'
#' # test the two possible priors & plot
#' model.calibrated <- bmc.prior(mc.obs, stations, kth = 5)
#' model.fromdata <- bmc.prior(mc.obs, stations, kth = 5, support = "data")
#' data <- model.calibrated[[2]]
#' di <- seq(0, max(data$d.kth))
#' params.cal <- model.calibrated[[1]]
#' params.dat <- model.fromdata[[1]]
#' plot(data$d.kth, data$mc)
#' lines(di, params.cal$c1*di^params.cal$c2+params.cal$c3, col="orange")
#' lines(di, params.dat$c1*di^params.dat$c2+params.dat$c3, col="red")
bmc.prior <- function(mc.obs, stations, kth = 4, support = "calibrated") {
  grid.n <- nrow(mc.obs)
  sta.n <- nrow(stations)
  pt.grid <- matrix(c(mc.obs$lon, mc.obs$lat), nrow=grid.n, ncol=2)
  pt.sta <- matrix(c(stations$lon, stations$lat), nrow=sta.n, ncol=2)

  d_m <- sapply(1:grid.n, function(i) geosphere::distHaversine(pt.grid[i, ], pt.sta))
  d.kth <- sapply(1:grid.n, function(i) sort(d_m[,i])[kth] * 1e-3)  # in km
  dat2fit <- data.frame(d = d.kth, mc = mc.obs$mc)

  if(support == "calibrated") {
    params <- bmc.prior.default(kth)
    residual <- dat2fit$mc - (params$c1 * dat2fit$d ^ params$c2 + params$c3)
    params$c3 <- params$c3 + mean(residual, na.rm = T)
    residual.calib <- dat2fit$mc - params$c1 * dat2fit$d ^ params$c2 + params$c3
    params$sigma <- sd(residual.calib, na.rm = T)
  }
  if(support == "data") {
    params <- NULL
    tryCatch({
      fit <- nls(mc ~ c1 * d ^ c2 + c3, data = dat2fit, start = list(c1 = 1, c2 = 0.5, c3 = -1))
      c1 = coef(fit)["c1"]; c2 = coef(fit)["c2"]; c3 = coef(fit)["c3"]
      sigma <- sd(resid(fit), na.rm=T)
      params <- list(c1 = c1, c2 = c2, c3 = c3, sigma = sigma, kth = kth, support = "data")
    }, error = function(e) cat("ERROR :",conditionMessage(e),
                               "-> Calibration used instead", "\n"))
    if(is.null(params)) params <- bmc.prior(mc.obs, stations, kth = 4, support = "calibrated")
  }

  return(list(params, data.frame(mc = mc.obs$mc, d.kth = d.kth)))
}

#' BMC Method
#'
#' Runs all the steps of the Bayesian Magnitude of Completeness (BMC) method (Mignan et al.,
#' 2011) and produces a geographical map data frame of the completeness magnitude (observed,
#' predicted, and posterior) and associated uncertainties (observed, predicted, and
#' posterior).
#'
#' It is a wrap-up of other functions. See Examples for a description of the different steps.
#' The BMC method has been used in
#'
bmc <- function(seism, stations, support = "fast", mbin = 0.1, box = NULL, dbin = NULL,
                nmin = NULL, kth = 4) {
  if(support == "fast") {
    cat("Compute observed mc map (optimized)")
    mc.obs <- mc.geogr(seism, "mode", "circle.opt", mbin = mbin, box = box, dbin = dbin,
                       nmin = 4, R = NULL, stations = stations, kth = kth)

    cat("Compute predicted mc map (calibrated default prior model)")
    params <- bmc.prior(mc.obs, stations, kth = kth, support = "calibrated")
    grid.n <- nrow(mc.obs)
    pt.grid <- matrix(c(mc.obs$lon, mc.obs$lat), nrow = grid.n, ncol = 2)
    pt.sta <- matrix(c(stations$lon, stations$lat), nrow = sta.n, ncol = 2)
    d_m <- sapply(1:grid.n, function(i) geosphere::distHaversine(pt.grid[i, ], pt.sta))
    d.kth <- sapply(1:grid.n, function(i) sort(d_m[,i])[kth] * 1e-3)  # in km
    mc.pred <- (params$c1 * d.kth ^ params$c2 + params$c3)

    cat("Compute posterior mc map (combining both observed & predicted mc maps)")

  }
}




