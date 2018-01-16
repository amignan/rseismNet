#' Completeness Magnitude Spatial Mapping
#'
#' Map the completeness magnitude \emph{\out{m<sub>c</sub>}} in space for a seismicity data
#' frame, following one of two possible geographical mappings: "grid" (computes
#' \emph{\out{m<sub>c</sub>}} in cells of size dbin) or "bmc" (computes mc according to the
#' Bayesian Magnitude of Completeness method of Mignan et al., Bull. Seismol. Soc. Am., 2011).
#'
#' Describe gridding...
#'
#' @param seism an earthquake catalog data frame of parameters:
#' * \code{lon} the earthquake longitude
#' * \code{lat} the earthquake latitude
#' * \code{m}   the earthquake magnitude
#' @param method the method to be used: \code{"mode"}, \code{"mbass"}, or \code{"gft"}
#' (see Details of function \code{mc.val})
#' @param mapping the mapping to be used: \code{"grid"} (see Details)
#' @param mbin the magnitude binning value (if not provided, \code{mbin = 0.1})
#' @param box a vector of the minimum longitude, maximum longitude, minimum latitude and
#' maximum latitude, in this order (if not provided, \code{box} is calculated from \code{seism})
#' @param dbin the spatial binning value (if not provided, \code{dbin} is calculated such that
#' the map is made of 10 longitudinal cells based on \code{box})
#' @return The completeness magnitude mc data frame:
#' * \code{lon} the longitude of the cell center
#' * \code{lat} the latitude of the cell center
#' * \code{mc}  the completeness magnitude \emph{\out{m<sub>c</sub>}} in the cell
#' @examples
#' url <- "http://service.scedc.caltech.edu/ftp/catalogs/"
#' cat <- "hauksson/Socal_DD/hs_1981_2011_06_comb_K2_A.cat_so_SCSN_v01"
#' dat <- scan(paste(url, cat, sep = ""), what = "character", sep = "\n")
#' lat <- as.numeric(substr(dat, start=35, stop=42))
#' lon <- as.numeric(substr(dat, start=44, stop=53))
#' m <- as.numeric(substr(dat, start=63, stop=67))
#' seism <- data.frame(lon = lon, lat = lat, m = m)
#' mc.obs <- mc.geogr(seism, "mode", "grid")
#' image(matrix(mc.obs$mc, nrow=length(unique(mc.obs$lon)), ncol=length(unique(mc.obs$lat))))
mc.geogr <- function(seism, method, mapping, mbin = 0.1, box = NULL, dbin = NULL, nmin = 4,
                     R = NULL, bmc.prior = NULL) {
  if(is.null(box)) box <- c(floor(min(seism$lon)), ceiling(max(seism$lon)),
                            floor(min(seism$lat)), ceiling(max(seism$lat)))
  if(is.null(dbin)) dbin <- (box[2]-box[1])/10

  grid <- expand.grid(lon = seq(box[1], box[2], dbin), lat = seq(box[3], box[4], dbin))
  grid.n <- nrow(grid)
  if (mapping == "grid") {
    ind.cell <- sapply(1:grid.n, function(i) which(seism$lon >= grid$lon[i] - dbin / 2 &
                                                   seism$lon < grid$lon[i] + dbin / 2 &
                                                   seism$lat >= grid$lat[i] - dbin / 2 &
                                                   seism$lat < grid$lat[i] + dbin / 2))
  }
  if (mapping == "circle.cst") {
    ind.circle <- sapply(1:grid.n, function(i) NULL)
  }
  if (mapping == "circle.opt") {
    ind.circle <- sapply(1:grid.n, function(i) NULL)
  }

  mc.cell <- sapply(1:grid.n, function(i) if(length(unlist(ind.cell[i])) >= nmin)
  {mc.val(seism$m[unlist(ind.cell[i])], method, mbin)} else {NA})
  return(data.frame(grid, mc=unlist(mc.cell)))
}

bmc.prior.generic <- function(kth) {
  params <- NULL
  if(kth == 3) params <- list(c1 = 4.81, c2 = 0.0883, c3 = -4.36,
                              sigma = 0.19, kth = 3, support = "generic")
  if(kth == 4) params <- list(c1 = 5.96, c2 = 0.0803, c3 = -5.80,
                              sigma = 0.18, kth = 4, support = "generic")
  if(kth == 5) params <- list(c1 = 9.42, c2 = 0.0598, c3 = -9.60,
                              sigma = 0.18, kth = 5, support = "generic")
  return(params)
}

#' BMC Prior
#'
#' Define the prior model of the Bayesian Magnitude of Completeness (BMC) method (Mignan et
#' al., 2011) by fitting a function of the form xx with mc the observed completeness magnitude
#' in any given cell and d the distance to the nearest kth seismic station
#'
#' @param mc.obs a data frame of the observed mc map, defined by the function \code{mc.geogr}
#' @param stations the seismic station data frame of parameters:
#' * \code{lon} the station longitude
#' * \code{lat} the station latitude
#' @param kth the kth seismic station used for distance calculation (if not provided,
#' \code{kth = 4})
#' @param support the information supporting the prior model: \code{"calibrated"} or
#' \code{"data"} (read Details)
#' @seealso \code{bmc}; \code{bmc.prior.generic}; \code{mc.geogr}
#' @examples
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
#' url <- "http://service.scedc.caltech.edu/station/weblist.php"
#' dat <- scan(url, what = "character", sep = "\n", skip = 7)
#' network <- substr(dat, start = 1, stop = 2)
#' sta.name <- substr(dat, start = 5, stop = 9)
#' sta.lat <- as.numeric(substr(dat, start = 52, stop = 59))
#' sta.lon <- as.numeric(substr(dat, start = 61, stop = 70))
#' sta.on <- as.numeric(substr(dat, start = 78, stop = 81))
#' sta.off <- as.numeric(substr(dat, start = 89, stop = 92))
#' stations <- data.frame(lon = sta.lon, lat = sta.lat)
#' stations <- subset(stations, (network == "CI" & sta.off > min(yr) & sta.on < max(yr)))
#' stations <- subset(stations, (duplicated(sta.name) == F))
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
    params <- bmc.prior.generic(kth)
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
    if(is.null(params)) params <- bmc.prior(mc.obs, stations, kth = 4)
  }

  return(list(params, data.frame(mc.obs, d.kth = d.kth)))
}


