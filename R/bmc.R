#' Distance between Geographical Points
#'
#' Calculate the distance d (in kilometers) between two geographical points or between one
#' geographical point and a vector of geographical points.
#'
#' The \code{method = "fast"} assumes that the surface is flat. The
#' \code{method = "haversine"}, slower, calculates the haversine distance by using the
#' function \code{geosphere::distHaversine}. For regional earthquake catalogues, the error
#' of the fast method is insignificant.
#'
#' @param pt.ref a list of 2 parameters:
#' * \code{lon} the reference location longitude
#' * \code{lat} the reference location latitude
#' @param pt.list a list or data frame of 2 parameters:
#' * \code{lon} the other location longitude(s)
#' * \code{lat} the other location latitude(s)
#' @param method the method to be used to evaluate d: \code{"fast"} or \code{"haversine"}
#' (see Details)
#' @return the value(s) of the distance d between \code{pt.ref} and \code{pt.list}.
d.geogr2km <- function(pt.ref, pt.list, method) {
  if(method == "fast") {
    rad_earth <- 6378.1   #km
    lat_km <- rad_earth*pi/180
    lon_km <- rad_earth*cos(pt.ref$lat*pi/180)*pi/180
    d <- sqrt(((pt.ref$lon - pt.list$lon) * lon_km) ^ 2 +
              ((pt.ref$lat - pt.list$lat) * lat_km) ^ 2)
  }
  if(method == "haversine") {
    pt1.mat <- matrix(c(pt.ref$lon, pt.ref$lat), nrow = 1, ncol = 2)
    pt2.mat <- matrix(c(pt.list$lon, pt.list$lat), nrow = nrow(pt.list), ncol = 2)
    d_m <- geosphere::distHaversine(pt1.mat, pt2.mat)
    d <- d_m*1e-3
  }
  return(d)
}


#' Completeness Magnitude Mapping
#'
#' Map the completeness magnitude \out{m<sub>c</sub>} for a seismicity data
#' frame \code{seism}, following one of three possible geographical mappings (see Details).
#'
#' \code{"mapping = grid"} computes \out{m<sub>c</sub>} for earthquakes located in
#' each cell of bin size \code{dbin}.
#'
#' \code{"mapping = circle.cst"} computes \out{m<sub>c</sub>} for earthquakes located
#' in cylinders centered on each cell and with radius \code{R} (see the smoothing impact of
#' increasing radius in Mignan et al., 2011).
#'
#' \code{"mapping = circle.opt"} computes \out{m<sub>c</sub>} for earthquakes located
#' in cylinders centered on each cell and with variable radius estimated from the seismic
#' network spatial density (Mignan et al., 2011). This method aims at minimizing
#' \out{m<sub>c</sub>} spatial heterogeneities by using a smaller radius in the dense
#' parts of the network where \out{m<sub>c</sub>} changes faster.
#'
#' The completeness magnitude \out{m<sub>c</sub>} is evaluated for \code{n.bootstrap}
#' bootstraps, with \code{mc.obs} the mean value and \code{sigma.obs} the standard error. For
#' \code{n.bootstrap = 0}, \code{sigma.obs = NA}.
#'
#' The function calls \code{geosphere::distHaversine} to calculate distances between
#' geographical points.
#'
#' @param seism an earthquake catalog data frame of parameters:
#' * \code{lon} the earthquake longitude
#' * \code{lat} the earthquake latitude
#' * \code{m}   the earthquake magnitude
#' * \code{...} other earthquake parameters
#' @param method the method to be used to evaluate \out{m<sub>c</sub>}: \code{"mode"},
#' \code{"mbass"}, or \code{"gft"} (see Details of function \code{mc.val})
#' @param mapping the mapping to be used: \code{"grid"}, \code{"circle.cst"}, or
#' \code{"circle.opt"} (see Details)
#' @param mbin the magnitude binning value (if not provided, \code{mbin = 0.1})
#' @param box a vector of the minimum longitude, maximum longitude, minimum latitude and
#' maximum latitude, in this order (if not provided, \code{box} is calculated from
#' the geographical limits of \code{seism})
#' @param dbin the spatial binning value (if not provided, \code{dbin} is calculated such that
#' the map is made of 10 longitudinal cells based on \code{box})
#' @param nmin the minimum number of earthquakes required to evaluate
#' \out{m<sub>c</sub>} (if not provided, \code{nmin = 4} for \code{method = "mode"},
#' otherwise \code{nmin = 50}; see explanation in Mignan et al., 2011)
#' @param R the cylinder radius in kilometers for \code{"mapping = circle.cst"} (if not provided,
#' \code{R = 30})
#' @param stations the seismic network data frame for \code{"mapping = circle.opt"},
#' of parameters:
#' * \code{lon} the seismic station longitude
#' * \code{lat} the seismic station latitude
#' * \code{...} other station attributes
#' @param kth the \out{k<sup>th</sup>} nearest seismic station used for
#' distance calculation in \code{"mapping = circle.opt"} (if not provided, \code{kth = 4})
#' @param params the BMC prior parameter list to be used for \code{"mapping = circle.opt"}
#' (if not provided, uses the default parameters from function \code{bmc.prior.default}):
#' * \code{c1}, \code{c2}, \code{c3} the empirical parameters
#' * \code{sigma} the standard deviation
#' * \code{kth}   the \out{k<sup>th</sup>} nearest seismic station used for
#' distance calculation
#' * \code{support} the information supporting the prior model
#' @param dist.calc the method to be used to evaluate distances (if not provided,
#' \code{dist.calc = "fast"}; read Details of function \code{d.geogr2km})
#' @param n.bootstrap the number of bootstraps (if not provided, \code{n.bootstrap = 0})
#' @return The data frame of 4 parameters:
#' * \code{lon} the longitude of the cell center
#' * \code{lat} the latitude of the cell center
#' * \code{mc.obs}  the completeness magnitude \out{m<sub>c</sub>} observed in the cell
#' * \code{sigma.obs}  the \out{m<sub>c</sub>} standard error in the cell
#' @references Mignan, A., Werner, M.J., Wiemer, S., Chen, C.-C., Wu, Y.-M. (2011),
#' Bayesian Estimation of the Spatially Varying Completeness Magnitude of Earthquake
#' Catalogs, Bull. Seismol. Soc. Am., 101, 1371-1385,
#' \href{https://pubs.geoscienceworld.org/bssa/article-lookup/101/3/1371}{doi: 10.1785/0120100223}
#' @seealso \code{bmc.prior}; \code{bmc.prior.default}; \code{d.geogr2km}; \code{mc.val}
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
#' # map mc with cylinder smoothing (this may take a few minutes!)
#' mc.cst <- mc.geogr(seism, "mbass", "circle.cst", dbin = 0.1, R = 20)
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
  R = 30, stations = NULL, kth = 4, params = NULL, dist.calc = "fast", n.bootstrap = 0) {
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
#   r <- sapply(1:grid.n, function(i) d.geogr2km(grid[i,], seism, method = dist.calc))
#   sapply not used to avoid reaching virtual memory limit for very large matrices
    r <- array(NA, dim = c(nrow(seism), grid.n))
    for(i in 1:grid.n) r[,i] <- d.geogr2km(grid[i,], seism, method = dist.calc)
    ind.cell <- sapply(1:grid.n, function(i) which(r[,i] <= R))
  }
  if (mapping == "circle.opt") {
    #   r <- sapply(1:grid.n, function(i) d.geogr2km(grid[i,], seism, method = dist.calc))
    #   sapply not used to avoid reaching virtual memory limit for very large matrices
    r <- array(NA, dim = c(nrow(seism), grid.n))
    for(i in 1:grid.n) r[,i] <- d.geogr2km(grid[i,], seism, method = dist.calc)
    d <- sapply(1:grid.n, function(i) d.geogr2km(grid[i,], stations, method = dist.calc))

    d.kth <- sapply(1:grid.n, function(i) sort(d[,i])[kth])
    if(is.null(params)) {
      params <- bmc.prior.default(kth)
    }
    mc.pred <- (params$c1 * d.kth ^ params$c2 + params$c3)
    dlow <- (((mc.pred - params$sigma) - params$c3) / params$c1) ^ (1 / params$c2)
    dhigh <- (((mc.pred + params$sigma) - params$c3) / params$c1) ^ (1 / params$c2)
    R <- (dhigh-dlow)/2
    Rmin <- sqrt(2) * dbin / 2 * 111
    R[which(R <= Rmin)] <- Rmin
    ind.cell <- sapply(1:grid.n, function(i) which(r[,i] <= R))
  }

  if(n.bootstrap > 0) {
    mc.unc <- sapply(1:grid.n, function(i) if(length(unlist(ind.cell[i])) >= nmin) {
      sapply(1:n.bootstrap, function(j) mc.val(
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
#' \out{k<sup>th</sup>} values (only \code{kth = 3}, \code{4} and \code{5}
#' allowed, otherwise returns \code{NULL}).
#'
#' The BMC default model is defined as the prior model evaluated for the Taiwan earthquake data
#' (Mignan et al., 2011), as it represents the best constrained data set so far (read more
#' on this in Mignan and Chouliaras, 2014). It often provides a better prior once calibrated
#' to other data than a new fit to the data (see example given for the function
#' \code{bmc.prior}). It is also used for rapid mapping of the optimal
#' \out{m<sub>c</sub>} map (see function \code{mc.geogr}).
#'
#' @param kth the \out{k<sup>th</sup>} nearest seismic station used for
#' distance calculation (if not provided, \code{kth = 4})
#' @return the BMC prior parameter list
#' * \code{c1}, \code{c2}, \code{c3} the empirical parameters
#' * \code{sigma} the standard error
#' * \code{kth}   the \out{k<sup>th</sup>} nearest seismic station used for
#' distance calculation
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
#' al., 2011) by fitting a function of the form
#' \ifelse{html}{\out{m<sub>c</sub> = c<sub>1</sub>d^c<sub>2</sub>+c<sub>3</sub>}}{\eqn{m_c = c_1 d^{c_2} +c_3}}
#' the observed completeness magnitude per cell and d the distance from the cell center to
#' the \out{k<sup>th</sup>} nearest seismic station.
#'
#' \code{support = "calibrated"} uses the default BMC prior model defined by the function
#' \code{bmc.prior.default} and calibrates it to the \code{mc.obs} data by shifting the
#' residual average to 0, substracting it from \out{c<sub>3</sub>} (see e.g., Mignan et al., 2013; Mignan
#' and Chouliaras, 2014).
#'
#' \code{support = "data"} directly fits the prior function to the \code{mc.obs} data
#' by using the Nonlinear Least Squares function \code{stats::nls}. If the estimation fails,
#' \code{support = "calibrated"} is used instead.
#'
#' @param mc.obs a data frame of the observed \out{m<sub>c</sub>} map defined by the
#' function \code{mc.geogr}
#' @param stations the seismic network data frame of parameters:
#' * \code{lon} the seismic station longitude
#' * \code{lat} the seismic station latitude
#' * \code{...} other station attributes
#' @param kth the \out{k<sup>th</sup>} nearest seismic station used for distance
#' calculation (if not provided, \code{kth = 4})
#' @param support the information supporting the prior model: \code{"calibrated"} or
#' \code{"data"} (if not provided, \code{support = "calibrated"} - see Details)
#' @param dist.calc the method to be used to evaluate distances (if not provided,
#' \code{dist.calc = "fast"}; read Details of function \code{d.geogr2km})
#' @return A list of:
#' @return the BMC prior model parameter list:
#' * \code{c1}, \code{c2}, \code{c3} the empirical parameters
#' * \code{sigma} the standard deviation
#' * \code{kth}   the \out{k<sup>th</sup>} nearest seismic station used for distance calculation
#' * \code{support}   the information supporting the prior model
#' @return the input data frame:
#' * \code{mc}   the completeness magnitude value per cell
#' * \code{d.kth}   the distance to the \out{k<sup>th</sup>} nearest seismic station per
#' cell (in km)
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
#' @seealso \code{bmc}; \code{bmc.prior.default}; \code{d.geogr2km}; \code{mc.geogr}
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
#' # map the observed mc (this may take a few minutes)
#' mc.obs <- mc.geogr(seism, "mode", "grid", dbin = 0.1)
#'
#' # test the two possible priors & plot
#' model.calibrated <- bmc.prior(mc.obs, stations, kth = 5)
#' model.fromdata <- bmc.prior(mc.obs, stations, kth = 5, support = "data")
#' data <- model.calibrated[[2]]
#' di <- seq(0, max(data$d.kth))
#' params.cal <- model.calibrated[[1]]
#' params.dat <- model.fromdata[[1]]
#' plot(data$d.kth, data$mc.obs)
#' lines(di, params.cal$c1*di^params.cal$c2+params.cal$c3, col="orange")
#' lines(di, params.dat$c1*di^params.dat$c2+params.dat$c3, col="red")
bmc.prior <- function(mc.obs, stations, kth = 4, support = "calibrated", dist.calc = "fast") {
  grid.n <- nrow(mc.obs)
  d <- sapply(1:grid.n, function(i) d.geogr2km(grid[i,], stations, method = dist.calc))
  d.kth <- sapply(1:grid.n, function(i) sort(d[,i])[kth])
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

  return(list(params, data.frame(mc.obs, d.kth = d.kth)))
}

#' BMC Bayesian Map
#'
#' Compute the map of the posterior completeness magnitude and posterior standard error
#' based on the maps of the observed and predicted completeness magnitudes.
#'
#' This is the final step of the Bayesian Magnitude of Completeness (BMC) method (Mignan et
#' al., 2011).
#'
#' @param mc.obs a data frame of 4 parameters defined by the function \code{mc.geogr}:
#' * \code{lon} the longitude of the cell center
#' * \code{lat} the latitude of the cell center
#' * \code{mc.obs}  the completeness magnitude \out{m<sub>c</sub>} in the cell
#' * \code{sigma.obs}  the \out{m<sub>c</sub>} standard error in the cell
#' @param mc.pred a vector of \out{m<sub>c</sub>} predicted values per cell
#' @param sigma.pred a vector of the BMC prior model standard deviation, repeated for all cells
#' @return The data frame of 8 parameters:
#' * \code{lon} the longitude of the cell center
#' * \code{lat} the latitude of the cell center
#' * \code{mc.obs}  the observed \out{m<sub>c</sub>} in the cell
#' * \code{sigma.obs}  the observed standard error in the cell
#' * \code{mc.pred}  the predicted \out{m<sub>c</sub>} in the cell
#' * \code{sigma.pred}  the prior model standard deviation in the cell
#' * \code{mc.post}  the posterior \out{m<sub>c</sub>} in the cell
#' * \code{sigma.post}  the posterior standard error in the cell
#' @references Mignan, A., Werner, M.J., Wiemer, S., Chen, C.-C., Wu, Y.-M. (2011),
#' Bayesian Estimation of the Spatially Varying Completeness Magnitude of Earthquake
#' Catalogs, Bull. Seismol. Soc. Am., 101, 1371-1385,
#' \href{https://pubs.geoscienceworld.org/bssa/article-lookup/101/3/1371}{doi: 10.1785/0120100223}
#' @seealso \code{bmc}; \code{bmc.prior}; \code{bmc.prior.default}; \code{mc.geogr}
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
#'
#' # reduce catalogue size for faster computation
#' seism <- subset(seism, yr >= 2008)
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
#' # map the observed mc & predicted mc (quick & dirty)
#' mc.obs <- mc.geogr(seism, "mode", "grid", dbin = 0.1, n.bootstrap = 100)
#' prior <-  bmc.prior(mc.obs, stations)
#' mc.pred <- (prior[[1]]$c1 * prior[[2]]$d.kth ^ prior[[1]]$c2 + prior[[1]]$c3)
#' sigma.pred <- rep(prior[[1]]$sigma, nrow(mc.obs))
#' res <- bmc.bayes(mc.obs, mc.pred, sigma.pred)
#'
#' #display the 6 BMC maps (mc.obs, mc.pred, mc.post, sigma.obs, sigma.pred, sigma.post)
#' image(matrix(res$mc.obs, nrow=length(unique(res$lon)), ncol=length(unique(res$lat))))
#' image(matrix(res$mc.pred, nrow=length(unique(res$lon)), ncol=length(unique(res$lat))))
#' image(matrix(res$mc.post, nrow=length(unique(res$lon)), ncol=length(unique(res$lat))))
#' image(matrix(res$sigma.obs, nrow=length(unique(res$lon)), ncol=length(unique(res$lat))))
#' image(matrix(res$sigma.pred, nrow=length(unique(res$lon)), ncol=length(unique(res$lat))))
#' image(matrix(res$sigma.post, nrow=length(unique(res$lon)), ncol=length(unique(res$lat))))
bmc.bayes <- function(mc.obs, mc.pred, sigma.pred) {
  grid.n <- nrow(mc.obs)
  mc.post <- sapply(1:grid.n, function(i) (mc.pred[i] * mc.obs$sigma.obs[i] ^ 2 +
    mc.obs$mc.obs[i] * sigma.pred[i] ^ 2) / (sigma.pred[i] ^ 2 + mc.obs$sigma.obs[i] ^ 2))
  sigma.post <- sapply(1:grid.n, function(i) sqrt((sigma.pred[i] ^ 2 * mc.obs$sigma.obs[i] ^ 2)
    / (sigma.pred[i] ^ 2 + mc.obs$sigma.obs[i] ^ 2)))

  indna <- which(is.na(mc.post) == T)
  mc.post[indna] <- mc.pred[indna]
  sigma.post[indna] <- sigma.pred[indna]
  return(data.frame(mc.obs, mc.pred, sigma.pred,
                    mc.post = unlist(mc.post), sigma.post = unlist(sigma.post)))
}


#' BMC Wrapper
#'
#' Run all the steps of the Bayesian Magnitude of Completeness (BMC) method (Mignan et al.,
#' 2011) and produce a spatial data frame of the completeness magnitude (observed,
#' predicted, and posterior) and associated uncertainties (observed, predicted, and
#' posterior).
#'
#' It is a wrap-up of other functions. See Examples of the function \code{bmc.bayes} for
#' a possible break-down of the different steps of the BMC method.
#'
#' The \code{support = "fast"} approach is the only one provided for the BMC wrapper so far.
#' It consists in estimating the optimal observed completeness magnitude \out{m<sub>c</sub>}
#' by directly using the default BMC prior model. The model is then calibrated to the
#' optimal observed \out{m<sub>c</sub>}. Finally, the Bayesian method is applied. This fast
#' approach was successfully tested in a number of regions (e.g., Kraft et al., 2013;
#' Mignan et al., 2013; Mignan and Chouliaras, 2014; Tormann et al., 2014; Panzera et al., 2017).
#'
#' @param seism an earthquake catalog data frame of parameters:
#' * \code{lon} the earthquake longitude
#' * \code{lat} the earthquake latitude
#' * \code{m}   the earthquake magnitude
#' * \code{...} other earthquake parameters
#' @param stations the seismic network data frame of parameters:
#' * \code{lon} the seismic station longitude
#' * \code{lat} the seismic station latitude
#' * \code{...} other station attributes
#' @param support the information supporting the BMC method: only  \code{"fast"} available
#' so far)
#' @param mbin the magnitude binning value (if not provided, \code{mbin = 0.1})
#' @param box a vector of the minimum longitude, maximum longitude, minimum latitude and
#' maximum latitude, in this order (if not provided, \code{box} is calculated from
#' the geographical limits of \code{seism})
#' @param dbin the spatial binning value (if not provided, \code{dbin} is calculated such that
#' the map is made of 10 longitudinal cells based on \code{box})
#' @param kth the \out{k<sup>th</sup>} nearest seismic station used for distance calculation
#' (if not provided, \code{kth = 4})
#' @return The data frame of 8 parameters:
#' * \code{lon} the longitude of the cell center
#' * \code{lat} the latitude of the cell center
#' * \code{mc.obs}  the observed \out{m<sub>c</sub>} in the cell
#' * \code{sigma.obs}  the observed standard error in the cell
#' * \code{mc.pred}  the predicted \out{m<sub>c</sub>} in the cell
#' * \code{sigma.pred}  the prior model standard deviation in the cell
#' * \code{mc.post}  the posterior \out{m<sub>c</sub>} in the cell
#' * \code{sigma.post}  the posterior standard error in the cell
#' @references Kraft, T., Mignan, A., Giardini, D. (2013), Optimization of a large-scale
#' microseismic monitoring network in northern Switzerland, Geophys. J. Int., 195, 474-490,
#' \href{https://academic.oup.com/gji/article/195/1/474/601101}{doi: 10.1093/gji/ggt225}
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
#' @references Panzera, F., Mignan, A., Vogfjord, K.S. (2017), Spatiotemporal evolution of
#' the completeness magnitude of the Icelandic earthquake catalogue from 1991 to 2013, J.
#' Seismol., 21, 615-630,
#' \href{https://link.springer.com/article/10.1007/s10950-016-9623-3}{doi: 10.1007/s10950-016-9623-3}
#' @references Tormann, T., Wiemer, S., Mignan, A. (2014), Systematic survey of
#' high-resolution b value imaging along Californian faults: inference on asperities, J.
#' Geophys. Res. Solid Earth, 119, 2029-2054,
#' \href{http://onlinelibrary.wiley.com/doi/10.1002/2013JB010867/full}{doi: 10.1002/2013JB010867}
#' @seealso \code{bmc.bayes}; \code{bmc.prior}; \code{bmc.prior.default}; \code{mc.geogr}
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
#'
#' # reduce catalogue size for faster computation
#' seism <- subset(seism, yr >= 2008)
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
#' # Apply the BMC method (this may take a few minutes)
#' res <- bmc(seism, stations, dbin = 0.1)
#'
#' #display the 6 BMC maps (mc.obs, mc.pred, mc.post, sigma.obs, sigma.pred, sigma.post)
#' image(matrix(res$mc.obs, nrow=length(unique(res$lon)), ncol=length(unique(res$lat))))
#' image(matrix(res$mc.pred, nrow=length(unique(res$lon)), ncol=length(unique(res$lat))))
#' image(matrix(res$mc.post, nrow=length(unique(res$lon)), ncol=length(unique(res$lat))))
#' image(matrix(res$sigma.obs, nrow=length(unique(res$lon)), ncol=length(unique(res$lat))))
#' image(matrix(res$sigma.pred, nrow=length(unique(res$lon)), ncol=length(unique(res$lat))))
#' image(matrix(res$sigma.post, nrow=length(unique(res$lon)), ncol=length(unique(res$lat))))
bmc <- function(seism, stations, support = "fast", mbin = 0.1, box = NULL, dbin = NULL, kth = 4) {
  if(support == "fast") {
    cat("Compute observed mc map (optimized)", "\n")
    mc.obs <- mc.geogr(seism, "mode", "circle.opt", mbin = mbin, box = box, dbin = dbin,
                       nmin = 4, R = NULL, stations = stations, kth = kth, n.bootstrap = 200)

    cat("Compute predicted mc map (calibrated default prior model)", "\n")
    prior <- bmc.prior(mc.obs, stations, kth = kth, support = "calibrated")
    mc.pred <- (prior[[1]]$c1 * prior[[2]]$d.kth ^ prior[[1]]$c2 + prior[[1]]$c3)
    sigma.pred <- rep(prior[[1]]$sigma, nrow(mc.obs))

    cat("Compute posterior mc map (combining both observed & predicted mc maps)", "\n")
    bmc.res <- bmc.bayes(mc.obs, mc.pred, sigma.pred)
    return(bmc.res)
  }
}




