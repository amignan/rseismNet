#' Map the completeness magnitude in space
#'
#' Map the completeness magnitude mc in space for a seismicity data frame, following
#' one of two possible geographical mappings: "grid" (computes mc in cells of size dbin)
#'  or "bmc" (computes mc according to the Bayesian Magnitude of Completeness method of
#'  Mignan et al., Bull. Seismol. Soc. Am., 2011).
#'
#' @param m Earthquake magnitude vector
#' @param lon Earthquake longitude vector
#' @param lat Earthquake latitude vector
#' @return The completeness magnitude mc data frame
#' @examples
#' mbin <- 0.1
#' url <- "http://service.scedc.caltech.edu/ftp/catalogs/"
#' dat <- "hauksson/Socal_DD/hs_1981_2016_comb_K4_A.cat_so_SCSN_v2q"
#' seism <- scan(paste(url, dat, sep = ""), what = "character", sep = "\n")
#' yr <- as.numeric(substr(seism, start=1, stop=4))
#' lat <- as.numeric(substr(seism, start=35, stop=42))
#' lon <- as.numeric(substr(seism, start=44, stop=53))
#' mag <- round(as.numeric(substr(seism, start=63, stop=67)), digits = log10(1/mbin))
#' seism <- data.frame(yr = yr, lon = lon,lat = lat, mag = mag)
#' seism <- subset(seism, yr >= 2010)
#' mc <- mc.geogr(seism, mbin, "mode", "grid", 0.1)
#' image(matrix(mc$mc, nrow=length(unique(mc$lon)), ncol=length(unique(mc$lat))))
mc.geogr <- function(seism, mbin, method, mapping, dbin = NULL, box = NULL, nmin = 4) {
  if(is.null(box)) box <- c(floor(min(seism$lon)), ceiling(max(seism$lon)),
                            floor(min(seism$lat)), ceiling(max(seism$lat)))
  if(is.null(dbin)) dbin <- (box[2]-box[1])/5

  grid <- expand.grid(lon = seq(box[1], box[2], dbin), lat = seq(box[3], box[4], dbin))
  grid.n <- nrow(grid)

  if (mapping == "grid") {
    ind.cell <- sapply(1:grid.n, function(i) which(seism$lon >= grid$lon[i] - dbin / 2 &
                                                   seism$lon < grid$lon[i] + dbin / 2 &
                                                   seism$lat >= grid$lat[i] - dbin / 2 &
                                                   seism$lat < grid$lat[i] + dbin / 2))
    mc.cell <- sapply(1:grid.n, function(i) if(length(unlist(ind.cell[i])) >= nmin) 
      { mc.val(seism$m[unlist(ind.cell[i])], mbin, method) } else { NA })
    return(data.frame(grid, mc=unlist(mc.cell)))
  }
  if (mapping == "bmc") {
    NULL
  }
}

