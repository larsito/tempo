#' Calculate syncrhony indices
#'
#' Function to calculate all or a subset of synchrony indices, including their
#' three term local quadrat variance (t3) version, and the decomposition of the
#' variance ratio index into Stotal, Strend, and Sdetrended.
#'
#' @param x temporal community data in a data frame. Species as columns, years
#'   as rows.
#' @param decompose If TRUE, the outputted synchrony indices will also contain
#'   the total variance ratio (Stotal) and its decomposition into parts
#'   attributable to long term trend (Strend) and actual synchrony (Sdetrended).
#'   See function decompostr.
#' @param indices A character vector that specifies wich of the synchrony
#'   indices should be calculated. By default, all available indices are
#'   calculated. If only a subset is desired, the names of the functions to
#'   calculate must be provided in the character vector (e.g. c("eta", "phi") to
#'   calculate only standard versions of eta and phi).
#' @return A data frame containing all indices that were specified to be
#' calculated by the indices and decompose arguments.
#' @seealso \code{\link{eta}}, \code{\link{eta_w}}, \code{\link{eta_t3}},
#'   \code{\link{eta_t3_w}}, \code{\link{phi}} , \code{\link{phi_t3}},
#'   \code{\link{varrat}}, \code{\link{varrat_t3}}, \code{\link{log_varrat}},
#'   \code{\link{log_varrat_t3}}, \code{\link{decompostr}}
#' @export calc_sync
calc_sync <- function(x,
                      decompose = TRUE,
                      indices = c(
                        "eta",
                        "eta_t3",
                        "eta_w",
                        "eta_t3_w",
                        "phi",
                        "phi_t3",
                        "varrat",
                        "varrat_t3",
                        "log_varrat",
                        "log_varrat_t3"
                      )) {
  if(any(colSums(x) == 0)) {
    warning("Species with only zero abundances were removed")
    x <- x[ , colSums(x) != 0]
  }
  if(any(apply(x, 2, sd) == 0)) {
    warning("At least one species has constant abundance through time. Eta can not be computed")
  }
  # fun_indices <- sapply(indices, match.fun, USE.NAMES = TRUE)
  fun_indices <- sapply(indices, get, envir = environment(syngenr))
  calc <- function(index) {
    sync <- index(x)
    sync <- data.frame(obs = sync)
    return(sync)
  }
  out <- lapply(fun_indices, function(f, x) {
    calc(f)
  }, x)
  out <- do.call("cbind", out)
  names(out) <- indices
  if(decompose == TRUE) {
    deco <- tempo:::decompostr(x)
    out <- data.frame(out, deco)
  }
  return(out)
}

#' Decomposition of the variance ratio synchrony measure
#'
#' \code{decompostr} decomposes the variance ratio into the part that is
#' produced by linear longterm trends in the data, and the part that is due to
#' actual synchronous or asynchronous fluctuations of species.
#'
#' @param x temporal community data in a data frame. Species as columns, years
#'   as rows.
#'
#' @return A data frame of one row with three variables. The total variance
#'   ratio, the part ot the total variance ration that is attributable to the
#'   long term trend, and the part that is attributable to actual synchrony
#'   between species (i.e. Stotal, Strend, Sdetrended)
#' @export decompostr
decompostr <- function(x) {
  time <- 1:nrow(x)
  repre <- apply(x, 2, function(x, y) { #save fitted and residuls of regressions
    mod <- lm(y ~ x)
    fit <- fitted(mod)
    resi <- resid(mod)
    out <- list(fit, resi)
    return(out)
  }, x = time)
  pred_mat <- t(do.call("rbind", unlist(lapply(repre, "[", 1),
                                        recursive = FALSE)))
  resi_mat <- t(do.call("rbind", unlist(lapply(repre, "[", 2),
                                        recursive = FALSE)))
  l <- list(dat_mat = x,
            pred_mat = pred_mat,
            resi_mat = resi_mat)
  vars <- sapply(l, function(x) {# calculate components
    sumvar <- sum(apply(x, 2, var))
    varsum <- var(rowSums(x))
    sumcovar <- varsum - sumvar
    out <- c(varsum = varsum,
             sumvar = sumvar,
             sumcovar = sumcovar)
    return(out)
  })
  s_tot <- vars[3, 1] / vars[2, 1]
  s_tr <- vars[3, 2] / vars[2, 1]
  s_resi <- vars[3, 3] / vars[2, 1]
  out <-
    data.frame(syn_total = s_tot,
               syn_trend = s_tr,
               syn_detrend = s_resi)
  return(out)
}

#' Function to calculate different versions of the eta synchrony measure (Gross
#' et al. 2014)
#'
#' \code{cor_algo} calculates standard and t3, as well as abundance weighted and
#' unweighted versions of eta. All these are based on correlations of abundances
#' between species. There is a series of helper functions, \code{eta},
#' \code{eta_w}, \code{eta_t3}, and \code{eta_t3_w}, that makes the different
#' versions of eta available in single functions, by predefined argument values.
#'
#' @param x temporal community data in a data frame. Species as columns, years
#'   as rows.
#' @param method the correlation method to be used. This can be either cor, or
#'   cor_t3, depending if correlations should be based on three term local
#'   variance.
#' @param weighted logical. If TRUE, the averaging over the species will be
#'   weighted by their abundance. The default is FALSE.
#' @return a single numeric value of the calculated version of eta.
#' @references Gross, K., B. J. Cardinale, J. W. Fox, A. Gonzalez, M. Loreau, H.
#'   W. Polley, P. B. Reich, and J. van Ruijven. 2014. Species richness and the
#'   temporal stability of biomass production: a new analysis of recent
#'   biodiversity experiments. Am Nat 183:1-12.
#' @references Bluthgen, N., N. K. Simons, K. Jung, D. Prati, S. C. Renner, S.
#'   Boch, M. Fischer, N. Holzel, V. H. Klaus, T. Kleinebecker, M. Tschapka, W.
#'   W. Weisser, and M. M. Gossner. 2016. Land use imperils plant and animal
#'   community stability through changes in asynchrony rather than diversity.
#'   Nat Commun 7:10697.
#' @export cor_algo
cor_algo <- function(x, method, weighted = FALSE) {
  n <- ncol(x)
  cor_sp <- vector("numeric", length = n)
  for (i in 1:n) {
    xi <- as.numeric(x[, i])
    xnoti <- rowSums(x[,-i, drop = FALSE])
    cor_sp[i] <- method(xi, xnoti)
  }
  if (weighted) {
    rel_ab <- colSums(x) / sum(x)
    sync <- sum(cor_sp * rel_ab)
  } else {
    sync <- mean(cor_sp)
  }
  return(sync)
}

#' Standard, weighted, and t3 versions of eta
#'
#' Short functions that make use of cor_algo to enable calculating all four
#' versions of eta (standard, based on t3, and both either weighted by species
#' abundance or not) in separate functions, making use of the more generic
#' \code{cor_algo}.
#' @param x temporal community data in a data frame. Species as columns, years
#'   as rows.
#' @return a single numeric value of the calculated version of eta.
#' @references Gross, K., B. J. Cardinale, J. W. Fox, A. Gonzalez, M. Loreau, H.
#'   W. Polley, P. B. Reich, and J. van Ruijven. 2014. Species richness and the
#'   temporal stability of biomass production: a new analysis of recent
#'   biodiversity experiments. Am Nat 183:1-12.
#' @references Bluthgen, N., N. K. Simons, K. Jung, D. Prati, S. C. Renner, S.
#'   Boch, M. Fischer, N. Holzel, V. H. Klaus, T. Kleinebecker, M. Tschapka, W.
#'   W. Weisser, and M. M. Gossner. 2016. Land use imperils plant and animal
#'   community stability through changes in asynchrony rather than diversity.
#'   Nat Commun 7:10697.
eta <- function(x) {
  tempo:::cor_algo(x, method = cor)
}

#' @rdname eta
eta_t3 <- function(x) {
  tempo:::cor_algo(x, method = tempo:::cor_t3)
}

#' @rdname eta
eta_w <- function(x) {
  tempo:::cor_algo(x, method = cor, weighted = TRUE)
}

#' @rdname eta
eta_t3_w <- function(x) {
  tempo:::cor_algo(x, method = tempo:::cor_t3, weighted = TRUE)
}

#' Simple mean pairwise correlation between species as a measure of synchrony
#'
#' @param x temporal community data in a data frame. Species as columns, years
#'   as rows.
#' @return a single numeric value of the calculated mean correlation.
#'   @export mean_cor
mean_cor <-  function(x) {
  cormat <- cor(x)
  mcor <- mean(cormat[lower.tri(cormat)], na.rm = TRUE)
  return(mcor)
}

#' Varince ratio index of synchrony
#'
#' \code{varrat} Calculates the varince ratio index of synchrony, its log
#' transformed version, as well as t3 versions of both.
#'
#' @param x temporal community data in a data frame. Species as columns, years
#'   as rows.
#' @return a single numeric value of the variance ration index of the temporal
#' community.
#' @references Schluter, D. 1984. A Variance Test for Detecting Species
#'   Associations, with Some Example Applications. ECOLOGY 65:998-998.
#' @references Leps, J., M. Majekova, A. Vitova, J. Dolezal, and F. de Bello.
#'   2018. Stabilizing effects in temporal fluctuations: management, traits, and
#'   species richness in high-diversity communities. ECOLOGY 99:360-371.
#' @export varrat
varrat <- function(x, log_trans = FALSE, t3 = FALSE) {
  if (t3 == FALSE) {
    all_cov <- cov(x, use = "pairwise.complete.obs")
    col_var <- diag(all_cov)
  } else {
    all_cov <- tempo:::var_cov_mat_t3(comdat)
    col_var <- diag(all_cov)
  }
  com_var <- sum(all_cov)
  pop_var <- sum(col_var)
  vr <- com_var / pop_var
  if (log_trans == TRUE) {
    vr <- log(vr)
  }
  return(vr)
}

#' Single functions for variance ratio, the log transformed version, and the
#' respective t3 versions
#'
#' Set of functions that uses function \code{varrat} and prespecified arguments
#' within it to separately calculate the variance ratio as its log, its t3
#' version, and the log of the t3 version.
#'
#' @param x temporal community data in a data frame. Species as columns, years
#'   as rows.
#' @return a single numeric value of the variance ratio index.
#' @references Schluter, D. 1984. A Variance Test for Detecting Species
#'   Associations, with Some Example Applications. ECOLOGY 65:998-998.
#' @references Leps, J., M. Majekova, A. Vitova, J. Dolezal, and F. de Bello.
#'   2018. Stabilizing effects in temporal fluctuations: management, traits, and
#'   species richness in high-diversity communities. ECOLOGY 99:360-371.
#' @export log_varrat
#' @export varrat_t3
#' @export log_varrat_t3
log_varrat <- function(x) {
  vr <- tempo:::varrat(x, log_trans = TRUE, t3 = FALSE)
  return(vr)
}

#' @rdname log_varrat
varrat_t3 <- function(x) {
  vr <- tempo:::varrat(x, log_trans = TRUE, t3 = FALSE)
  return(vr)
}

#' @rdname log_varrat
log_varrat_t3 <- function(x) {
  vr <-tempo:::varrat(x, log_trans = TRUE)
  return(vr)
}

#' Standard and t3 version of phi synchrony measure.
#'
#' Functions \code{phi} and \code{phi_t3} calculate phi (Loreau and de
#' Mazancourt 2008) in its standard version, and based on three term local
#' quadrat variance.
#'
#' @param x temporal community data in a data frame. Species as columns, years
#'   as rows.
#' @return a single numeric value of phi.
#' @references Loreau, M., and C. de Mazancourt. 2008. Species Synchrony and Its
#'   Drivers: Neutral and Nonneutral Community Dynamics in Fluctuating
#'   Environments. The American Naturalist 172:E48-E66.
#' @export phi
#' @export phi_t3
phi <- function (x) {
  species.sd = apply(x, MARGIN = 2, FUN = sd)
  community.var = var(rowSums(x))
  return(community.var / sum(species.sd, na.rm = TRUE) ^ 2)
}

#' @rdname phi
phi_t3 <- function (x) {
  species.var = apply(x, MARGIN = 2, FUN = var_t3)
  species.sd = sqrt(species.var)
  community.var = tempo:::var_t3(rowSums(x))
  return(community.var / sum(species.sd, na.rm = TRUE) ^ 2)
}

#' Variance, covariance, and correlations based on three term local quadrat
#' variance (t3)
#'
#' Functions to calculate variance, covariance, and correlations based on three
#' term local quadrat variance (t3).These are required to calculate indices
#' based on t3.
#'
#' @param x,y vectors containing a continous variable. In the context of
#'   temporal community data, these will be the abundances of one species (for
#'   var_t3), or of two species (for cov_t3 and cor_t3) through time.
#'
#' @return For var_t3, the three term local quadrat variance of a single
#'   variable, returned in a vector of length one. For cov_t3 and cor_t3,
#'   covariance and correlation, respectively between two variables.
#' @export var_t3
var_t3 <- function(x) {
  ws = 1
  n <- length(x) # sample size
  neg <- rep(c(1, -1), ws + 1)[1:(ws + 2)]
  multip  <- c(1, seq(ws + 1, 1)) * neg
  xsq <- vector("numeric", n - ws - 2)
  for (i in 1:(n - ws - 1)) {
    # calculate the "within brackets" component of the 3 term local variance for
    # xi, then the actual 3 term local quadratic variance of xi
    xisq <- x[i:(i + ws + 1)] * multip
    xsq[i] <- sum(xisq) ** 2 / 6
  }
  v <- sum(xsq) / (n - 2) # three term local quadratic variance of x
  return(v)
}

#' @rdname var_t3
#' @export cov_t3
cov_t3 <- function(x, y) {
  cov <- (tempo:::var_t3(x + y) - tempo:::var_t3(x) - tempo:::var_t3(y)) / 2
  return(cov)
}

#' @rdname var_t3
#' @export cor_t3
cor_t3 <- function(x, y) {
  cor <- tempo:::cov_t3(x, y) / sqrt(tempo:::var_t3(x) * tempo:::var_t3(y))
  return(cor)
}

#' Variance-covariance matrix based on three term local quadrat variance
#'
#' In the context of temporal communities, this will be the variance-covariance
#' matrix between the abundances of all species in a temporal community over
#' time.
#'
#' @param x temporal community data in a data frame. Species as columns, years
#'   as rows.
#' @return A matrix containing the t3 based variance-covariance matrix.
#' @export var_cov_mat_t3
var_cov_mat_t3 <- function(x) {
  nc <- ncol(x)
  S <- matrix(NA, nc, nc)
  for (i in 1:nc) {
    S[i, i] <- tempo:::var_t3(x[, i])
    for (j in (i + 1):nc) {
      if (i < nc) {
        S[i, j] <- tempo:::cov_t3(x[, i], x[, j])
        S[j, i] <- S[i, j]
      }
    }
  }
  return(S)
}
