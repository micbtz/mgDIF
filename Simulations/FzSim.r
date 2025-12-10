
# auxiliary functions used for simulations



sim_multidif <- function(r, itempar, ngroups, nitems, samplsize, percdif, dif_level, dif_type, dif_effect_unif, dif_effect, mean_gr, sd_gr, tau_mcp)
{
  set.seed(r)
  
  source("../functions_IndLik.r")
  Rcpp::sourceCpp("../IndLik.cpp")
  
  genData <- function(n, itempar, mean, sd)
  {
    nitems <- nrow(itempar)
    discrm <- itempar[, 1]
    diff <- itempar[, 2]
    rn <- runif(n*nitems)
    rnmat <- matrix(rn, n, nitems)
    latent <- rnorm(n, mean, sd)
    lp <- outer(latent, diff, FUN = "-")
    discrmMAT <- matrix(rep(discrm, each = n), nrow = n)
    lp <- lp*discrmMAT
    pr <- plogis(lp)
    dat <- (pr>rnmat)*1
    dat
  }
  
  if (dif_type=="u") dif_eff <- dif_effect_unif[dif_level, ]
  if (dif_type=="nu") dif_eff <- dif_effect[dif_level, ]
  
  repitems <- nitems*ngroups
  set.seed(r)
  difitems <- sample(1:repitems, size = percdif*repitems)
  difitems
  difitems_mat <- matrix(1:repitems %in% difitems, nrow = nitems)
  difitems_pos <- which(difitems_mat, arr.ind = T)
  
  itempar_gr <- list()
  for (i in 1:ngroups) itempar_gr[[i]] <- itempar
  if (percdif>0)
  {
    for (i in 1:nrow(difitems_pos))
    {
      it <- difitems_pos[i, 1]
      gr <- difitems_pos[i, 2]
      itempar_gr[[gr]][it, ] <- itempar_gr[[gr]][it, ] + dif_eff
    }
  }
  
  set.seed(r)
  datasets <- list()
  for (i in 1:ngroups) 
    datasets[[i]] <- genData(n = samplsize, itempar = itempar_gr[[i]], mean = mean_gr[i], sd = sd_gr[i])
  
  library(mirt)
  mods <- list()
  for (i in 1:ngroups)
    mods[[i]] <- mirt(datasets[[i]], 1, itemtype = "2PL", SE = TRUE, 
                      technical = list(NCYCLES=2000))
  
  library(equateMultiple)
  mod2pl <- modIRT(est.mods = mods)
  
  eqLIK1 <- multiec_lik_unc(mods = mod2pl, base = 1, se = TRUE, start = NULL, X = NULL, iter.max = 100000, trace = TRUE, check.all = FALSE) 
  
  lambda <- seq(1, 100, by = 0.5)
  out <- list()
  for (i in 1:length(lambda))
  {
    eqLIK5 <- multiec_lik_unc(mods = mod2pl, base = 1, se = FALSE, start = eqLIK1, X = NULL, iter.max = 100000, trace = TRUE, check.all = TRUE, lambda_mcp = lambda[i], tau_mcp = tau_mcp)
    out[[i]] <- eqLIK5[c("A", "B", "coef", "loglik", "delta_a", "delta_b", "npar", "AIC", "BIC")]
  }
  aic <- lapply(out, function(x) x$AIC)
  bic <- lapply(out, function(x) x$BIC)
  sel_aic <- which.min(aic)
  sel_bic <- which.min(bic)
  
  res_diftest <- try(dif.test(est.mods = mods, reference = 1, method = "Haebara", purification = TRUE))
  if(inherits(res_diftest, "try-error")) res_diftest <- NA
  else
  {
    res_diftest$equatings <- c()
    res_diftest$var_trasf <- c()
  }
  
  eqLIK1$varFull <- c()
  
  return(list(difitems = difitems, difitems_mat = difitems_mat, eqLIK1 = eqLIK1, lambda = lambda, out_aic = out[[sel_aic]], out_bic = out[[sel_bic]], sel_aic = sel_aic, sel_bic = sel_bic, res_diftest = res_diftest, difitems_pos = difitems_pos))
  
}


sim_multidif_pos <- function(r, itempar_groups, ngroups, samplsize, mean_gr, sd_gr, cf)
{
  set.seed(r)
  
  source("../functions_IndLik.r")
  Rcpp::sourceCpp("../IndLik.cpp")
  
  genData <- function(n, itempar, mean, sd)
  {
    nitems <- nrow(itempar)
    discrm <- itempar[, 1]
    diff <- itempar[, 2]
    rn <- runif(n*nitems)
    rnmat <- matrix(rn, n, nitems)
    latent <- rnorm(n, mean, sd)
    lp <- outer(latent, diff, FUN="-")
    discrmMAT <- matrix(rep(discrm, each = n), nrow = n)
    lp <- lp * discrmMAT
    pr <- plogis(lp)
    dat <- (pr > rnmat)*1
    dat
  }
  
  itempar_gr <- list()
  for (g in 1:ngroups) itempar_gr[[g]] <- itempar_groups[[g]]
  
  for (g in 1:ngroups)
  {
    ni <- nrow(itempar_gr[[g]])
    itempar_gr[[g]][, 1] <- itempar_gr[[g]][, 1] + cf * (1:ni)
    itempar_gr[[g]][, 2] <- itempar_gr[[g]][, 2] + cf * (1:ni)
  }
  
  set.seed(r)
  count <- 0
  zero_perfect_scores <- TRUE
  datasets <- list()
  for (g in 1:ngroups) 
  {
    datasets[[g]] <- genData(n = samplsize, itempar = itempar_gr[[g]], mean = mean_gr[g], sd = sd_gr[g])
    colsums <- colSums(datasets[[g]])
    zero_scores <- colsums == 0
    perfect_scores <- colsums == samplsize
    datasets[[g]] <- datasets[[g]][,(!zero_scores & !perfect_scores)]
  }
    library(mirt)
    mods <- list()
    for (i in 1:ngroups)
      mods[[i]] <- mirt(datasets[[i]], 1, itemtype = "2PL", SE = TRUE, 
                        technical = list(NCYCLES=2000))
    
    library(equateMultiple)
    mod2pl <- modIRT(est.mods = mods)
    
    eqLIK1 <- multiec_lik_unc(mods = mod2pl, base = 1, se = TRUE, start = NULL, X = NULL, iter.max = 100000, trace = TRUE, check.all = FALSE) 
    
    Xpos <- list()
    for (g in 1:ngroups)
    {
      ni <- ncol(datasets[[g]])
      Xpos[[g]] <- matrix(1:ni)
      rownames(Xpos[[g]]) <- colnames(datasets[[g]])
    }
    X <- list(Xpos)
    
    eqLIK2 <- multiec_lik_unc(mods = mod2pl, base = 1, se = TRUE, start = NULL, X = X, iter.max = 100000, trace = TRUE, check.all = FALSE) 
    
    eqLIK1$varFull <- c()
    eqLIK2$varFull <- c()
  return(list(eqLIK1 = eqLIK1, eqLIK2 = eqLIK2))
}


numDif <- function(x)
{
  if (any(is.na(x$res_diftest)))
    out <- NA
  else
  {
    difitems <- unique(x$difitems_pos[, 1])
    out <- length(difitems)
  }
  out
}

numNotDif <- function(x)
{
  if (any(is.na(x$res_diftest)))
    out <- NA
  else
  {
    difitems <- unique(x$difitems_pos[, 1])
    out <- 30 - length(difitems)
  }
  out
}


tp <- function(x)
{
  if (any(is.na(x$res_diftest)))
    out <- NA
  else
  {
    difitems <- unique(x$difitems_pos[, 1])
    detected <- which(x$res_diftest$test[, 2] < 0.05)
    out <- sum(detected %in% difitems)
  }
  out
}

fp <- function(x)
{
  if (any(is.na(x$res_diftest)))
    out <- NA
  else
  {
    difitems <- unique(x$difitems_pos[, 1])
    detected <- which(x$res_diftest$test[, 2] < 0.05)
    out <- sum(!detected %in% difitems)
  }
}


