
library(numDeriv)
library(Matrix)




multiec_lik_unc <- function(mods, base, se, start, iter.max, trace, X, check.all, lambda_mcp = 0, tau_mcp=3)
{
  if (trace)
    cat("Computation of equating coefficients . . . . \n")
  itms_all <- lapply(mods, function(x) names(x$coef)) # item labels
  itms <- unique(unlist(itms_all))
  itms <- sort(itms)
  num.forms <- length(mods)
  modsnames <- names(mods)
  tab <- data.frame(itms = itms)
  for (k in 1:num.forms) tab <- merge(tab, mods[[k]]$coef, by.x = "itms", by.y = 0, all = TRUE, suffixes = c(k - 1, k))
  colnames(tab)[-1] <- modsnames
  tab$itms <- as.character(tab$itms)
  rownames(tab) <- tab$itms
  
  tab_vars <- data.frame(itms = itms)
  for (k in 1:num.forms) tab_vars <- merge(tab_vars, diag(mods[[k]]$var), by.x = "itms", by.y = 0, all = TRUE, suffixes = c(k - 1, k))
  colnames(tab_vars)[-1] <- modsnames
  
  itmvar <- lapply(mods, FUN = function(x) x$var)
  itmvar <- lapply(itmvar, equateMultiple:::delGussng2)
  for (i in 1:length(itmvar)) rownames(itmvar[[i]]) <- colnames(itmvar[[i]]) <- paste(rownames(itmvar[[i]]), i, sep = ".")
  
  itmp <- 2
  if (sum(substr(tab$itms, 1, 6) == "Dscrmn") == 0)
    itmp = 1
  if (sum(substr(tab$itms, 1, 6) == "Gussng") > 0)
    itmp = 3
  
  bj1T <- as.matrix(tab[substr(tab$itms, 1, 6) == "Dffclt", ][, -1])
  if (itmp > 1)
    aj1T <- as.matrix(tab[substr(tab$itms, 1, 6) == "Dscrmn", ][, -1])
  else aj1T <- NULL
  if (itmp == 3)
    cj1T <- as.matrix(tab[substr(tab$itms, 1, 6) == "Gussng", ][, -1])
  else cj1T <- NULL
  
  var_bj1T <- as.matrix(tab_vars[substr(tab$itms, 1, 6) == "Dffclt", ][, -1])
  if (itmp > 1)
    var_aj1T <- as.matrix(tab_vars[substr(tab$itms, 1, 6) == "Dscrmn", ][, -1])
  else var_aj1T <- NULL
  
  if (itmp >= 2) 
  {
    if (is.null(start))
      ini <- c(rep(1, num.forms - 1), rep(0, num.forms - 1))
    if (inherits(start, "mlteqc"))
      ini <- c(start$A[-1], start$B[-1])
    namesAB <- c(paste("A", (1:num.forms)[-base], sep = "."), paste("B", (1:num.forms)[-base], sep = "."))
    names(ini) <- namesAB
  }
  if (itmp == 1) 
  {
    if (is.null(start))
      ini <- rep(1, num.forms - 1)
    if (inherits(start, "mlteqc"))
      ini <- start$B[-1]
    namesB <- paste("B", (1:num.forms)[-base], sep = ".")
    names(ini) <- namesB
  }
  
  if (!is.null(X))
  {
    Xmat <- list()
    itms_nam <- substr(rownames(bj1T), 8, 100)
    for (i in 1:length(X)) # loop on variable in X
    {
      Xmat[[i]] <- data.frame(itms_nam = itms_nam)
      for (k in 1:num.forms) Xmat[[i]] <- merge(Xmat[[i]], X[[i]][[k]], by.x = "itms_nam", by.y = 0, all = TRUE, suffixes = c(k - 1, k))
      rownames(Xmat[[i]]) <- Xmat[[i]]$itms_nam
      Xmat[[i]]$itms_nam <- c()
      colnames(Xmat[[i]]) <- modsnames
      Xmat[[i]] <- as.matrix(Xmat[[i]])
    }
    if (itmp >= 2) 
    {
      namesXvar <- c(paste("Dscrmn", names(X), sep = "."), paste("Dffclt", names(X), sep = "."))
      ini_X <- rep(0, length(X) * 2)
      names(ini_X) <- namesXvar
    }
    if (itmp == 1) 
    {
      namesXvar <- paste("Dffclt", names(X), sep = ".")
      ini_X <- rep(0, length(X))
      names(ini_X) <- namesXvar
    }
    ini <- c(ini, ini_X)
  }
  else Xmat <- NULL
  
  if (itmp >= 2)
    opt <- try(optim(par = ini, fn = ProfLik_unc_Rcpp, gr = derProfLik_unc_Rcpp, notbase = (1:num.forms)[-base] - 1, numforms = num.forms, aj1T = aj1T, bj1T = bj1T, var_aj1T = var_aj1T, var_bj1T = var_bj1T, Xmat = Xmat, hessian = se, method = "BFGS", control = list(maxit = iter.max)))
  if (itmp == 1)
    opt <- try(optim(par = ini, fn = ProfLik_unc_1PL_Rcpp, gr = derProfLik_unc_1PL_Rcpp, notbase = (1:num.forms)[-base] - 1, numforms = num.forms, bj1T = bj1T, var_bj1T = var_bj1T, Xmat = Xmat, hessian = se, method = "BFGS", control = list(maxit = iter.max)))
  
  A <- B <- se.A <- se.B <- varAB <- as <- bs <- se.as <- se.bs <- partial <- 
    coef <- se.coef <- loglik <- delta_a <- delta_b <- npar <- AIC <- BIC <- NULL
  
  if (!isa(opt, "try-error")) 
  {
    npar <- length(ini)
    aj1T_adj <- aj1T
    bj1T_adj <- bj1T
    delta_a <- aj1T*0
    delta_b <- bj1T*0
    if (check.all)
    {
      gamma_mcp_a <- var_aj1T*(1+tau_mcp)
      gamma_mcp_b <- var_bj1T*(1+tau_mcp)
      notconv <- TRUE
      count <- 0
      while (notconv)
      {
        count <- count + 1
        val <- opt$val
        if (itmp >= 2)
        res <- res_ProfLik_Rcpp(par = opt$par, notbase = (1:num.forms)[-base] - 1, numforms = num.forms, aj1T = aj1T_adj, bj1T = bj1T_adj, var_aj1T = var_aj1T, var_bj1T = var_bj1T, Xmat = Xmat)
        if (itmp == 1)
        res <- res_ProfLik_1PL_Rcpp(par = opt$par, notbase = (1:num.forms)[-base] - 1, numforms = num.forms, bj1T = bj1T_adj, var_bj1T = var_bj1T, Xmat = Xmat)
        if (itmp >= 2) 
        {
          drv_a <- fn_mcp_der(res = res$res_a, delta = delta_a, lambda_mcp = lambda_mcp, gamma = gamma_mcp_a, var = var_aj1T)
          min_drv_a<-min_pos(drv_a)
        }
        drv_b <- fn_mcp_der(res = res$res_b, delta = delta_b, lambda_mcp = lambda_mcp, gamma_mcp = gamma_mcp_b, var = var_bj1T)
        min_drv_b<-min_pos(drv_b)
        minval<-min(min_drv_a$val,min_drv_b$val)
        if (minval>0) notconv <- FALSE # none negative directional derivative, any parameter updated
        else
        {
          sel <- which.min(c(min_drv_a$val,min_drv_b$val))
          
          if (sel == 1) # selected to update a discrimination parameter
          {
            selit <- min_drv_a$ind
            delta_a_new <- mcp_sol(res = res$res_a[selit],delta=delta_a[selit], lambda_mcp = lambda_mcp, gamma_mcp = gamma_mcp_a[selit], var = var_aj1T[selit])
            delta_a[selit] <- delta_a_new
            aj1T_adj[selit] <- aj1T[selit] - delta_a_new
            opt <- try(optim(par = opt$par, fn = ProfLik_unc_Rcpp, gr = derProfLik_unc_Rcpp, notbase = (1:num.forms)[-base] - 1, numforms = num.forms, aj1T = aj1T_adj, bj1T = bj1T_adj, var_aj1T = var_aj1T, var_bj1T = var_bj1T, Xmat = Xmat, hessian = se, method = "BFGS", control = list(maxit = iter.max)))
          }
          
          if (sel == 2) # selected to update a difficulty parameter
          {
            selit <- min_drv_b$ind
            delta_b_new <- mcp_sol(res = res$res_b[selit],delta=delta_b[selit], lambda_mcp = lambda_mcp, gamma_mcp = gamma_mcp_b[selit], var = var_bj1T[selit])
            delta_b[selit] <- delta_b_new
            bj1T_adj[selit] <- bj1T[selit] - delta_b_new
            if (itmp >= 2)
              opt <- try(optim(par = opt$par, fn = ProfLik_unc_Rcpp, gr = derProfLik_unc_Rcpp, notbase = (1:num.forms)[-base] - 1, numforms = num.forms, aj1T = aj1T_adj, bj1T = bj1T_adj, var_aj1T = var_aj1T, var_bj1T = var_bj1T, Xmat = Xmat, hessian = se, method = "BFGS", control = list(maxit = iter.max)))
            if (itmp == 1)
              opt <- try(optim(par = opt$par, fn = ProfLik_unc_1PL_Rcpp, gr = derProfLik_unc_1PL_Rcpp, notbase = (1:num.forms)[-base] - 1, numforms = num.forms, bj1T = bj1T_adj, var_bj1T = var_bj1T, Xmat = Xmat, hessian = se, method = "BFGS", control = list(maxit = iter.max)))
          }
          if (max(abs(val-opt$val))<0.001) notconv <- FALSE
          if (count == 5000)
          {
            message("max iterations reached")
            notconv <- FALSE
          }          
        }
      }
      npar <- npar + sum(delta_b != 0, na.rm = TRUE)
      if (itmp >= 2) 
      {
        npar <- npar + sum(delta_a != 0, na.rm = TRUE)
      }
    }
    
    loglik <- -opt$value
    AIC <- -2 * loglik + 2 * npar 
    n <- sum(!is.na(bj1T))
    if (itmp >= 2) n <- n * 2
    BIC <- -2 * loglik + log(n) * npar
    A <- rep(1, num.forms)
    B <- rep(0, num.forms)
    if (itmp >= 2) 
    {
      A[-base] <- opt$par[1:(num.forms - 1)]
      B[-base] <- opt$par[num.forms:(2 * (num.forms - 1))]
    }
    if (itmp == 1) 
      B[-base] <- opt$par[1:(num.forms - 1)]
    names(A) <- names(B) <- modsnames
    conv <- opt$convergence
    
    if (!is.null(X))
    {
      if (itmp >= 2) coef <- opt$par[(2 * num.forms - 1) : length(opt$par)]
      if (itmp == 1) coef <- opt$par[num.forms : length(opt$par)]
      names(coef) <- namesXvar
    }
    else coef <- NULL
    
    # computation of item parameters on a common scale
    if (itmp >= 2)
    {
      asbs <- calc_asbs_Rcpp(par = opt$par, notbase = (1:num.forms)[-base] - 1, numforms = num.forms, aj1T = aj1T_adj, bj1T = bj1T_adj, var_aj1T = var_aj1T, var_bj1T = var_bj1T, Xmat = Xmat)
      as <- as.vector(asbs$as)
      bs <- as.vector(asbs$bs)
      names(as) <- row.names(aj1T)
      names(bs) <- row.names(bj1T)
    }
    if (itmp == 1)
    {
      asbs <- calc_bs_1PL_Rcpp(par = opt$par, notbase = (1:num.forms)[-base] - 1, numforms = num.forms, bj1T = bj1T, var_bj1T = var_bj1T, Xmat = Xmat)
      bs <- as.vector(asbs$bs)
      names(bs) <- row.names(bj1T)
    }
    
    bj1Ts <- bj1T
    for (t in 1:num.forms) bj1Ts[, t] <- bj1Ts[, t] * A[t] + B[t]
    if (itmp >= 2)
    {
      aj1Ts <- aj1T
      for (t in 1:num.forms) aj1Ts[, t] <- aj1Ts[, t] / A[t]
    }
    if (itmp == 1)
      tabs <- bj1Ts
    if (itmp == 2)
      tabs <- rbind(bj1Ts, aj1Ts)
    if (itmp == 3)
      tabs <- rbind(bj1Ts, aj1Ts, cj1T)
    colnames(tabs) <- paste(colnames(tabs), modsnames[base], sep = ".as.")
    tabs <- tabs[, -base]
    tab <- cbind(tab, tabs)
    colnames(tab)[1] <- "Item"
    
    partial <- NULL
    varAB <- se.as <- se.bs <- se.coef <- NULL
    se.A <- se.B <- rep(NA, num.forms)
    names(se.A) <- names(se.B) <- modsnames
    
    if (se) 
    {
      if (trace) cat("Computation of standard errors ")
      if (trace) cat(" . ")
      derS_AB <- opt$hessian
      tabnomib <- outer(rownames(bj1T), 1:num.forms, paste, sep = ".")
      tabnomib[is.na(bj1T)] <- NA
      nomib <- tabnomib[!is.na(tabnomib)]
      nomiab <- nomib
      if (itmp >= 2)
      {
        gamma <- c(aj1T, bj1T)
        tabnomia <- outer(rownames(aj1T), 1:num.forms, paste, sep = ".")
        tabnomia[is.na(aj1T)] <- NA
        nomia <- tabnomia[!is.na(tabnomia)]
        nomiab <- c(nomia, nomib)
      }
      if (itmp == 1) gamma <- bj1T
      gamma <- gamma[!is.na(gamma)]
      if (itmp >= 2)
        derS_gamma <- jacobian(func = derAB_lik_unc, x = gamma, par = opt$par, base = base, num.forms = num.forms, aj1T = aj1T, bj1T = bj1T, var_aj1T = var_aj1T, var_bj1T = var_bj1T, Xmat = Xmat, method = "simple")
      if (itmp == 1)
        derS_gamma <- jacobian(func = derAB_lik_1PL_unc, x = gamma, par = opt$par, base = base, num.forms = num.forms, bj1T = bj1T, var_bj1T = var_bj1T, Xmat = Xmat, method = "simple")
      if (trace) cat(" . ")
      colnames(derS_gamma) <- nomiab
      derAB_gamma <- -Matrix::solve(derS_AB, derS_gamma)
      partial <- t(derAB_gamma) # derivatives of A and B equating coefficients with respect to the item parameters
      if (itmp == 1)
        partial <- cbind(matrix(0, nrow(partial), num.forms - 1), partial)
      namesAB <- c(paste("A", (1:num.forms)[-base], sep = "."), paste("B", (1:num.forms)[-base], sep = "."))
      colnames(partial)[1:(num.forms * 2-2)] <- namesAB
      
      var_gamma <- Matrix::bdiag(itmvar)
      
      VarNames_list <- sapply(itmvar, function(x) rownames(x))
      VarNames <- unlist(VarNames_list)
      rownames(var_gamma) <- colnames(var_gamma) <- VarNames
      
      sel <- colnames(derAB_gamma)
      if (trace) cat(" . ")
      vargamma_tderAB_gamma <- Matrix::tcrossprod(var_gamma[sel, sel], derAB_gamma)
      varAB <- derAB_gamma %*% vargamma_tderAB_gamma
      if (trace) cat(" . \n")
      
      stderr <- Matrix::diag(varAB)^0.5
      
      se.A <- se.B <- rep(0, num.forms)
      names(se.A) <- names(se.B) <- modsnames
      stderrs <- Matrix::diag(varAB)^0.5
      if (itmp >= 2) 
      {
        se.A[-base] <- stderrs[1:(num.forms - 1)]
        se.B[-base] <- stderrs[(num.forms):(2 * num.forms - 2)]
      }
      if (!is.null(X))
      {
        if (itmp >= 2) se.coef <- stderrs[(2 * num.forms - 1) : length(opt$par)]
        if (itmp == 1) se.coef <- stderrs[num.forms : length(opt$par)]
        names(se.coef) <- namesXvar
      }
      
      if (itmp == 1) 
      {
        se.B[-base] <- stderrs[1:(num.forms - 1)]
        tmp <- matrix(0, (2 * num.forms - 2), (2 * num.forms - 2))
        tmp[num.forms:(2 * num.forms - 2), num.forms:(2 * num.forms - 2)] <- as.matrix(varAB)[1 : (num.forms - 1), 1 : (num.forms - 1)]
        varAB <- tmp
        colnames(varAB) <- rownames(varAB) <- namesAB
      }
      varAB <- as.matrix(varAB)
    }
  }
  else conv <- 101
  out <- list(A = A, B = B, se.A = se.A, se.B = se.B, varAB = varAB, as = as, bs = bs, 
              se.as = se.as, se.bs = se.bs, tab = tab, varFull = itmvar, 
              partial = partial, itmp = itmp, method = "IndLik", 
              basename = modsnames[base], convergence = conv, 
              coef = coef, se.coef = se.coef, loglik = loglik, 
              delta_a = delta_a, delta_b = delta_b, npar = npar, AIC = AIC, BIC = BIC)
  class(out) <- "mlteqc"
  return(out)
  
} 



min_pos<-function(x)
{
  mins<-sapply(x,min,na.rm=TRUE)
  sel<-which.min(mins)
  val<-mins[sel]
  ind<-which(x[[sel]]==val,arr.ind = TRUE)
  list(sel=sel,ind=ind,val=val)
}

soft_thresholding_hetero <- function(zeta, var, lambda, weights)
{
  lambda_var <- lambda * var * weights
  rdc <- abs(zeta) - lambda_var
  pospart <- pmax(rdc, 0)
  sign(zeta) * pospart
}



mcp_sol<-function(res,delta,lambda_mcp,gamma_mcp,var)
{
  zeta<-res+delta
  cond1<-abs(zeta) <= lambda_mcp*gamma_mcp
  if (cond1) out <- soft_thresholding_hetero(zeta=zeta,var=var,lambda=lambda_mcp,weights=1)*gamma_mcp/(gamma_mcp-var)
  if (!cond1) out <- zeta
  out
}

fn_mcp_der<-function(res, delta, lambda_mcp, gamma_mcp, var)
{
  dir_der_pos<- dir_der_neg<- matrix(NA,nrow(delta),ncol(delta))
  res_var <- res/var
  
  cond1<- (abs(delta) <= lambda_mcp*gamma_mcp) & !is.na(res)

  if (any(cond1)) dir_der_pos[cond1]<- -res_var[cond1] + lambda_mcp*sign0p(delta[cond1]) - delta[cond1]/gamma_mcp[cond1]
  if (any(cond1)) dir_der_neg[cond1]<-  res_var[cond1] + lambda_mcp*sign0m(delta[cond1]) + delta[cond1]/gamma_mcp[cond1]

  if (any(!cond1)) dir_der_pos[!cond1]<- -res_var[!cond1]
  if (any(!cond1)) dir_der_neg[!cond1]<-  res_var[!cond1]
  
  list(dir_der_pos = dir_der_pos, dir_der_neg = dir_der_neg)
}




sign0p <- function(x)
{
  out <- sign(x)
  out[x == 0] <- 1
  out
}

sign0m <- function(x)
{
  out <- -sign(x)
  out[x == 0]  <- 1
  out
}






# this function is used for finding the second derivatives of the profile log-likelihood
# with respect to the equating coefficients and the item parameters
# x = item parameters
derAB_lik_unc <- function(x, par, base, num.forms, aj1T, bj1T, var_aj1T, var_bj1T, Xmat) 
{
  x1 <- x[1:(length(x) / 2)]
  x2 <- x[(length(x) / 2 + 1):length(x)]
  aj1T[!is.na(aj1T)] <- x1
  bj1T[!is.na(bj1T)] <- x2
  derProfLik_unc_Rcpp(par = par, notbase = (1:num.forms)[-base] - 1, numforms = num.forms, aj1T = aj1T, bj1T = bj1T, var_aj1T = var_aj1T, var_bj1T = var_bj1T, Xmat = Xmat)
}

derAB_lik_1PL_unc <- function(x, par, base, num.forms, bj1T, var_bj1T, Xmat) 
{
  bj1T[!is.na(bj1T)] <- x
  derProfLik_unc_1PL_Rcpp(par = par, notbase = (1:num.forms)[-base] - 1, numforms = num.forms, bj1T = bj1T, var_bj1T = var_bj1T, Xmat = Xmat)
}



