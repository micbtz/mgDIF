

# ========================================================================================
# DIF


load("eqLIK_b1E.RData")

hist(eqLIK_b1E$as)
min(eqLIK_b1E$as) # 0.5627793
max(eqLIK_b1E$as) # 1.897122
hist(eqLIK_b1E$bs)
mean(eqLIK_b1E$bs) # -0.2063936
sd(eqLIK_b1E$bs) # 0.9927318



set.seed(1)
nitems <- 30
discrm <- runif(nitems, 0.5, 2)
diff <- rnorm(nitems, 0, 1)
itempar <- cbind(discrm, diff)
rownames(itempar) <- paste0("I", formatC(1:nitems, flag = "0", width = 2))
plot(discrm, diff)

ngroups <- 10

mean_gr <- eqLIK_b1E$B[1:10]
sd_gr <- eqLIK_b1E$A[1:10]
Atrue <- eqLIK_b1E$A[2:10]
Btrue <- eqLIK_b1E$B[2:10]

dif_effect <- matrix(NA, 3, 2)
dif_effect[1, ] <- c(0.10, 0.2)
dif_effect[2, ] <- c(0.25, 0.5)
dif_effect[3, ] <- c(0.5, 1)

dif_effect_unif <- matrix(NA, 3, 2)
dif_effect_unif[1, ] <- c(0, 0.2)
dif_effect_unif[2, ] <- c(0, 0.5)
dif_effect_unif[3, ] <- c(0, 1)


settings <- expand.grid(samplsize = c(250, 500, 1000, 2000), percdif = c(0.01, 0.05, 0.1), dif_level = 1:3, dif_type = c("u", "nu"))
settings
dim(settings)


source("FzSim.r")
library(parallel)
no_cores <- detectCores() - 1

tau<-3

R <- 500 # number of replications

# res <- list()
for (s in 1:nrow(settings))
{
  print(s)
  samplsize <- settings$samplsize[s]
  percdif <- settings$percdif[s]
  dif_level <- settings$dif_level[s]
  dif_type <- settings$dif_type[s]
  
  cl <- makeCluster(no_cores)
  res_s <- parLapply(cl, 1:R, sim_multidif, itempar = itempar, ngroups = ngroups, 
                     nitems = nitems, samplsize = samplsize, percdif = percdif, 
                     dif_level = dif_level, dif_type = dif_type, 
                     dif_effect_unif = dif_effect_unif, dif_effect = dif_effect, 
                     mean_gr = mean_gr, sd_gr = sd_gr, tau_mcp=tau)
  stopCluster(cl)
  
  res[[s]] <- res_s
}

S<-72

# equating coefficients estimated ignoring DIF
estA <- lapply(res, function(x) sapply(x, function(x) x$eqLIK1$A))
estB <- lapply(res, function(x) sapply(x, function(x) x$eqLIK1$B))
# bias of equating coefficients estimated ignoring DIF
biasA <- lapply(estA, function(x) rowMeans(x[-1, ] - Atrue))
biasB <- lapply(estB, function(x) rowMeans(x[-1, ] - Btrue))
# average bias for each setting
aveBiasA <- sapply(biasA, mean)
aveBiasB <- sapply(biasB, mean)
round(aveBiasA, 3)
round(aveBiasB, 3)

# DIF detection

# equating coefficients estimated accounting for DIF
# lambda selected using AIC
estAadj_aic <- lapply(res, function(x) sapply(x, function(x) x$out_aic$A))
estBadj_aic <- lapply(res, function(x) sapply(x, function(x) x$out_aic$B))
biasAadj_aic <- lapply(estAadj_aic, function(x) rowMeans(x[-1, ]-Atrue))
biasBadj_aic <- lapply(estBadj_aic, function(x) rowMeans(x[-1, ]-Btrue))
aveBiasAadj_aic <- sapply(biasAadj_aic, mean)
aveBiasBadj_aic <- sapply(biasBadj_aic, mean)
aveBiasAadj_aic
aveBiasBadj_aic
plot(aveBiasA, aveBiasAadj_aic)
abline(0, 1)
plot(aveBiasB, aveBiasBadj_aic)
abline(0, 1)

# lambda selected using BIC
estAadj_bic <- lapply(res, function(x) sapply(x, function(x) x$out_bic$A))
estBadj_bic <- lapply(res, function(x) sapply(x, function(x) x$out_bic$B))
biasAadj_bic <- lapply(estAadj_bic, function(x) rowMeans(x[-1, ]-Atrue))
biasBadj_bic <- lapply(estBadj_bic, function(x) rowMeans(x[-1, ]-Btrue))
aveBiasAadj_bic <- sapply(biasAadj_bic, mean)
aveBiasBadj_bic <- sapply(biasBadj_bic, mean)
aveBiasAadj_bic
aveBiasBadj_bic
dev.new()
plot(aveBiasA, aveBiasAadj_bic)
abline(0, 1)
plot(aveBiasB, aveBiasBadj_bic)
abline(0, 1)

# tuning parameters selected
lambda_aic <- lapply(res, function(x) sapply(x, function(x) seq(1, 100, by = 0.5)[x$sel_aic]))
lambda_bic <- lapply(res, function(x) sapply(x, function(x) seq(1, 100, by = 0.5)[x$sel_bic]))
lambda_aic
lambda_bic
sort(sapply(lambda_bic, max))
sort(sapply(lambda_bic, min))

# matrices indicating items with DIF
difitems_mat <- lapply(res, function(x) lapply(x, function(x) x$difitems_mat))


# number of parameters

npar <- lapply(res, function(x) sapply(x, function(x) x$out_bic$npar))
npar

# DIF detection
any_dif <- lapply(difitems_mat, function(x) lapply(x, rowSums)) # for each setting and for each replication, for each item, the number of groups with DIF
lapply(any_dif, function(x) sapply(x, table))
lapply(any_dif, function(x) max(sapply(x, max))) # for each setting, maximum number of groups with DIF for each item across all replications

lapply(any_dif, function(x) table(do.call("c", x))/(R*nitems*ngroups))

round(dbinom(0:10, 10, 0.01), 3)
round(dbinom(0:10, 10, 0.05), 3)
round(dbinom(0:10, 10, 0.1), 3)

1-dbinom(0, 10, 0.01)
1-dbinom(0, 10, 0.05)
1-dbinom(0, 10, 0.1)



# detection of DIF in difficulty parameters

detect <- list()
for (s in 1:S)
{
  print(s)
  detect_s <- list()
  for (r in 1:R)
  {
    nas <- rep(NA, 30)
    mat_s_r <- data.frame(nDIF = nas, nDIFdet = nas, nNOTDIFdet = nas, meanDIFdet = nas, meanNOTDIFdet = nas, rev = nas)
    delta <- res[[s]][[r]]$out_bic$delta_b
    dif <- difitems_mat[[s]][[r]]
    delta_dif <- delta
    delta_dif[!dif] <- NA # delta of dif items
    mat_s_r$nDIF <- rowSums(dif) # number of groups with DIF for each item
    mat_s_r$nDIFdet <- rowSums(delta_dif!=0, na.rm = T) # number of groups detected for each item with DIF (true positives)
    delta_notdif <- delta
    delta_notdif[dif] <- NA # delta of non dif items
    mat_s_r$nNOTDIFdet <- rowSums(delta_notdif!=0, na.rm = T) # number of groups detected for each item without DIF (false positives)
    delta_dif[delta_dif==0] <- NA
    delta_notdif[delta_notdif==0] <- NA
    mat_s_r$meanDIFdet <- rowMeans(delta_dif, na.rm = T) # delta mean of true positives
    mat_s_r$meanNOTDIFdet <- rowMeans(delta_notdif, na.rm = T) # delta mean of false positives
    
    # condition: number of DIF groups>=5, number of false positive > number of true positives, delta mean of false positive > delta mean of true positives
    cond <- mat_s_r$nDIF>=5 & mat_s_r$nNOTDIFdet>mat_s_r$nDIFdet & abs(mat_s_r$meanNOTDIFdet)>abs(mat_s_r$meanDIFdet)
    if (any(cond, na.rm = TRUE))
    {
      print(c(s,r))
      mat_s_r[cond, ] <- mat_s_r[cond, c(1, 3, 2, 5, 4)] # reverse of false positive and true positives
      mat_s_r[cond, c(4,5)] <- - mat_s_r[cond, c(4,5)] # reverse the sign of mean delta
      mat_s_r[cond, ]$rev <- 1
    }
    detect_s[[r]] <- mat_s_r
  }
  detect[[s]] <- detect_s
}



res_detect <- matrix(NA, 72, 6)
for (s in 1:S)
{
  print(s)
  detect_s <- do.call("rbind", detect[[s]])
  res_detect[s, ] <- colSums(detect_s, na.rm = TRUE)
  res_detect[s, 4:5] <- colMeans(abs(detect_s[, 4:5]), na.rm = TRUE)
}
res_detect

res_detect1 <- res_detect
res_detect1[, 2] <- res_detect1[, 2] / res_detect1[, 1] # rate of cases detected out of the DIF cases (true positives)
res_detect1[, 3] <- res_detect1[, 3] / (30 * 10 * R - res_detect1[, 1]) # rate of cases detected out of the non-DIF cases (false positives)

res_detect1 <- as.data.frame(res_detect1)
colnames(res_detect1) <- c("nDIF", "rateTruePositives", "rateFalsePositives", "MeanDeltaTP", "MeanDeltaFP")
res_detect1 <- cbind(settings, res_detect1)

res_detect1$percdif <- as.factor(res_detect1$percdif)
levels(res_detect1$percdif) <- c("1%", "5%", "10%")
res_detect1$dif_level <- as.factor(res_detect1$dif_level)
levels(res_detect1$dif_level) <- c("low", "medium", "high")
levels(res_detect1$dif_type) <- c("uniform", "not uniform")

# rateTruePositives: number of true positives out of the number of items*groups with DIF (computed for each item and each group)
# rateFalsePositives: number of false positives out of the number of items*groups without DIF (computed for each item and each group)
# MeanDeltaTP: Mean of delta of true positives (only deltas>0)
# MeanDeltaFP: Mean of delta of false positives (only deltas>0) 


# detection of DIF in discrimination parameters

detect_a <- list()
for (s in 1:S)
{
  print(s)
  detect_s <- list()
  for (r in 1:R)
  {
    nas <- rep(NA, 30)
    mat_s_r <- data.frame(nDIF = nas, nDIFdet = nas, nNOTDIFdet = nas, meanDIFdet = nas, meanNOTDIFdet = nas, rev = nas)
    delta <- res[[s]][[r]]$out_bic$delta_a
    if (settings$dif_type[s]=="u") dif <- matrix(FALSE, 30, 10)
    else dif <- difitems_mat[[s]][[r]]
    delta_dif <- delta
    delta_dif[!dif] <- NA # delta of dif items
    mat_s_r$nDIF <- rowSums(dif) # number of groups with DIF for each item
    mat_s_r$nDIFdet <- rowSums(delta_dif!=0, na.rm = T) # number of groups detected for each item with DIF (true positives)
    delta_notdif <- delta
    delta_notdif[dif] <- NA # delta of non dif items
    mat_s_r$nNOTDIFdet <- rowSums(delta_notdif!=0, na.rm = T) # number of groups detected for each item without DIF (false positives)
    delta_dif[delta_dif==0] <- NA
    delta_notdif[delta_notdif==0] <- NA
    mat_s_r$meanDIFdet <- rowMeans(delta_dif, na.rm = T) # delta mean of true positives
    mat_s_r$meanNOTDIFdet <- rowMeans(delta_notdif, na.rm = T) # delta mean of false positives
    
    # condizione: number of DIF groups>=5, number of false positive > number of true positives, delta mean of false positive > delta mean of true positives
    cond <- mat_s_r$nDIF>=5 & mat_s_r$nNOTDIFdet>mat_s_r$nDIFdet & abs(mat_s_r$meanNOTDIFdet)>abs(mat_s_r$meanDIFdet)
    if (any(cond, na.rm = TRUE))
    {
      mat_s_r[cond, ] <- mat_s_r[cond, c(1, 3, 2, 5, 4)] # reverse of false positive and true positives
      mat_s_r[cond, c(4,5)] <- - mat_s_r[cond, c(4,5)] # reverse the sign of mean delta
      mat_s_r[cond, ]$rev <- 1
    }
    detect_s[[r]] <- mat_s_r
  }
  detect_a[[s]] <- detect_s
}



res_detect_a <- matrix(NA, 72, 6)
for (s in 1:S)
{
  print(s)
  detect_s <- do.call("rbind", detect_a[[s]])
  res_detect_a[s, ] <- colSums(detect_s, na.rm = TRUE)
  res_detect_a[s, 4:5] <- colMeans(abs(detect_s[, 4:5]), na.rm = TRUE)
}
res_detect_a

res_detect1_a <- res_detect_a
res_detect1_a[, 2] <- res_detect1_a[, 2]/res_detect1_a[, 1] # rate of cases detected out of the DIF cases
res_detect1_a[, 3] <- res_detect1_a[, 3]/(30*10*R-res_detect1_a[, 1]) # rate of cases detected out of the non-DIF cases

res_detect1_a <- as.data.frame(res_detect1_a)
colnames(res_detect1_a) <- c("nDIF", "rateTruePositives", "rateFalsePositives", "MeanDeltaTP", "MeanDeltaFP")
res_detect1_a <- cbind(settings, res_detect1_a)

res_detect1_a$percdif <- as.factor(res_detect1_a$percdif)
levels(res_detect1_a$percdif) <- c("1%", "5%", "10%")
res_detect1_a$dif_level <- as.factor(res_detect1_a$dif_level)
levels(res_detect1_a$dif_level) <- c("low", "medium", "high")
levels(res_detect1_a$dif_type) <- c("uniform", "not uniform")



# computation for each item (independently of which group, diff or discr)
# before for each item*group separately for diff and discr

res_detect_ab <- matrix(NA, 72, 2)
for (s in 1:S)
{
  print(s)
  detect_b_s <- do.call("rbind", detect[[s]])
  detect_a_s <- do.call("rbind", detect_a[[s]])
  isdif <- detect_b_s$nDIF>0 | detect_a_s$nDIF>0
  isdet <- (detect_a_s$nDIFdet + detect_a_s$nNOTDIFdet + detect_b_s$nDIFdet + detect_b_s$nNOTDIFdet)>0
  res_detect_ab[s, 1] <- sum(isdif & isdet) / sum(isdif)
  res_detect_ab[s, 2] <- sum(!isdif & isdet) / sum(!isdif)
}
res_detect_ab

res_detect_ab <- as.data.frame(res_detect_ab)
colnames(res_detect_ab) <- c("tpr", "fpr")
# tpr: number of dif items detected as dif in at least one group in either difficulty or discrimination parameters out of the number of dif items
# fpr: number of non dif items detected as dif in at least one group in either difficulty or discrimination parameters out of the number of non dif items

res_detect_ab <- cbind(settings, res_detect_ab)

res_detect_ab$percdif <- as.factor(res_detect_ab$percdif)
levels(res_detect_ab$percdif) <- c("1%", "5%", "10%")
res_detect_ab$dif_level <- as.factor(res_detect_ab$dif_level)
levels(res_detect_ab$dif_level) <- c("low", "medium", "high")
levels(res_detect_ab$dif_type) <- c("uniform", "not uniform")




library(ggplot2)

# difficulty and discrimination, all groups together

dev.new(noRStudioGD = TRUE, width = 8.5, height = 3.8)
ggplot(res_detect_ab, 
       aes(x = samplsize, y = tpr, col = percdif, 
           shape = dif_level)) + 
  geom_point() + 
  geom_line() + 
  facet_wrap(~dif_type, ncol = 4) + 
  xlab("sample size") + 
  ylab("true positive rate") + 
  labs(colour = "DIF percentage") + 
  labs(shape = "DIF level") + 
  ylim(0, 1) + 
  guides(color = guide_legend(override.aes = list(shape="", lwd = 1))) + 
  theme_bw()

savePlot(paste("tpr_all_tau",tau,sep=""), type="pdf")



ggplot(res_detect_ab, 
       aes(x = samplsize, y = fpr, col = percdif, 
           shape = dif_level)) + 
  geom_point() + 
  geom_line() + 
  facet_wrap(~dif_type, ncol = 4) + 
  xlab("sample size") + 
  ylab("false positive rate") + 
  labs(colour = "DIF percentage") + 
  labs(shape = "DIF level") + 
  ylim(0, 1) + 
  guides(color = guide_legend(override.aes = list(shape="", lwd = 1))) + 
  theme_bw()

savePlot(paste("fpr_all_tau",tau,sep=""), type="pdf")



res_detect1$type <- "Difficulty"
res_detect1_a$type <- "Discrimination"
res_detect1 <- rbind(res_detect1, res_detect1_a[37:72, ])

# difficulty parameters

dev.new(noRStudioGD = TRUE, width = 11, height = 3.8)
ggplot(res_detect1, 
       aes(x = samplsize, y = rateTruePositives, col = percdif, 
           shape = dif_level)) + 
  geom_point() + 
  geom_line() + 
  facet_wrap(~dif_type + type) + 
  xlab("sample size") + 
  ylab("true positive rate") + 
  labs(colour = "DIF percentage") + 
  labs(shape = "DIF level") + 
  ylim(0, 1) + 
  guides(color = guide_legend(override.aes = list(shape="", lwd = 1))) + 
  theme_bw()

savePlot(paste("tpr_tau",tau,sep=""), type="pdf")


ggplot(res_detect1, 
       aes(x = samplsize, y = rateFalsePositives, col = percdif, 
           shape = dif_level)) + 
  geom_point() + 
  geom_line() + 
  facet_wrap(~dif_type + type) + 
  xlab("sample size") + 
  ylab("false positive rate") + 
  labs(colour = "DIF percentage") + 
  labs(shape = "DIF level") + 
  #ylim(0, 1) + 
  guides(color = guide_legend(override.aes = list(shape="", lwd = 1))) + 
  theme_bw()

savePlot(paste("fpr_tau",tau,sep=""), type="pdf")



ggplot(res_detect1, 
       aes(x = samplsize, y = MeanDeltaTP, col = percdif, 
           shape = dif_level)) + 
  geom_point() + 
  geom_line() + 
  facet_wrap(~dif_type + type) + 
  xlab("sample size") + 
  ylab("") + 
  ylim(0, 1) + 
  labs(colour = "DIF percentage") + 
  labs(shape = "DIF level") + 
  guides(color = guide_legend(override.aes = list(shape="", lwd = 1))) + 
  theme_bw()

savePlot(paste("deltatp_tau",tau,sep=""), type="pdf")


ggplot(res_detect1, 
       aes(x = samplsize, y = MeanDeltaFP, col = percdif, 
           shape = dif_level)) + 
  geom_point() + 
  geom_line() + 
  facet_wrap(~dif_type + type, ncol = 4) + 
  xlab("sample size") + 
  ylab("") + 
  ylim(0, 1) + 
  labs(colour = "DIF percentage") + 
  labs(shape = "DIF level") + 
  guides(color = guide_legend(override.aes = list(shape="", lwd = 1))) + 
  theme_bw()

savePlot(paste("deltafp_tau",tau,sep=""), type="pdf")
















# some Wald tests are NA
NAflag <- lapply(res, function(x) sapply(x, function(x) any(is.na(x$res_diftest))))
lapply(NAflag, table)


# number of items with DIF in at least one group (considered only if wald test works)
num_difitems <- lapply(res, function(x) sum(sapply(x, numDif), na.rm = T))
# number of items without DIF in at least one group (considered only if wald test works)
num_nondifitems <- lapply(res, function(x) sum(sapply(x, numNotDif), na.rm = T))
# mapply(" + ", num_difitems, num_nondifitems)
# number of true positives according to the Wald test
num_true_pos <- lapply(res, function(x) sum(sapply(x, tp), na.rm = T))
# number of false positives according to the Wald test
num_false_pos <- lapply(res, function(x) sum(sapply(x, fp), na.rm = T))
true_pos_rate <- mapply(FUN="/", num_true_pos, num_difitems)
false_pos_rate <- mapply(FUN="/", num_false_pos, num_nondifitems)



resWald <- cbind(settings, true_pos_rate)

resWald$percdif <- as.factor(resWald$percdif)
levels(resWald$percdif) <- c("1%", "5%", "10%")
resWald$dif_level <- as.factor(resWald$dif_level)
levels(resWald$dif_level) <- c("low", "medium", "high")
levels(resWald$dif_type) <- c("uniform", "not uniform")



library(ggplot2)

dev.new(noRStudioGD = TRUE, width = 8.5, height = 3.8)

ggplot(resWald, 
       aes(x = samplsize, y = true_pos_rate, col = percdif, 
           shape = dif_level)) + 
  geom_point() + 
  geom_line() + 
  facet_wrap(~dif_type) + 
  xlab("sample size") + 
  ylab("true positive rate") + 
  labs(colour = "prop DIF") + 
  labs(shape = "DIF level") + 
  ylim(0, 1) + 
  guides(color = guide_legend(override.aes = list(shape="", lwd = 1))) + 
  theme_bw()
savePlot("tpr_wald", type="pdf")

ggplot(resWald, 
       aes(x = samplsize, y = false_pos_rate, col = percdif, 
           shape = dif_level)) + 
  geom_point() + 
  geom_line() + 
  facet_wrap( ~ dif_type) + 
  xlab("sample size") + 
  ylab("false positive rate") + 
  labs(colour = "prop DIF") + 
  labs(shape = "DIF level") + 
  ylim(0, 1) + 
  guides(color = guide_legend(override.aes = list(shape="", lwd = 1))) + 
  theme_bw()





# ========================================================================================
# no DIF

# load("C:/Users/Battauz/OneDrive - Università degli Studi di Udine/PCufficio/Documenti/RICERCA/2402 LikelihoodEquatingDIF/TIMSS19/eqLIK_b1E.RData")
load("../TIMSS19/eqLIK_b1E.RData")

set.seed(1)
nitems <- 30
discrm <- runif(nitems, 0.5, 2)
diff <- rnorm(nitems, 0, 1)
itempar <- cbind(discrm, diff)
rownames(itempar) <- paste0("I", formatC(1:nitems, flag = "0", width = 2))
plot(discrm, diff)

ngroups <- 10

mean_gr <- eqLIK_b1E$B[1:10]
sd_gr <- eqLIK_b1E$A[1:10]
Atrue <- eqLIK_b1E$A[2:10]
Btrue <- eqLIK_b1E$B[2:10]

settings <- expand.grid(samplsize = c(250, 500, 1000, 2000), percdif = 0)
settings


source("FzSim.r")
library(parallel)
no_cores <- detectCores() - 1

tau<-3

R <- 500 # number of replications

res <- list()
for (s in 1:nrow(settings))
{
  print(s)
  samplsize <- settings$samplsize[s]
  percdif <- settings$percdif[s]
  
  cl <- makeCluster(no_cores)
  res_s <- parLapply(cl, 1:R, sim_multidif, itempar = itempar, ngroups = ngroups, 
                     nitems = nitems, samplsize = samplsize, percdif = percdif, 
                     dif_level = NULL, dif_type = "none", 
                     dif_effect_unif = "none", dif_effect = "none", 
                     mean_gr = mean_gr, sd_gr = sd_gr, tau_mcp=tau)
  stopCluster(cl)
  
  res[[s]] <- res_s
}



# equating coefficients estimates
estA <- lapply(res, function(x) sapply(x, function(x) x$eqLIK1$A))
estB <- lapply(res, function(x) sapply(x, function(x) x$eqLIK1$B))
# bias of equating coefficients estimated ignoring DIF
biasA <- lapply(estA, function(x) rowMeans(x[-1, ] - Atrue))
biasB <- lapply(estB, function(x) rowMeans(x[-1, ] - Btrue))
# average bias for each setting
aveBiasA <- sapply(biasA, mean)
aveBiasB <- sapply(biasB, mean)
round(aveBiasA, 3)
round(aveBiasB, 3)

npar <- lapply(res, function(x) sapply(x, function(x) x$out_bic$npar))
npar

delta_b<-lapply(res, function(x) lapply(x, function(x) x$out_bic$delta_b))
delta_a<-lapply(res, function(x) lapply(x, function(x) x$out_bic$delta_a))

NumDetectedb<-lapply(delta_b, function(x) sapply(x, function(x) sum(x!=0)))
NumDetecteda<-lapply(delta_a, function(x) sapply(x, function(x) sum(x!=0)))

MeanDetectedb<-lapply(delta_b, function(x) sapply(x, function(x) mean(x[x!=0])))
MeanDetecteda<-lapply(delta_a, function(x) sapply(x, function(x) mean(x[x!=0])))

rateFPb<-sapply(NumDetectedb,sum)/(R*30*10)
rateFPa<-sapply(NumDetecteda,sum)/(R*30*10)

meanFPb<-sapply(MeanDetectedb,mean,na.rm=T)
meanFPa<-sapply(MeanDetecteda,mean,na.rm=T)

res_detect<-data.frame(samplsize=settings$samplsize,type=rep(c("Difficulty","Discrimination"),each=4),
                       rateFP=c(rateFPb,rateFPa),meanFP=c(meanFPb,meanFPa))



library(ggplot2)

dev.new(noRStudioGD = TRUE, width = 8.5, height = 3.8)

ggplot(res_detect, 
       aes(x = samplsize, y = rateFP)) + 
  geom_point() + 
  geom_line() + 
  facet_wrap(~ type) + 
  xlab("sample size") + 
  ylab("false positive rate") + 
  # ylim(0, 1) + 
  # guides(color = guide_legend(override.aes = list(shape="", lwd = 1))) + 
  theme_bw()

savePlot(paste("fpr_NODIF_tau",tau,sep=""), type="pdf")


ggplot(res_detect, 
       aes(x = samplsize, y = meanFP)) + 
  geom_point() + 
  geom_line() + 
  facet_wrap(~ type) + 
  xlab("sample size") + 
  ylab("") + 
  #ylim(0, 1) + 
  # guides(color = guide_legend(override.aes = list(shape="", lwd = 1))) + 
  theme_bw()

savePlot(paste("deltafp_noDIF_tau",tau,sep=""), type="pdf")









