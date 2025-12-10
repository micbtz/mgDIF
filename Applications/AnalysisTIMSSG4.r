

load("DataApplication.RData")
# objects in the file:
# data_840: list of dataframes including data from USA, one for each booklet
# data_b1E: list of dataframes, one from each European country, booklet 1
# X: item position in all the booklets
# Xsc: item position in all the booklets (accounting for science)
# X2: item position in all the booklets (linear and quadratic effects)

library(mirt)
library(equateMultiple)
source("../functions_IndLik.r")
Rcpp::sourceCpp("../IndLik.cpp")


# ===============================================================================
# Application 1: Position effect, data from USA
# =============================================================================== 

# item parameter estimation separately for each booklet
mods_840 <- list()
for (i in 1:length(data_840))
  mods_840[[i]] <- mirt(data_840[[i]][, -(1:3)], 1, itemtype = "2PL", SE = TRUE, 
                        technical = list(NCYCLES = 2000))

# extract item parameter estimates and covariance matrices
names <- sapply(data_840, function(x) paste("T", x$IDCNTRY[1], x$IDBOOK[1], sep = "_"))
mod2pl_840 <- modIRT(est.mods = mods_840, names = names)

# common items between booklets
linkp(mod2pl_840)

# estimation of equating coefficients with method multiple mean-geometric mean and multiple mean-mean
eqMGM_840 <- multiec(mod2pl_840, base = 1, method = "mean-gmean")
summary(eqMGM_840)

eqMM_840 <- multiec(mod2pl_840, base = 1, method = "mean-mean")
summary(eqMM_840)

# estimation of equating coefficients with the likelihood-based method assuming independence
eqLIK_840 <- multiec_lik_unc(mods = mod2pl_840, base = 1, se = TRUE, start = NULL, X = NULL, iter.max = 100000, trace = TRUE, check.all = FALSE) 
summary(eqLIK_840)

### inclusion of a linear position effect

# accounting for science items in computing the position
eqLIK_possc_840 <- multiec_lik_unc(mods = mod2pl_840, base = 1, se = TRUE, start = NULL, X = Xsc, iter.max = 100000, trace = TRUE, check.all = FALSE)
summary(eqLIK_possc_840)

# not accounting for science items in computing the position
eqLIK_pos_840 <- multiec_lik_unc(mods = mod2pl_840, base = 1, se = TRUE, start = NULL, X = X, iter.max = 100000, trace = TRUE, check.all = FALSE)
summary(eqLIK_pos_840)

# comparison of log-likelihood
eqLIK_pos_840$loglik # not accounting for science (chosen)
eqLIK_possc_840$loglik # accounting for science

# coefficients and standard error of the linear position effect
eqLIK_pos_840$coef
eqLIK_pos_840$se.coef
eqLIK_pos_840$coef/eqLIK_pos_840$se.coef

# likelihood ratio test: linear position effect against any position effect
eqLIK_840$loglik
eqLIK_pos_840$loglik
pchisq(2*(eqLIK_pos_840$loglik-eqLIK_840$loglik), 2, lower.tail = FALSE)

### inclusion of a quadratic position effect

eqLIK_pos2_840 <- multiec_lik_unc(mods = mod2pl_840, base = 1, se = TRUE, start = NULL, X = X2, iter.max = 100000, trace = TRUE, check.all = FALSE)
summary(eqLIK_pos2_840)

# comparison of log-likelihood
eqLIK_840$loglik # without position effect
eqLIK_pos_840$loglik # linear position effect
eqLIK_pos2_840$loglik # quadratic position effect

# coefficients and standard error of the quadratic position effect
eqLIK_pos2_840$coef
eqLIK_pos2_840$se.coef
eqLIK_pos2_840$coef/eqLIK_pos2_840$se.coef

# likelihood ratio test: quadratic position effect against any position effect
pchisq(2*(eqLIK_pos2_840$loglik-eqLIK_840$loglik), 4, lower.tail = FALSE)
# likelihood ratio test: quadratic position effect against linear position effect
pchisq(2*(eqLIK_pos2_840$loglik-eqLIK_pos_840$loglik), 2, lower.tail = FALSE)

# Figure of the quadratic position effect
pos <- 1:30
posEff <- data.frame(pos = pos, eff = c(eqLIK_pos2_840$coef[1] * pos + eqLIK_pos2_840$coef[2] * pos^2, eqLIK_pos2_840$coef[3] * pos + eqLIK_pos2_840$coef[4] * pos^2), type = rep(c("discrimination parameters", "difficulty parameters"), each = 30))
library(ggplot2)
dev.new(width = 8, height = 4, noRStudioGD = TRUE)
ggplot(posEff, 
       aes(x = pos, y = eff)) +
  geom_point()+
  facet_wrap(~ type) + 
  xlab("item position")+
  ylab("")+
  theme_bw()
savePlot("poseff", type = "pdf")


### delta parameters

lambda <- seq(1, 50, by = 0.5)

# without position effects
out_840 <- list()
# parameter estimation for each lambda
for (i in 1:length(lambda))
{
  eqLIK_lambda <- multiec_lik_unc(mods = mod2pl_840, base = 1, se = FALSE, start = eqLIK_840, X = NULL, iter.max = 100000, trace = TRUE, check.all = TRUE, lambda_mcp = lambda[i], tau_mcp = 3)
  out_840[[i]] <- eqLIK_lambda[c("A", "B", "coef", "loglik", "delta_a", "delta_b", "npar", "AIC", "BIC")]
}
bic_840 <- lapply(out_840, function(x) x$BIC)
plot(lambda, bic_840)
sel_bic_840 <- which.min(bic_840)
sel_bic_840
lambda[sel_bic_840]
out_840[sel_bic_840]

sela_840 <- which(out_840[[sel_bic_840]]$delta_a!= 0, arr.ind = TRUE)
selb_840 <- which(out_840[[sel_bic_840]]$delta_b!= 0, arr.ind = TRUE)
sela_840 # non-zero delta^a estimates
selb_840 # non-zero delta^b estimates


# with quadratic position effect
out1_840 <- list()
for (i in 1:length(lambda))
{
  eqLIK_lambda <- multiec_lik_unc(mods = mod2pl_840, base = 1, se = FALSE, start = eqLIK_840, X = X2, iter.max = 100000, trace = TRUE, check.all = TRUE, lambda = lambda[i], tau_mcp = 3)
  out1_840[[i]] <- eqLIK_lambda[c("A", "B", "coef", "loglik", "delta_a", "delta_b", "npar", "AIC", "BIC")]
}
bic1_840 <- sapply(out1_840, function(x) x$BIC)
plot(lambda, bic1_840)
sel_bic1_840 <- which.min(bic1_840)
sel_bic1_840
lambda[sel_bic1_840]
out1_840[sel_bic1_840]

sela1_840 <- which(out1_840[[sel_bic1_840]]$delta_a!= 0, arr.ind = TRUE)
selb1_840 <- which(out1_840[[sel_bic1_840]]$delta_b!= 0, arr.ind = TRUE)
sela1_840 # non-zero delta^a estimates
selb1_840 # non-zero delta^b estimates


# with linear position effect

out2_840 <- list()
for (i in 1:length(lambda))
{
  eqLIK_lambda <- multiec_lik_unc(mods = mod2pl_840, base = 1, se = FALSE, start = eqLIK_840, X = X, iter.max = 100000, trace = TRUE, check.all = TRUE, lambda = lambda[i], tau_mcp = 3)
  out2_840[[i]] <- eqLIK_lambda[c("A", "B", "coef", "loglik", "delta_a", "delta_b", "npar", "AIC", "BIC")]
}
bic2_840 <- lapply(out2_840, function(x) x$BIC)
plot(lambda, bic2_840)
sel_bic2_840 <- which.min(bic2_840)
sel_bic2_840
lambda[sel_bic2_840]
out2_840[sel_bic2_840]

sela2_840 <- which(out2_840[[sel_bic2_840]]$delta_a!= 0, arr.ind = TRUE)
selb2_840 <- which(out2_840[[sel_bic2_840]]$delta_b!= 0, arr.ind = TRUE)
sela2_840 # non-zero delta^a estimates
selb2_840 # non-zero delta^b estimates






# ====================================================================
# Application 2: European countries, booklet 1
# ==================================================================== 

# item parameter estimation separately for each country
mods_b1E <- list()
for (i in 1:length(data_b1E))
  mods_b1E[[i]] <- mirt(data_b1E[[i]][, -(1:3)], 1, itemtype = "2PL", SE = TRUE, 
                        technical = list(NCYCLES = 2000))

# extract item parameter estimates and covariance matrices
mod2pl_b1E <- modIRT(est.mods = mods_b1E, names = names(data_b1E))

# commom items between countries (booklet 1)
linkp(mod2pl_b1E)

# estimation of equating coefficients with method multiple mean-geometric mean and multiple mean-mean
eqMGM_b1E <- multiec(mod2pl_b1E, base = 1, method = "mean-gmean")
summary(eqMGM_b1E)

eqMM_b1E <- multiec(mod2pl_b1E, base = 1, method = "mean-mean")
summary(eqMM_b1E)

# estimation of equating coefficients with the likelihood-based method assuming independence
eqLIK_b1E <- multiec_lik_unc(mods = mod2pl_b1E, base = 1, se = TRUE, start = NULL, X = NULL, iter.max = 100000, trace = TRUE, check.all = FALSE) 
summary(eqLIK_b1E)


### delta parameters

lambda <- seq(1, 40, by = 0.5)
out_b1E <- list()
for (i in 1:length(lambda))
{
  eqLIK_lambda <- multiec_lik_unc(mods = mod2pl_b1E, base = 1, se = FALSE, start = eqLIK_b1E, X = NULL, iter.max = 100000, trace = TRUE, check.all = TRUE, lambda_mcp = lambda[i],tau_mcp = 3)
  out_b1E[[i]] <- eqLIK_lambda[c("A", "B", "coef", "loglik", "delta_a", "delta_b", "npar", "AIC", "BIC")]
}
bic_b1E <- lapply(out_b1E, function(x) x$BIC)
plot(lambda, bic_b1E)
sel_bic_b1E <- which.min(bic_b1E)
sel_bic_b1E
lambda[sel_bic_b1E]
abline(v = lambda[sel_bic_b1E])

# estimate of the parameters with the selected lambda
eqLIK_lambdas_b1E <- multiec_lik_unc(mods = mod2pl_b1E, base = 1, se = FALSE, start = NULL, X = NULL, iter.max = 100000, trace = TRUE, check.all = TRUE, lambda_mcp = lambda[sel_bic_b1E], tau_mcp = 3)
summary(eqLIK_lambdas_b1E)
eqLIK_lambdas_b1E$delta_a
eqLIK_lambdas_b1E$delta_b
sela <- which(eqLIK_lambdas_b1E$delta_a!= 0, arr.ind = TRUE)
selb <- which(eqLIK_lambdas_b1E$delta_b!= 0, arr.ind = TRUE)
sela # non-zero delta^a estimates
selb # non-zero delta^b estimates

write.csv2(eqLIK_lambdas_b1E$delta_b,file = "deltab_countries.csv")
write.csv2(eqLIK_lambdas_b1E$delta_a,file = "deltaa_countries.csv")


### Figures

# comparison of equating coefficients
EQappl2 <- data.frame(unadjusted = c(eqLIK_b1E$A, eqLIK_b1E$B), adjusted = c(eqLIK_lambdas_b1E$A, eqLIK_lambdas_b1E$B), EQ = rep(c("A", "B"), each = 18))
library(ggplot2)
rn <- range(EQappl2[, 1:2])
dev.new(noRStudioGD = TRUE, width = 6, height = 4)
ggplot(EQappl2, 
       aes(x = unadjusted, y = adjusted, shape = EQ)) +
  geom_point()+
  xlim(rn) +
  ylim(rn) +
  geom_abline(intercept = 0, slope = 1)+
  labs(shape = "equating coefficients") + 
  theme_bw()
savePlot("comp", type = "pdf")


# delta parameter estimates for all lambda
deltab <- sapply(out_b1E, function(x) x$delta_b)
deltaa <- sapply(out_b1E, function(x) x$delta_a)

# prepare results for figures
deltab <- as.data.frame(deltab)
deltab$items <- substr(rownames(out_b1E[[1]]$delta_b), 8, 20)
deltab$groups <- rep(colnames(out_b1E[[1]]$delta_b), each = nrow(out_b1E[[1]]$delta_b))
deltab$groups[deltab$groups == "SlovakRepublic"] <- "Slovak Republic"
deltab$groups[deltab$groups == "CzechRepublic"] <- "Czech Republic"
deltab_long <- reshape(deltab, direction = "long", v.names = "delta", varying = list(colnames(deltab)[1:(ncol(deltab)-2)]), times = lambda, timevar = "lambda")
deltab_long[1:10, ]

deltaa <- as.data.frame(deltaa)
deltaa$items <- substr(rownames(out_b1E[[1]]$delta_a), 8, 20)
deltaa$groups <- rep(colnames(out_b1E[[1]]$delta_a), each = nrow(out_b1E[[1]]$delta_a))
deltaa$groups[deltaa$groups == "SlovakRepublic"] <- "Slovak Republic"
deltaa$groups[deltaa$groups == "CzechRepublic"] <- "Czech Republic"
deltaa_long <- reshape(deltaa, direction = "long", v.names = "delta", varying = list(colnames(deltaa)[1:(ncol(deltaa)-2)]), times = lambda, timevar = "lambda")
deltaa_long[1:10, ]


library(ggplot2)

mycol <- c("cyan", "firebrick3", "springgreen3", "red", "orange", "green", "brown", 
           "mediumslateblue", "magenta", "purple", "turquoise3", "darkolivegreen1", 
           "dimgray", "royalblue", "thistle3", "salmon", "deepskyblue", "pink")


dev.new(noRStudioGD = TRUE, width = 8, height = 8)
ggplot(deltab_long, 
       aes(x = lambda, y = delta, col = groups)) +
  geom_line(lwd = 0.4)+
  facet_wrap(~items, ncol = 4)+
  geom_vline(xintercept = lambda[sel_bic_b1E], lty = 2)+
  xlab(expression(lambda))+
  ylab(expression(delta))+
  labs(colour = "Countries")+
  guides(color = guide_legend(override.aes = list(lwd = 1)))+
  scale_colour_manual(values = mycol)+
  theme_bw()+
  theme(strip.text = element_text(size = 9))

savePlot("regPathAllDiff", type = "pdf")

dev.new()
ggplot(deltab_long, 
       aes(x = lambda, y = delta, col = items)) +
  geom_line()+
  facet_wrap(~groups)+
  geom_vline(xintercept = lambda[sel_bic_b1E])+
  xlab(expression(lambda))+
  ylab(expression(delta))+
  theme_bw()



dev.new(noRStudioGD = TRUE, width = 8, height = 8)
ggplot(deltaa_long, 
       aes(x = lambda, y = delta, col = groups)) +
  geom_line(lwd = 0.4)+
  facet_wrap(~items, ncol = 4)+
  geom_vline(xintercept = lambda[sel_bic_b1E], lty = 2)+
  xlab(expression(lambda))+
  ylab(expression(delta))+
  labs(colour = "Countries")+
  guides(color = guide_legend(override.aes = list(lwd = 1)))+
  scale_colour_manual(values = mycol)+
  theme_bw()+
  theme(strip.text = element_text(size = 9))

savePlot("regPathAllDisc", type = "pdf")



dev.new()
ggplot(deltaa_long, 
       aes(x = lambda, y = delta, col = items)) +
  geom_line()+
  facet_wrap(~groups)+
  geom_vline(xintercept = lambda[sel_bic_b1E])+
  xlab(expression(lambda))+
  ylab(expression(delta))+
  theme_bw()




tmp1 <- deltab_long[deltab_long$items == "ME71041" | deltab_long$items == "ME71162", ]
tmp2 <- deltaa_long[deltaa_long$items == "ME71041" | deltaa_long$items == "ME71162", ]
tmp1$type <- "Difficulty"
tmp2$type <- "Discrimination"
delta_ab <- rbind(tmp1, tmp2)

mycol <- c("cyan", "firebrick3", "springgreen3", "red", "orange", "green", "brown", 
           "mediumslateblue", "magenta", "purple", "turquoise3", "darkolivegreen1", 
           "dimgray", "royalblue", "thistle3", "salmon", "deepskyblue", "pink")

dev.new(noRStudioGD = TRUE, width = 8.5, height = 6)

ggplot(delta_ab, 
       aes(x = lambda, y = delta, col = groups)) +
  geom_line(lwd = 0.4)+
  facet_grid(type~items)+
  geom_vline(xintercept = lambda[sel_bic_b1E], lty = 2)+
  xlab(expression(lambda))+
  ylab(expression(delta))+
  xlim(0,13)+
  labs(colour = "Countries")+
  # theme(legend.position = "bottom")+
  guides(color = guide_legend(override.aes = list(lwd = 1)))+
  scale_colour_manual(values = mycol)+
  theme_bw()

savePlot("regPath", type = "pdf")


library(equateIRT)
difts <- dif.test(est.mods = mods_b1E, purification = TRUE, method = "Haebara", signif.level = 0.0001)
difts

difts <- dif.test(est.mods = mods_b1E, purification = FALSE, method = "Haebara", signif.level = 0.0001)
difts





# ====================================================================
# Application 3: European countries and gender, booklet 1
# ==================================================================== 

# load data with gender
load("asg.RData")
asg<-do.call("rbind",asg)
rownames(asg)<-paste(asg$IDCNTRY,asg$IDSTUD,sep="_")

# split data on the basis of gender
data_b1E_gender<-list()

for (i in (1:length(data_b1E)))
{
  ids<-paste(data_b1E[[i]]$IDCNTRY,data_b1E[[i]]$IDSTUD,sep="_")
  gender<-asg[ids,]$ASBG01 # 1 girl, 2 boy
  tmp<-split(data_b1E[[i]],gender)
  data_b1E_gender[[(i-1)*2+1]]<-tmp[[1]]
  data_b1E_gender[[(i-1)*2+2]]<-tmp[[2]]
  names(data_b1E_gender)[[(i-1)*2+1]]<-paste(names(data_b1E)[[i]],"F",sep=" ")
  names(data_b1E_gender)[[(i-1)*2+2]]<-paste(names(data_b1E)[[i]],"M",sep=" ")
}
names(data_b1E_gender)
lapply(data_b1E_gender,dim)


#analysis of DIF for gender only

dataF<-do.call("rbind",data_b1E_gender[1:18*2-1])
dataM<-do.call("rbind",data_b1E_gender[1:18*2])

mods_gender <- list()
library(mirt)
mods_gender[[1]] <- mirt(dataF[, -(1:3)], 1, itemtype = "2PL", SE = TRUE, 
                         technical = list(NCYCLES = 20000))
mods_gender[[2]] <- mirt(dataM[, -(1:3)], 1, itemtype = "2PL", SE = TRUE, 
                         technical = list(NCYCLES = 20000))


# extract item parameter estimates and covariance matrices
mod2pl_gender <- modIRT(est.mods = mods_gender, names = c("F","M"))

# estimation of equating coefficients with the likelihood-based method assuming independence
eqLIK_gender <- multiec_lik_unc(mods = mod2pl_gender, base = 1, se = TRUE, start = NULL, X = NULL, iter.max = 100000, trace = TRUE, check.all = FALSE) 
summary(eqLIK_gender)


lambda <- seq(1, 20, by = 0.5)
out_gender <- list()
for (i in 1:length(lambda))
{
  eqLIK_lambda <- multiec_lik_unc(mods = mod2pl_gender, base = 1, se = FALSE, start = eqLIK_gender, X = NULL, iter.max = 100000, trace = TRUE, check.all = TRUE, lambda_mcp = lambda[i],tau_mcp = 3)
  out_gender[[i]] <- eqLIK_lambda[c("A", "B", "coef", "loglik", "delta_a", "delta_b", "npar", "AIC", "BIC")]
}
bic_gender <- sapply(out_gender, function(x) x$BIC)
plot(lambda, bic_gender)
sel_bic_gender <- which.min(bic_gender)
sel_bic_gender
lambda[sel_bic_gender]
abline(v = lambda[sel_bic_gender])

# estimate of the parameters with the selected lambda
eqLIK_lambdas_gender <- multiec_lik_unc(mods = mod2pl_gender, base = 1, se = FALSE, start = NULL, X = NULL, iter.max = 100000, trace = TRUE, check.all = TRUE, lambda_mcp = lambda[sel_bic_gender], tau_mcp = 3)
summary(eqLIK_lambdas_gender)
eqLIK_lambdas_gender$delta_a
eqLIK_lambdas_gender$delta_b
sela_gender <- which(eqLIK_lambdas_gender$delta_a!= 0, arr.ind = TRUE)
selb_gender <- which(eqLIK_lambdas_gender$delta_b!= 0, arr.ind = TRUE)
sela_gender # non-zero delta^a estimates
selb_gender # non-zero delta^b estimates

# delta parameter estimates for all lambda
deltab_b1E_gender <- sapply(out_gender, function(x) x$delta_b)
deltaa_b1E_gender <- sapply(out_gender, function(x) x$delta_a)

# prepare results for figures
deltab_b1E_gender <- as.data.frame(deltab_b1E_gender)
deltab_b1E_gender$items <- substr(rownames(out_gender[[1]]$delta_b), 8, 20)
deltab_b1E_gender$groups <- rep(colnames(out_gender[[1]]$delta_b), each = nrow(out_gender[[1]]$delta_b))
deltab_b1E_gender_long <- reshape(deltab_b1E_gender, direction = "long", v.names = "delta", varying = list(colnames(deltab_b1E_gender)[1:(ncol(deltab_b1E_gender)-2)]), times = lambda, timevar = "lambda")
deltab_b1E_gender_long[1:10, ]

deltaa_b1E_gender <- as.data.frame(deltaa_b1E_gender)
deltaa_b1E_gender$items <- substr(rownames(out_gender[[1]]$delta_a), 8, 20)
deltaa_b1E_gender$groups <- rep(colnames(out_gender[[1]]$delta_a), each = nrow(out_gender[[1]]$delta_a))
deltaa_b1E_gender_long <- reshape(deltaa_b1E_gender, direction = "long", v.names = "delta", varying = list(colnames(deltaa_b1E_gender)[1:(ncol(deltaa_b1E_gender)-2)]), times = lambda, timevar = "lambda")
deltaa_b1E_gender_long[1:10, ]


library(ggplot2)

mycol <- c("cyan", "firebrick3", "springgreen3", "red", "orange", "green", "brown", 
           "mediumslateblue", "magenta", "purple", "turquoise3", "darkolivegreen1", 
           "dimgray", "royalblue", "thistle3", "salmon", "deepskyblue", "pink")
# mycol<-rep(mycol,each=2)

dev.new(noRStudioGD = TRUE, width = 8, height = 8)
ggplot(deltab_b1E_gender_long, 
       aes(x = lambda, y = delta, col = groups)) +
  geom_line(lwd = 0.7)+
  facet_wrap(~items, ncol = 4)+
  geom_vline(xintercept = lambda[sel_bic_gender], lty = 2)+
  xlab(expression(lambda))+
  ylab(expression(delta))+
  labs(colour = "Countries")+
  guides(color = guide_legend(override.aes = list(lwd = 1)))+
  scale_colour_manual(values = mycol)+
  theme_bw()+
  theme(strip.text = element_text(size = 9))

# savePlot("regPathAllDiff", type = "pdf")

dev.new()
ggplot(deltab_b1E_gender_long, 
       aes(x = lambda, y = delta, col = items)) +
  geom_line()+
  facet_wrap(~groups)+
  geom_vline(xintercept = lambda[sel_bic_gender])+
  xlab(expression(lambda))+
  ylab(expression(delta))+
  theme_bw()



dev.new(noRStudioGD = TRUE, width = 8, height = 8)
ggplot(deltaa_b1E_gender_long, 
       aes(x = lambda, y = delta, col = groups)) +
  geom_line(lwd = 0.7)+
  facet_wrap(~items, ncol = 4)+
  geom_vline(xintercept = lambda[sel_bic_b1E], lty = 2)+
  xlab(expression(lambda))+
  ylab(expression(delta))+
  labs(colour = "Countries")+
  guides(color = guide_legend(override.aes = list(lwd = 1)))+
  scale_colour_manual(values = mycol)+
  theme_bw()+
  theme(strip.text = element_text(size = 9))

# savePlot("regPathAllDisc", type = "pdf")



dev.new()
ggplot(deltaa_b1E_gender_long, 
       aes(x = lambda, y = delta, col = items)) +
  geom_line()+
  facet_wrap(~groups)+
  geom_vline(xintercept = lambda[sel_bic_b1E])+
  xlab(expression(lambda))+
  ylab(expression(delta))+
  theme_bw()




library(equateIRT)
difts_gender <- dif.test(est.mods = mods_gender, purification = TRUE, method = "Haebara", signif.level = 0.05)
difts_gender






















