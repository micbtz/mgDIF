





load("eqLIK_b1E.RData")

ngroups<-10

set.seed(1)
nitems<-200
discrm<-runif(nitems,0.5,2)
diff<-rnorm(nitems,0,1)
plot(discrm,diff)
itempar<-cbind(discrm,diff)
rownames(itempar)<-paste0("I", formatC(1:nitems,flag = "0", width = 3))

itempar_groups <- list()
for (g in 1:(ngroups-1))
{
  itempar_groups[[g]] <- itempar[(1:40) + (20 * (g - 1)), ]
}
itempar_groups[[ngroups]]<-itempar[c(181:200,1:20),]




mean_gr<-eqLIK_b1E$B[1:10]
sd_gr<-eqLIK_b1E$A[1:10]
Atrue<-eqLIK_b1E$A[2:10]
Btrue<-eqLIK_b1E$B[2:10]


settings<-expand.grid(samplsize=c(250,500,1000,2000),coef=c(0.01,0.02))
settings
dim(settings)



source("FzSim.r")
library(parallel)
no_cores <- detectCores() - 1

R <- 500 # number of replications

res<-list()
for (s in 1:nrow(settings))
{
  print(s)
  samplsize <- settings$samplsize[s]
  cf <- settings$coef[s]
  
  # tmp<-sim_multidif_pos(r=r, itempar_groups = itempar_groups, ngroups=ngroups,samplsize=samplsize, mean_gr = mean_gr, sd_gr = sd_gr, cf=cf)
  
  cl <- makeCluster(no_cores)
  res_s <- parLapply(cl, 1:R, sim_multidif_pos, itempar_groups=itempar_groups, ngroups=ngroups, 
                     samplsize=samplsize, mean_gr = mean_gr, sd_gr = sd_gr, cf=cf)
  stopCluster(cl)
  
  res[[s]]<-res_s
}



# equating coefficients estimated ignoring DIF
estA<-lapply(res,function(x) sapply(x, function(x) x$eqLIK1$A))
estB<-lapply(res,function(x) sapply(x, function(x) x$eqLIK1$B))
# bias of equating coefficients estimated ignoring DIF
biasA<-lapply(estA,function(x) rowMeans(x[-1,]-Atrue))
biasB<-lapply(estB,function(x) rowMeans(x[-1,]-Btrue))
# average bias for each setting
aveBiasA<-sapply(biasA,function(x) mean(abs(x)))
aveBiasB<-sapply(biasB,function(x) mean(abs(x)))
round(aveBiasA,3)
round(aveBiasB,3)
# mse of equating coefficients accounting for the position
rmseA<-lapply(estA,function(x) sqrt(rowMeans((x[-1,]-Atrue)^2)))
rmseB<-lapply(estB,function(x) sqrt(rowMeans((x[-1,]-Btrue)^2)))
aveRmseA<-sapply(rmseA,mean)
aveRmseB<-sapply(rmseB,mean)
sdA<-lapply(estA,sd)
sdB<-lapply(estB,sd)
aveSDA<-sapply(sdA,mean)
aveSDB<-sapply(sdB,mean)

# equating coefficients accounting for the position
estAc<-lapply(res,function(x) sapply(x, function(x) x$eqLIK2$A))
estBc<-lapply(res,function(x) sapply(x, function(x) x$eqLIK2$B))
# bias of equating coefficients accounting for the position
biasAc<-lapply(estAc,function(x) rowMeans(x[-1,]-Atrue))
biasBc<-lapply(estBc,function(x) rowMeans(x[-1,]-Btrue))
# average bias for each setting
aveBiasAc<-sapply(biasAc,function(x) mean(abs(x)))
aveBiasBc<-sapply(biasBc,function(x) mean(abs(x)))
round(aveBiasAc,3)
round(aveBiasBc,3)
# rmse of equating coefficients accounting for the position
rmseAc<-lapply(estAc,function(x) sqrt(rowMeans((x[-1,]-Atrue)^2)))
rmseBc<-lapply(estBc,function(x) sqrt(rowMeans((x[-1,]-Btrue)^2)))
aveRmseAc<-sapply(rmseAc,mean)
aveRmseBc<-sapply(rmseBc,mean)
sdAc<-lapply(estAc,sd)
sdBc<-lapply(estBc,sd)
aveSDAc<-sapply(sdAc,mean)
aveSDBc<-sapply(sdBc,mean)



bias_rmse<-cbind(settings,aveBiasA,aveBiasAc,aveRmseA,aveRmseAc,aveBiasB,aveBiasBc,aveRmseB,aveRmseBc)
bias_rmse_long<-reshape(bias_rmse,direction="long",varying = list(colnames(bias_rmse)[c(3,5,7,9)],colnames(bias_rmse)[c(4,6,8,10)]),
        times = c("aveBIASA","aveRMSEA","aveBIASB","aveRMSEB"))
colnames(bias_rmse_long)[4]<-"bias_rmse"
colnames(bias_rmse_long)[5]<-"bias_rmse_corr"
bias_rmse_long$eqc<-rep(c("A","B"),each=16)
bias_rmse_long$coef<-as.character(bias_rmse_long$coef)
bias_rmse_long$time<-substr(bias_rmse_long$time,4,7)





rn<-c(0,max(bias_rmse_long$bias_rmse))

library(ggplot2)

dev.new(noRStudioGD = TRUE, width=8, height=6)

ggplot(bias_rmse_long, 
       aes(x = bias_rmse, y = bias_rmse_corr, col = samplsize, 
           shape = coef)) + 
  geom_point() + 
  geom_abline(intercept = 0, slope = 1)+
  facet_grid(time~eqc) + 
  xlim(rn) +
  ylim(rn) +
  xlab("unadjusted") + 
  ylab("adjusted") + 
  labs(colour = "sample size") + 
  labs(shape = "position effect") + 
  theme_bw()
savePlot("BiasRmse",type="pdf")


# nel calcolo della posizione, togliere quelli eliminati


coefs<-lapply(res,function(x) sapply(x, function(x) x$eqLIK2$coef))
sapply(coefs,rowMeans)

round(sapply(coefs[1:4],rowMeans)-0.01,5)
round(sapply(coefs[5:8],rowMeans)-0.02,5)





