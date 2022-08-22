#---- Fit SEIR model to San Juan serotype case data -----#
# libraries
rm(list=ls())
library(cmdstanr)
library(rstan)
library(data.table)
library(ggplot2)
library(DirichletReg)


# read in Puerto Rico data
setwd('/Users/wangyr/Desktop')
df <- read.csv('san_juan_testing_data.csv')
dates <- unique(df$week_start_date)
Ncases <- matrix(NA, nrow=length(dates), ncol=4)
for(s in 1:4) Ncases[,s] <- df[,3+s]


# inputs for model
data <- list()
data$Ncases <- Ncases
data$NcasesTot <- df$total_cases
data$Nt <- nrow(data$Ncases)
data$Nt2 <- data$Nt*7
data$pop <- c(rep(2327000,35),
              rep(2345000,52),
              rep(2363000,52),
              rep(2381000,52),
              rep(2400000,52),
              rep(2418000,52),
              rep(2437000,52),
              rep(2456000,52),
              rep(2475000,52),
              rep(2494000,52),
              rep(2508000,52),
              rep(2505000,52),
              rep(2502000,52),
              rep(2499000,52),
              rep(2496000,52),
              rep(2493000,52),
              rep(2490000,52),
              rep(2487000,52),
              rep(2484000,52),
              rep(2481000,52),
              rep(2478000,52),
              rep(2475000,52),
              rep(2472000,52),
              rep(2469000,17))
data$sigma <- 1/4
data$omega <- 1/6
data$kappa <- 1/6
data$gamma <- 1/80/365
data$K <- 26
data$phi <-1
data$zeta <- 1/2/365
data$alpha <- c(20,60,rep(5,8),rep(0.1,16))/2
data$delta <- 0.3
indL <- seq(1,data$Nt*7,7)
indU <- seq(7,data$Nt*7+6,7)
data$indL <- indL
data$indU <- indU
data$ind <- sort(rep(seq(1,data$Nt),7))
data$pSample <- rowSums(data$Ncases)/data$NcasesTot
data$pSample[is.na(data$pSample)] <- 0


# chain starting values
alpha <- c(20,60,rep(5,8),rep(0.1,16))/1
initI <- rdirichlet(100, alpha)
ii <- as.vector(initI[1,])
inits_chain1 <- list(log_rho=-1.6, init=ii, b=1, log_var=-3, log_phiNB=2, log_B=matrix(0.3,ncol=4,nrow=data$Nt))


# fit model
check_cmdstan_toolchain(fix=T)
set_cmdstan_path()
mod <- cmdstan_model('scenarios_CD.stan', pedantic=T)
fit <- mod$sample(data=data, chains=1, parallel_chains=1, iter_sampling=1500, refresh=10, iter_warmup=500, save_warmup=TRUE, init=list(inits_chain1))
stanfit <- rstan::read_stan_csv(fit$output_files())

# Check convergence
chains <- rstan::extract(stanfit)
traceplot1 <- traceplot(stanfit, pars=c('init','lp__'), inc_warmup=T)
traceplot2 <- traceplot(stanfit, pars=c('rho','b','log_var','log_phiNB'),inc_warmup=T)

# Initial condition estimates
pars <- c('S','R','MoP1','MoP2','MoP3','MoP4','MoS1','MoS2','MoS3','MoS4','E11','E12','E13','E14','I11','I12','I13','I14',
          'E21','E22','E23','E24','I21','I22','I23','I24')
init <- data.frame(pars=pars,med=NA,ciL=NA,ciU=NA)
for(i in 1:nrow(init)){
  init[i,2:4] <- quantile(chains$init[,i], c(0.5,0.025,0.975))
}
init$sero <- c(0,0,rep(seq(1,4),6))
initplot <- ggplot(init, aes(pars,med, col=factor(sero)))+ geom_point()+
  theme_minimal()+ geom_linerange(aes(ymin=ciL,ymax=ciU))+ xlab('Pars')+ 
  ylab('Proportion')


# plot time-series fit
fitS <- data.frame(date=as.Date(dates),t=rep(seq(1:data$Nt),4),sero=sort(rep(seq(1,4),data$Nt)),obs=NA, med=NA,ciL=NA,ciU=NA)
fitTot <- data.frame(date=as.Date(dates),t=seq(1,data$Nt),obs=data$NcasesTot, med=NA,ciL=NA,ciU=NA)
for(i in 1:data$Nt){
  for(s in 1:4){
    fitS[fitS$t==i & fitS$sero==s,5:7] <- quantile(chains$predCases[,i,s], c(0.5,0.025,0.975))
    fitS[fitS$t==i & fitS$sero==s,4] <- data$Ncases[i,s]
  }
  fitTot[i,4:6] <- quantile(chains$predCasesTot[,i], c(0.5,0.025,0.975))
}

# Time-series plots
fitPlotS <- ggplot(fitS, aes(date,obs))+ geom_point(alpha=0.4, col='grey40')+
  theme_minimal()+ theme(legend.position='none')+ theme(text=element_text(size=16))+
  geom_line(aes(date,med,col=factor(sero)))+ geom_ribbon(aes(date, ymin=ciL,ymax=ciU,fill=factor(sero)),alpha=0.3)+
  xlab('Time')+ ylab('Cases')+ theme(legend.position='none')+ facet_wrap(~sero,scale='free_y')
fitPlotS
fitPlotTot <- ggplot(fitTot, aes(date,obs))+ geom_point(alpha=0.4, col='grey40')+
  theme_minimal()+ geom_line(aes(date,med), col='blue')+ xlab('Time')+
  geom_ribbon(aes(ymin=ciL,ymax=ciU), fill='dodgerblue', alpha=0.5)+ ylab('Cases')+
  theme(text=element_text(size=16))
fitPlotTot

# Epidemiological parameter estimates
BtS <- data.frame(date=as.Date(dates),t=rep(seq(1:data$Nt),4),sero=sort(rep(seq(1,4),data$Nt)),med=NA,ciL=NA,ciU=NA)
Rt <- data.frame(date=as.Date(dates),t=rep(seq(1:data$Nt),4),sero=sort(rep(seq(1,4),data$Nt)),med=NA,ciL=NA,ciU=NA)
lam <- data.frame(date=seq(1:data$Nt2),t=rep(seq(1:data$Nt2),4),sero=sort(rep(seq(1,4),data$Nt)),med=NA,ciL=NA,ciU=NA)
for(i in 1:data$Nt){
  for(s in 1:4){
    BtS[BtS$t==i & BtS$sero==s,4:6] <- quantile(chains$B[,i,s], c(0.5,0.025,0.975))
    Rt[Rt$t==i & Rt$sero==s,4:6] <- quantile(chains$Rt[,i,s], c(0.5,0.025,0.975))
    
  }
}
for(i in 1:data$Nt2) for(s in 1:4){
  lam[lam$t==i & lam$sero==s,4:6] <- quantile(chains$lam[,i,s], c(0.5,0.025,0.975))
}
BtPlot <- ggplot(BtS, aes(date,med,col=factor(sero)))+ geom_point()+
  facet_wrap(~sero, scales='free_y')+ theme_minimal()+ theme(text=element_text(size=16))+
  geom_linerange(aes(ymin=ciL,ymax=ciU,fill=factor(sero)), alpha=0.3)+
  theme(legend.position='none')+ ylab('Bt')+ xlab('Time')
BtPlot
RtPlot <- ggplot(Rt, aes(date,med,col=factor(sero)))+ geom_point()+
  facet_wrap(~sero)+ theme_minimal()+ ylim(0,NA)+ theme(text=element_text(size=16))+
  geom_linerange(aes(ymin=ciL,ymax=ciU,fill=factor(sero)), alpha=0.3)+
  theme(legend.position='none')+ ylab('Rt')+ xlab('Time')+
  geom_hline(yintercept=1, linetype='dashed')
RtPlot
lamPlot <- ggplot(lam, aes(date,med))+ geom_line(aes(col=factor(sero)))+
  facet_wrap(~sero, scales='free_y')+ theme_minimal()+ ylim(0,NA)+ theme(text=element_text(size=16))+
  geom_ribbon(aes(ymin=ciL,ymax=ciU,fill=factor(sero)), alpha=0.3)+
  theme(legend.position='none')+ ylab('FOI')+ xlab('days')
lamPlot

setwd('/Users/wangyr/Desktop/result')
saveRDS(data,'Data.RDS')
saveRDS(chains, 'Chains.RDS')
png(filename='initplot.png', width=20, height=14, res=400, units='cm')
plot(initplot)
dev.off()
png(filename='FitTotalCases.png', width=20, height=14, res=400, units='cm')
plot(fitPlotTot)
dev.off()
png(filename='FitSeroCases.png', width=20, height=14, res=400, units='cm')
plot(fitPlotS)
dev.off()
png(filename='RtPlot.png', width=20, height=15, res=400, units='cm')
plot(RtPlot)
dev.off()
png(filename='BtPlot.png', width=20, height=15, res=400, units='cm')
plot(BtPlot)
dev.off()
png(filename='FOIPlot.png', width=20, height=15, res=400, units='cm')
plot(lamPlot)
dev.off()
png(filename='traceplot1.png', width=30, height=14, res=400, units='cm')
plot(traceplot1)
dev.off()
png(filename='traceplot2.png', width=20, height=14, res=400, units='cm')
plot(traceplot2)
dev.off()

