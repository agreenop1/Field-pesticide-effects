setwd("C:/Users/arrgre/iCloudDrive/PHD WORK/Field pest control/Paper_1/Analysis")
#packages
library (dplyr)
library(rstan)
library(rstanarm)
library(brms)
library(ggplot2)
library(bayesplot)
library(broom)
library(broom)
library(boot)

#ALL DATA
ALL <-read.csv("ALL_FEEDING.csv",header=T)

#####RAW FEEDINF DATA FOR EACH SPECIES
Alive<-subset(ALL,ALL$ID_2=="ALIVE")

#SUBSET EACH SPECIES
A_pleb <- subset (Alive, Alive$species_name=="Amara_plebja")
write.csv(data.frame(A_pleb$mean_mass),"A_pleb_m.csv")

H_affi <- subset (Alive, Alive$species_name=="Harpalus_affinis") 
write.csv(data.frame(H_affi$mean_mass),"H_affi_m.csv")

H_axyr <- subset (Alive, Alive$species_name=="Harmonia_axyridis")
write.csv(data.frame(H_axyr$mean_mass),"H_axyr_m.csv")

N_brev <- subset (Alive, Alive$species_name=="Nebria_brevicollis") 
write.csv(data.frame(N_brev$mean_mass),"N_brev_m.csv")

B_bull <- subset (Alive, Alive$species_name=="Badister_bullatus")
write.csv(data.frame(B_bull$mean_mass),"B_bull_m.csv")

P_cogn <- subset (Alive, Alive$species_name=="Philonthus_cognatus")
write.csv(data.frame(P_cogn$mean_mass),"P_cogn_m.csv")

P_cupr <- subset (Alive, Alive$species_name=="Poecilus _cupreus")
write.csv(data.frame(P_cupr$mean_mass),"P_cupr_m.csv")

A_dors <- subset (Alive, Alive$species_name=="Anchomenus_dorsalis")
write.csv(data.frame(A_dors$mean_mass),"A_dors_m.csv")

P_madi <- subset (Alive, Alive$species_name=="Pterostichus_madidus")
write.csv(data.frame(P_madi$mean_mass),"P_madi_m.csv")

P_mela <- subset (Alive, Alive$species_name=="Pterostichus_melanarius")
write.csv(data.frame(P_mela$mean_mass),"P_mela_m.csv")

H_rufi <- subset (Alive, Alive$species_name=="Harpalus_rufipes")
write.csv(data.frame(H_rufi$mean_mass),"H_rufi_m.csv")

C_sept <- subset (Alive, Alive$species_name=="Cocinellidae_septempunctata")
write.csv(data.frame(C_sept$mean_mass),"C_sept_m.csv")

#CUSTOM CODE FOR BETA BINOMIAL DISTRIBUTION BASED ON BRMS TUTORIAL
beta_binomial2 <- custom_family(
  "beta_binomial2", dpars = c("mu", "phi"),
  links = c("logit", "log"), lb = c(NA, 0),
  type = "int", vars = "vint1[n]"
)
stan_funs <- "
  real beta_binomial2_lpmf(int y, real mu, real phi, int T) {
return beta_binomial_lpmf(y | T, mu * phi, (1 - mu) * phi);
}
int beta_binomial2_rng(real mu, real phi, int T) {
return beta_binomial_rng(T, mu * phi, (1 - mu) * phi);
}
"
stanvars <- stanvar(scode = stan_funs, block = "functions")

Pred <- data.frame(Species=0,Response=0,Mean=0,lb=0,ub=0)
#MODEL OUT FUNCTION
mod_out <-function (x){
  a1<- posterior_samples(x)
  if(is.null(a1$b_TreatmentRecovery)){
    #only resistance treatment
    out<-cbind(round((inv.logit (a1$b_Intercept)*20),2),
             round((inv.logit (a1$b_Intercept+a1$b_TreatmentResistance)*20),2))}else{
    #both treatments
    out<-cbind(round((inv.logit (a1$b_Intercept)*20),2),
               round((inv.logit (a1$b_Intercept+a1$b_TreatmentResistance)*20),2),
               round((inv.logit (a1$b_Intercept+a1$b_TreatmentRecovery)*20),2))           
             }
  return(out)
}
####################################################################
#########Harpalus affinis###########################################
####################################################################
#prior 10
H_aff10 <- brm( Prey_E | vint(Prey_D) ~ 1 + Treatment+(1|Block)
              ,data = H_affi,stanvars = stanvars,
              family = beta_binomial2, 
              prior = c(set_prior ("normal(0,10)",class="b"),
                        set_prior ("normal(0,10)",class="Intercept")),
              chains = 4,
              iter = 4000,control = list(adapt_delta =0.99 ,max_treedepth=15) ,cores=4,
              warmup = 1000)

#prior 5
H_aff2 <- brm( Prey_E | vint(Prey_D) ~ 1 + Treatment+(1|Block)
               ,data = H_affi,stanvars = stanvars,
               family = beta_binomial2, 
               prior = c(set_prior ("normal(0,5)",class="b"),
                         set_prior ("normal(0,5)",class="Intercept")),
               chains = 4,
               iter = 4000,control = list(adapt_delta =0.99 ,max_treedepth=15) ,cores=4,
               warmup = 1000)

#prior 2.5
H_aff1.4 <- brm( Prey_E | vint(Prey_D) ~ 1 + Treatment+(1|Block)
               ,data = H_affi,stanvars = stanvars,
               family = beta_binomial2, 
               prior = c(set_prior ("normal(0,2.5)",class="b"),
                         set_prior ("normal(0,2.5)",class="Intercept")),
               chains = 4,
               iter = 4000,control = list(adapt_delta =0.99 ,max_treedepth=15) ,cores=4,
               warmup = 1000)

#prior 1
H_aff1 <- brm( Prey_E | vint(Prey_D) ~ 1 + Treatment+(1|Block)
                 ,data = H_affi,stanvars = stanvars,
                 family = beta_binomial2, 
                 prior = c(set_prior ("normal(0,1)",class="b"),
                           set_prior ("normal(0,1)",class="Intercept")),
                 chains = 4,
                 iter = 4000,control = list(adapt_delta =0.99 ,max_treedepth=15) ,cores=4,
                 warmup = 1000)


#posterior distribution from the priors on the count scale 
H_aff_p <-data.frame(ABase=0,Aexp_1=0,Aexp_5=0,
                     BBase=0,Bexp_1=0,Bexp_5=0,
                     CBase=0,Cexp_1=0,Cexp_5=0,
                     DBase=0,Dexp_1=0,Dexp_5=0)


#DIFFERENT PRIOR MODEL OUTPUTS
H_aff_p[c(1:12000),1:3]<- mod_out(H_aff1.4)

H_aff_p[c(1:12000),4:6]<- mod_out(H_aff2)

H_aff_p[c(1:12000),7:9]<- mod_out(H_aff10)

H_aff_p[c(1:12000),10:12]<- mod_out(H_aff1)
###
write.csv (H_aff_p,"H_affR.csv")

#THESE FUNCTIONS ALLOW ESTIMATES AND PRIOR PREDICTIVE CHECKS TO BE RECOVERED FROM BETA BINOMIAL
expose_functions(H_aff1, vectorize = TRUE)
log_lik_beta_binomial2 <- function(i, draws) {
  mu <- draws$dpars$mu[, i]
  phi <- draws$dpars$phi
  trials <- draws$data$vint1[i]
  y <- draws$data$Y[i]
  beta_binomial2_lpmf(y, mu, phi, trials)
}
predict_beta_binomial2 <- function(i, draws, ...) {
  mu <- draws$dpars$mu[, i]
  phi <- draws$dpars$phi
  trials <- draws$data$vint1[i]
  beta_binomial2_rng(mu, phi, trials)
}
fitted_beta_binomial2 <- function(draws) {
  mu <- draws$dpars$mu
  trials <- draws$data$vint1
  trials <- matrix(trials, nrow = nrow(mu), ncol = ncol(mu), byrow = TRUE)
  mu * trials
}

####################################################################
#########Harmonia axyridis##########################################
####################################################################
H_axyr$orle <- paste0(seq(1:nrow(H_axyr)),"A")
##
#prior 10
H_ax10 <- brm(Prey_E|trials(Prey_D)~ 1 + Treatment+(1|Block)+(1|orle)
             ,data = H_axyr,
             family = "binomial", 
             prior = c(set_prior ("normal(0,10)",class="b"),
                       set_prior ("normal(0,10)",class="Intercept")),
             chains = 4,cores=4,
             iter = 4000,control = list(adapt_delta =0.99 ,max_treedepth=15) ,
             warmup = 1000)
#prior 5
H_ax2 <- brm(Prey_E|trials(Prey_D)~ 1 + Treatment+(1|Block)+(1|orle)
            ,data = H_axyr,
            family = "binomial", 
            prior = c(set_prior ("normal(0,5)",class="b"),
                      set_prior ("normal(0,5)",class="Intercept")),
            chains = 4,
            iter = 4000,control = list(adapt_delta =0.99 ,max_treedepth=15) ,
            warmup = 1000)
#prior 2.5
H_ax1.4 <- brm(Prey_E|trials(Prey_D)~ 1 + Treatment+(1|Block)+(1|orle)
              ,data = H_axyr,
              family = "binomial", 
              prior = c(set_prior ("normal(0,2.5)",class="b"),
                        set_prior ("normal(0,2.5)",class="Intercept")),
              chains = 4,
              iter = 4000,control = list(adapt_delta =0.99 ,max_treedepth=15) ,
              warmup = 1000)
#prior 1
H_ax1 <- brm(Prey_E|trials(Prey_D)~ 1 + Treatment+(1|Block)+(1|orle)
               ,data = H_axyr,
               family = "binomial", 
               prior = c(set_prior ("normal(0,1)",class="b"),
                         set_prior ("normal(0,1)",class="Intercept")),
               chains = 4,
               iter = 4000,control = list(adapt_delta =0.99 ,max_treedepth=15) ,
               warmup = 1000)
##posterior distribution from the priors on the count scale 
H_ax_p <-data.frame(ABase=0,Aexp_1=0,Aexp_5=0,
                     BBase=0,Bexp_1=0,Bexp_5=0,
                     CBase=0,Cexp_1=0,Cexp_5=0,
                     DBase=0,Dexp_1=0,Dexp_5=0)

#DIFFERENT PRIOR MODEL OUTPUTS
H_ax_p[c(1:12000),1:3]<-   mod_out(H_ax1.4)
H_ax_p[c(1:12000),4:6]<-   mod_out(H_ax2)
H_ax_p[c(1:12000),7:9]<-   mod_out(H_ax10)
H_ax_p[c(1:12000),10:12]<- mod_out(H_ax1)


write.csv(H_ax_p,"H_axyR.csv")
###################################################################
#####################Nebira brevicollis############################
###################################################################
N_brev$orle <- paste0(seq(1:nrow(N_brev)),"A")
#prior 5
n_br2 <- brm(Prey_E|trials(Prey_D)~ 1 + Treatment+(1|Block)+(1|orle)
              ,data = N_brev,
              family = "binomial", 
              prior = c(set_prior ("normal(0,5)",class="b"),
                        set_prior ("normal(0,5)",class="Intercept")),
              chains = 4,cores=4,
              iter = 4000,control = list(adapt_delta =0.99 ,max_treedepth=15) ,
              warmup = 1000)
#prior 10
n_br10 <- brm(Prey_E|trials(Prey_D)~ 1 + Treatment+(1|Block)+(1|orle)
               ,data = N_brev,
               family = binomial (link=logit), 
               prior = c(set_prior ("normal(0,10)",class="b"),
                         set_prior ("normal(0,10)",class="Intercept")),
               chains = 4,cores=4,
               iter = 4000,control = list(adapt_delta =0.99 ,max_treedepth=15) ,
               warmup = 1000)
#prior 2.5
n_br1.4 <- brm(Prey_E|trials(Prey_D)~ 1 + Treatment+(1|Block)+(1|orle)
                ,data = N_brev,
                family = "binomial", 
                prior = c(set_prior ("normal(0,2.5)",class="b"),
                          set_prior ("normal(0,2.5)",class="Intercept")),
                chains = 4,cores=4,
                iter = 4000,control = list(adapt_delta =0.99 ,max_treedepth=15) ,
                warmup = 1000)
#prior 1
n_br1 <- brm(Prey_E|trials(Prey_D)~ 1 + Treatment+(1|Block)+(1|orle)
               ,data = N_brev,
               family = "binomial", 
               prior = c(set_prior ("normal(0,1)",class="b"),
                         set_prior ("normal(0,1)",class="Intercept")),
               chains = 4,cores=4,
               iter = 4000,control = list(adapt_delta =0.99 ,max_treedepth=15) ,
               warmup = 1000)
##posterior distribution from the priors on the count scale 
n_br_p <-data.frame(ABase=0,Aexp_1=0,Aexp_5=0,
                    BBase=0,Bexp_1=0,Bexp_5=0,
                    CBase=0,Cexp_1=0,Cexp_5=0,
                    DBase=0,Dexp_1=0,Dexp_5=0)

#DIFFERENT PRIOR MODEL OUTPUTS
n_br_p[c(1:12000),1:3]<-   mod_out(n_br1.4)
n_br_p[c(1:12000),4:6]<-   mod_out(n_br2)
n_br_p[c(1:12000),7:9]<-   mod_out(n_br10)
n_br_p[c(1:12000),10:12]<- mod_out(n_br1)


write.csv(n_br_p,"N_brevR.csv")
########################################################################
#####################Badister bullatus##################################
########################################################################
B_bull$orle <- paste0(seq(1:nrow(B_bull)),"A")
B_bull<- subset (B_bull,B_bull$Treatment!="Recovery")
#prior 1
b_1b <- brm(Prey_E|trials(Prey_D)~ 1 + Treatment+(1|Block)
              ,data = B_bull,
              family = "binomial", 
              prior = c(set_prior ("normal(0,1)",class="b"),
                        set_prior ("normal(0,1)",class="Intercept")),
              chains = 4,
              iter = 4000,control = list(adapt_delta =0.99 ,max_treedepth=11) ,cores=4,
              warmup = 1000)
#prior 2.5
b_1.4b <- brm(Prey_E|trials(Prey_D)~ 1 + Treatment+(1|Block)
              ,data = B_bull,
              family = "binomial", 
              prior = c(set_prior ("normal(0,2.5)",class="b"),
                        set_prior ("normal(0,2.5)",class="Intercept")),
              chains = 4,
              iter = 4000,control = list(adapt_delta =0.99 ,max_treedepth=11) ,cores=4,
              warmup = 1000)
#prior 5
b_2b <- brm(Prey_E|trials(Prey_D)~ 1 + Treatment+(1|Block)
              ,data = B_bull,
              family = "binomial", 
              prior = c(set_prior ("normal(0,5)",class="b"),
                        set_prior ("normal(0,5)",class="Intercept")),
              chains = 4,
              iter = 4000,control = list(adapt_delta =0.99 ,max_treedepth=11) ,cores=4,
              warmup = 1000)
#prior 10
b_10b <- brm(Prey_E|trials(Prey_D)~ 1 + Treatment+(1|Block)
            ,data = B_bull,
            family = "binomial", 
            prior = c(set_prior ("normal(0,10)",class="b"),
                      set_prior ("normal(0,10)",class="Intercept")),
            chains = 4,
            iter = 4000,control = list(adapt_delta =0.99 ,max_treedepth=11) ,cores=4,
            warmup = 1000)
##posterior distribution from the priors on the count scale 
B_bul_p <-data.frame(ABase=0,Aexp_1=0,
                    BBase=0,Bexp_1=0,
                    CBase=0,Cexp_1=0,
                    DBase=0,Dexp_1=0)

#DIFFERENT PRIOR MODEL OUTPUTS
B_bul_p[c(1:12000),1:2]<-   mod_out(b_1.4b)
B_bul_p[c(1:12000),3:4]<-   mod_out(b_2b)
B_bul_p[c(1:12000),5:6]<-   mod_out(b_10b)
B_bul_p[c(1:12000),7:8]<-   mod_out(b_1b)


write.csv(B_bul_p,"B_bullR.csv")

#########################################################################
#####################Philonthus cognatus#################################
#########################################################################
P_cogn$orle<-paste0(seq(1:nrow(P_cogn)),"A")
#prior 1
P_con1<- brm(Prey_E|vint(Prey_D)~ 1 + Treatment+(1|Block)
               ,data = P_cogn,
               family = beta_binomial2, stanvars = stanvars,
               prior = c(set_prior ("normal(0,1)",class="b"),
                         set_prior ("normal(0,1)",class="Intercept")),
               chains = 4,cores=4,
               iter = 4000,control = list(adapt_delta =0.99 ,max_treedepth=15) ,
               warmup = 1000)
#prior 2.5
P_con1.4<- brm(Prey_E|vint(Prey_D)~ 1 + Treatment+(1|Block)
           ,data = P_cogn,
           family = beta_binomial2, stanvars = stanvars,
           prior = c(set_prior ("normal(0,2.5)",class="b"),
                     set_prior ("normal(0,2.5)",class="Intercept")),
           chains = 4,cores=4,
           iter = 4000,control = list(adapt_delta =0.99 ,max_treedepth=15) ,
           warmup = 1000)
#prior 5
P_con2<- brm(Prey_E|vint(Prey_D)~ 1 + Treatment+(1|Block)
             ,data = P_cogn,
             family = beta_binomial2, stanvars = stanvars,
             prior = c(set_prior ("normal(0,5)",class="b"),
                       set_prior ("normal(0,5)",class="Intercept")),
             chains = 4,cores=4,
             iter = 4000,control = list(adapt_delta =0.99 ,max_treedepth=15) ,
             warmup = 1000)
#prior 10
P_con10<- brm(Prey_E|vint(Prey_D)~ 1 + Treatment+(1|Block)
             ,data = P_cogn,
             family = beta_binomial2, stanvars = stanvars,
             prior = c(set_prior ("normal(0,10)",class="b"),
                       set_prior ("normal(0,10)",class="Intercept")),
             chains = 4,cores=4,
             iter = 4000,control = list(adapt_delta =0.99 ,max_treedepth=15) ,
             warmup = 1000)

##posterior distribution from the priors on the count scale 
n_Pc_p <-data.frame(ABase=0,Aexp_1=0,Aexp_5=0,
                    BBase=0,Bexp_1=0,Bexp_5=0,
                    CBase=0,Cexp_1=0,Cexp_5=0,
                    DBase=0,Dexp_1=0,Dexp_5=0)

#DIFFERENT PRIOR MODEL OUTPUTS
n_Pc_p[c(1:12000),1:3]<-   mod_out(P_con1.4)
n_Pc_p[c(1:12000),4:6]<-   mod_out(P_con2)
n_Pc_p[c(1:12000),7:9]<-   mod_out(P_con10)
n_Pc_p[c(1:12000),10:12]<- mod_out(P_con1)

write.csv(n_Pc_p,"P_cognR.csv")

####################################################################
#################Poecilus cupreus###################################
####################################################################
P_cupr$orle <-paste0(seq(1:nrow(P_cupr)),"A") 
#prior 1
P_cupo1<- brm(Prey_E|vint(Prey_D)~ 1 + Treatment+(1|Block)
                ,data = P_cupr,
                family = beta_binomial2,stanvars=stanvars, 
                prior = c(set_prior ("normal(0,1)",class="b"),
                          set_prior ("normal(0,1)",class="Intercept")),
                chains = 4,cores=4,
                iter = 4000,control = list(adapt_delta =0.99 ,max_treedepth=15) ,
                warmup = 1000)
#prior 2.5
P_cupo1.4<- brm(Prey_E|vint(Prey_D)~ 1 + Treatment+(1|Block)
             ,data = P_cupr,
             family = beta_binomial2,stanvars=stanvars, 
             prior = c(set_prior ("normal(0,2.5)",class="b"),
                       set_prior ("normal(0,2.5)",class="Intercept")),
             chains = 4,cores=4,
             iter = 4000,control = list(adapt_delta =0.99 ,max_treedepth=15) ,
             warmup = 1000)
#prior 5
P_cupo2<- brm(Prey_E|vint(Prey_D)~ 1 + Treatment+(1|Block)
             ,data = P_cupr,
             family = beta_binomial2,stanvars=stanvars, 
             prior = c(set_prior ("normal(0,5)",class="b"),
                       set_prior ("normal(0,5)",class="Intercept")),
             chains = 4,cores=4,
             iter = 4000,control = list(adapt_delta =0.99 ,max_treedepth=15) ,
             warmup = 1000)
#prior 10
P_cupo10<- brm(Prey_E|vint(Prey_D)~ 1 + Treatment+(1|Block)
             ,data = P_cupr,
             family = beta_binomial2,stanvars=stanvars, 
             prior = c(set_prior ("normal(0,10)",class="b"),
                       set_prior ("normal(0,10)",class="Intercept")),
             chains = 4,cores=4,
             iter = 4000,control = list(adapt_delta =0.99 ,max_treedepth=15) ,
             warmup = 1000)

##posterior distribution from the priors on the count scale 
n_Pcu_p <-data.frame(ABase=0,Aexp_1=0,Aexp_5=0,
                     BBase=0,Bexp_1=0,Bexp_5=0,
                     CBase=0,Cexp_1=0,Cexp_5=0,
                     DBase=0,Dexp_1=0,Dexp_5=0)

#DIFFERENT PRIOR MODEL OUTPUTS
n_Pcu_p[c(1:12000),1:3]<-   mod_out(P_cupo1.4)
n_Pcu_p[c(1:12000),4:6]<-   mod_out(P_cupo2)
n_Pcu_p[c(1:12000),7:9]<-   mod_out(P_cupo10)
n_Pcu_p[c(1:12000),10:12]<- mod_out(P_cupo1)

write.csv(n_Pcu_p,"P_cuprR.csv")
#####################Anchomenus dorsalis#########################
#####################Anchomenus dorsalis#########################
#####################Anchomenus dorsalis#########################
A_dors2 <- subset (A_dors,A_dors$Treatment!="Recovery")
A_dors2$orle<-paste0(seq(1:nrow(A_dors1)),"A") 
#prior 1
A_d1<- brm(Prey_E|vint(Prey_D)~ 1 + Treatment+(1|Block)
             ,data = A_dors2,
             family = beta_binomial2,  stanvars = stanvars,
             prior = c(set_prior ("normal(0,1)",class="b"),
                       set_prior ("normal(0,1)",class="Intercept")),
             chains = 4,cores=4,
             iter = 4000,control = list(adapt_delta =0.99 ,max_treedepth=15) ,
             warmup = 1000)
#prior 2.5
A_d1.4<- brm(Prey_E|vint(Prey_D)~ 1 + Treatment+(1|Block)
           ,data = A_dors2,
           family = beta_binomial2,  stanvars = stanvars,
           prior = c(set_prior ("normal(0,2.5)",class="b"),
                     set_prior ("normal(0,2.5)",class="Intercept")),
           chains = 4,cores=4,
           iter = 4000,control = list(adapt_delta =0.99 ,max_treedepth=15) ,
           warmup = 1000)
#prior 5
A_d2<- brm(Prey_E|vint(Prey_D)~ 1 + Treatment+(1|Block)
           ,data = A_dors2,
           family = beta_binomial2,  stanvars = stanvars,
           prior = c(set_prior ("normal(0,5)",class="b"),
                     set_prior ("normal(0,5)",class="Intercept")),
           chains = 4,cores=4,
           iter = 4000,control = list(adapt_delta =0.99 ,max_treedepth=15) ,
           warmup = 1000)
#prior 10
A_d10<- brm(Prey_E|vint(Prey_D)~ 1 + Treatment+(1|Block)
           ,data = A_dors2,
           family = beta_binomial2,  stanvars = stanvars,
           prior = c(set_prior ("normal(0,10)",class="b"),
                     set_prior ("normal(0,10)",class="Intercept")),
           chains = 4,cores=4,
           iter = 4000,control = list(adapt_delta =0.99 ,max_treedepth=15) ,
           warmup = 1000)

##posterior distribution from the priors on the count scale 
n_Ad_p <-data.frame(ABase=0,Aexp_1=0,
                     BBase=0,Bexp_1=0,
                     CBase=0,Cexp_1=0,
                    DBase=0,Dexp_1=0)

#DIFFERENT PRIOR MODEL OUTPUTS
n_Ad_p[c(1:12000),1:2]<-   mod_out(A_d1.4)
n_Ad_p[c(1:12000),3:4]<-   mod_out(A_d2)
n_Ad_p[c(1:12000),5:6]<-   mod_out(A_d10)
n_Ad_p[c(1:12000),7:8]<-   mod_out(A_d1)


write.csv(n_Ad_p,"A_dorsR.csv")
#####################Pterostichus madidus##########################
#####################Pterostichus madidus##########################
#####################Pterostichus madidus########################## 
P_madi$orle<-paste0(seq(1:nrow(P_madi)),"A")
#prior 2.5
P_md1.4<- brm(Prey_E|vint(Prey_D)~ 1 + Treatment+(1|Block)
              ,data = P_madi,
              family = beta_binomial2,stanvars=stanvars,
              prior = c(set_prior ("normal(0,2.5)",class="b"),
                        set_prior ("normal(0,2.5)",class="Intercept")),
              chains = 4,cores=4,
              iter = 4000,control = list(adapt_delta =0.99 ,max_treedepth=15) ,
              warmup = 1000)
#prior 1
P_md1<- brm(Prey_E|vint(Prey_D)~ 1 + Treatment+(1|Block)
              ,data = P_madi,
              family = beta_binomial2,stanvars=stanvars,
              prior = c(set_prior ("normal(0,1)",class="b"),
                        set_prior ("normal(0,1)",class="Intercept")),
              chains = 4,cores=4,
              iter = 4000,control = list(adapt_delta =0.99 ,max_treedepth=15) ,
              warmup = 1000)
pp_check(P_md1.4b,nsamples=100)
#prior 5
P_md2<- brm(Prey_E|vint(Prey_D)~ 1 + Treatment+(1|Block)
               ,data = P_madi,
               family = beta_binomial2,stanvars=stanvars,
               prior = c(set_prior ("normal(0,5)",class="b"),
                         set_prior ("normal(0,5)",class="Intercept")),
               chains = 4,cores=4,
               iter = 4000,control = list(adapt_delta =0.99 ,max_treedepth=15) ,
               warmup = 1000)
#prior 10
P_md10<- brm(Prey_E|vint(Prey_D)~ 1 + Treatment+(1|Block)
               ,data = P_madi,
               family = beta_binomial2,stanvars=stanvars,
               prior = c(set_prior ("normal(0,10)",class="b"),
                         set_prior ("normal(0,10)",class="Intercept")),
               chains = 4,cores=4,
               iter = 4000,control = list(adapt_delta =0.99 ,max_treedepth=15) ,
               warmup = 1000)

##posterior distribution from the priors on the count scale 
n_Pmd_p <-data.frame(ABase=0,Aexp_1=0,Aexp_5=0,
                     BBase=0,Bexp_1=0,Bexp_5=0,
                     CBase=0,Cexp_1=0,Cexp_5=0,
                     DBase=0,Dexp_1=0,Dexp_5=0)

#DIFFERENT PRIOR MODEL OUTPUTS
n_Pmd_p[c(1:12000),1:3]<-   mod_out(P_md1.4)
n_Pmd_p[c(1:12000),4:6]<-   mod_out(P_md2)
n_Pmd_p[c(1:12000),7:9]<-   mod_out(P_md10)
n_Pmd_p[c(1:12000),10:12]<- mod_out(P_md1)


write.csv(n_Pmd_p,"P_madiR.csv")
######################################################################
#####################Pterostichus melanarius##########################
###################################################################### 
P_mela$orle<-paste0(seq(1:nrow(P_mela)),"A")

#prior 1
P_me1<- brm(Prey_E|vint(Prey_D)~ 1 + Treatment+(1|Block)
              ,data = P_mela,
              family = beta_binomial2,  stanvars = stanvars,
              prior = c(set_prior ("normal(0,1)",class="b"),
                        set_prior ("normal(0,1)",class="Intercept")),
              chains = 4,cores=4,
              iter = 4000,control = list(adapt_delta =0.99 ,max_treedepth=15) ,
              warmup = 1000)
#prior 1.4
P_me1.4<- brm(Prey_E|vint(Prey_D)~ 1 + Treatment+(1|Block)
           ,data = P_mela,
           family = beta_binomial2,  stanvars = stanvars,
           prior = c(set_prior ("normal(0,2.5)",class="b"),
                     set_prior ("normal(0,2.5)",class="Intercept")),
           chains = 4,cores=4,
           iter = 4000,control = list(adapt_delta =0.99 ,max_treedepth=15) ,
           warmup = 1000)

#prior 5
P_me2<- brm(Prey_E|vint(Prey_D)~ 1 + Treatment+(1|Block)
           ,data = P_mela,
           family = beta_binomial2,  stanvars = stanvars,
           prior = c(set_prior ("normal(0,5)",class="b"),
                     set_prior ("normal(0,5)",class="Intercept")),
           chains = 4,cores=4,
           iter = 4000,control = list(adapt_delta =0.99 ,max_treedepth=15) ,
           warmup = 1000)
#prior 10
P_me10<- brm(Prey_E|vint(Prey_D)~ 1 + Treatment+(1|Block)
            ,data = P_mela,
            family = beta_binomial2,  stanvars = stanvars,
            prior = c(set_prior ("normal(0,10)",class="b"),
                      set_prior ("normal(0,10)",class="Intercept")),
            chains = 4,cores=4,
            iter = 4000,control = list(adapt_delta =0.99 ,max_treedepth=15) ,
            warmup = 1000)
 
##posterior distribution from the priors on the count scale 
n_Pme_p <-data.frame(ABase=0,Aexp_1=0,Aexp_5=0,
                     BBase=0,Bexp_1=0,Bexp_5=0,
                     CBase=0,Cexp_1=0,Cexp_5=0,
                     DBase=0,Dexp_1=0,Dexp_5=0)

#DIFFERENT PRIOR MODEL OUTPUTS
n_Pme_p[c(1:12000),1:3]<-   mod_out(P_me1.4)
n_Pme_p[c(1:12000),4:6]<-   mod_out(P_me2)
n_Pme_p[c(1:12000),7:9]<-   mod_out(P_me10)
n_Pme_p[c(1:12000),10:12]<- mod_out(P_me1)

write.csv(n_Pme_p,"P_melaR.csv")

###############################################################
#####################Harpalus rufipes##########################
############################################################### 

H_rufi$orle<-paste0(seq(1:nrow(H_rufi)),"A")

#prior 2.5
H_rufs1.4<- brm(Prey_E|vint(Prey_D)~ 1 + Treatment+(1|Block)
            ,data = H_rufi,
            family = beta_binomial2,  stanvars = stanvars,
            prior = c(set_prior ("normal(0,2.5)",class="b"),
                      set_prior ("normal(0,2.5)",class="Intercept")),
            chains = 4,cores=4,
            iter = 4000,control = list(adapt_delta =0.99 ,max_treedepth=15) ,
            warmup = 1000)
#prior 1
H_rufs1<- brm(Prey_E|vint(Prey_D)~ 1 + Treatment+(1|Block)
                ,data = H_rufi,
                family = beta_binomial2,  stanvars = stanvars,
                prior = c(set_prior ("normal(0,1)",class="b"),
                          set_prior ("normal(0,1)",class="Intercept")),
                chains = 4,cores=4,
                iter = 4000,control = list(adapt_delta =0.99 ,max_treedepth=15) ,
                warmup = 1000)
#prior 5
H_rufs2<- brm(Prey_E|vint(Prey_D)~ 1 + Treatment+(1|Block)
                ,data = H_rufi,
                family = beta_binomial2,  stanvars = stanvars,
                prior = c(set_prior ("normal(0,5)",class="b"),
                          set_prior ("normal(0,5)",class="Intercept")),
                chains = 4,cores=4,
                iter = 4000,control = list(adapt_delta =0.99 ,max_treedepth=15) ,
                warmup = 1000)
#prior 10
H_rufs10<- brm(Prey_E|vint(Prey_D)~ 1 + Treatment+(1|Block)
              ,data = H_rufi,
              family = beta_binomial2,  stanvars = stanvars,
              prior = c(set_prior ("normal(0,10)",class="b"),
                        set_prior ("normal(0,10)",class="Intercept")),
              chains = 4,cores=4,
              iter = 4000,control = list(adapt_delta =0.99 ,max_treedepth=15) ,
              warmup = 1000)

##posterior distribution from the priors on the count scale 
n_Hrf_p <-data.frame(ABase=0,Aexp_1=0,Aexp_5=0,
                     BBase=0,Bexp_1=0,Bexp_5=0,
                     CBase=0,Cexp_1=0,Cexp_5=0,
                     DBase=0,Dexp_1=0,Dexp_5=0)

#DIFFERENT PRIOR MODEL OUTPUTS
n_Hrf_p[c(1:12000),1:3]<-   mod_out(H_rufs1.4)
n_Hrf_p[c(1:12000),4:6]<-   mod_out(H_rufs2)
n_Hrf_p[c(1:12000),7:9]<-   mod_out(H_rufs10)
n_Hrf_p[c(1:12000),10:12]<- mod_out(H_rufs1)


write.csv (n_Hrf_p,"H_rufiR.csv")

#####################Coccinellidae sept##########################
#####################Coccinellidae sept##########################
#####################Coccinellidae sept########################## 
C_sept$orle<-paste0(seq(1:nrow(C_sept)),"A")

#prior 1
c_s11<- brm(Prey_E|vint(Prey_D)~ 1 + Treatment+(1|Block)
             ,data = C_sept,
             family = beta_binomial2,  stanvars = stanvars,
             prior = c(set_prior ("normal(0,1)",class="b"),
                       set_prior ("normal(0,1)",class="Intercept")),
             chains = 4,
             iter = 4000,control = list(adapt_delta =0.99 ,max_treedepth=15) ,
             warmup = 1000)
#prior 2.5
c_s1.4<- brm(Prey_E|vint(Prey_D)~ 1 + Treatment+(1|Block)
             ,data = C_sept,
             family = beta_binomial2,  stanvars = stanvars,
             prior = c(set_prior ("normal(0,2.5)",class="b"),
                       set_prior ("normal(0,2.5)",class="Intercept")),
             chains = 4,
             iter = 4000,control = list(adapt_delta =0.99 ,max_treedepth=15) ,
             warmup = 1000)
#prior 5
c_s2<- brm(Prey_E|vint(Prey_D)~ 1 + Treatment+(1|Block)
           ,data = C_sept,
           family = beta_binomial2,  stanvars = stanvars,
           prior = c(set_prior ("normal(0,5)",class="b"),
                     set_prior ("normal(0,5)",class="Intercept")),
           chains = 4,
           iter = 4000,control = list(adapt_delta =0.99 ,max_treedepth=15) ,
           warmup = 1000)
#prior 10
c_s10<- brm(Prey_E|vint(Prey_D)~ 1 + Treatment+(1|Block)
           ,data = C_sept,
           family = beta_binomial2,  stanvars = stanvars,
           prior = c(set_prior ("normal(0,10)",class="b"),
                     set_prior ("normal(0,10)",class="Intercept")),
           chains = 4,
           iter = 4000,control = list(adapt_delta =0.99 ,max_treedepth=15) ,
           warmup = 1000)

##posterior distribution from the priors on the count scale 
n_c_s_p <-data.frame(ABase=0,Aexp_1=0,Aexp_5=0,
                     BBase=0,Bexp_1=0,Bexp_5=0,
                     CBase=0,Cexp_1=0,Cexp_5=0,
                     DBase=0,Dexp_1=0,Dexp_5=0)

#DIFFERENT PRIOR MODEL OUTPUTS
n_c_s_p[c(1:12000),1:3]<-   mod_out(c_s1.4)
n_c_s_p[c(1:12000),4:6]<-   mod_out(c_s2)
n_c_s_p[c(1:12000),7:9]<-   mod_out(c_s10)
n_c_s_p[c(1:12000),10:12]<- mod_out(c_s11)


write.csv(n_c_s_p,"C_septR.csv")
############################################################
#####################Amara plebja###########################
############################################################

A_pleb$orle<-paste0(seq(1:nrow(A_pleb)),"A")
#prior 2.5
a_pl1.4<- brm(Prey_E|vint(Prey_D)~ 1 + Treatment+(1|Block)
           ,data = A_pleb,
           family = beta_binomial2,  stanvars = stanvars,
           prior = c(set_prior ("normal(0,2.5)",class="b"),
                     set_prior ("normal(0,2.5)",class="Intercept")),
           chains = 4,cores=4,
           iter = 4000,control = list(adapt_delta =0.99 ,max_treedepth=15) ,
           warmup = 1000)
#prior 1
a_pl1<- brm(Prey_E|vint(Prey_D)~ 1 + Treatment+(1|Block)
             ,data = A_pleb,
             family = beta_binomial2, stanvars=stanvars,
             prior = c(set_prior ("normal(0,1)",class="b"),
                       set_prior ("normal(0,1)",class="Intercept")),
             chains = 4,cores=4,
             iter = 4000,control = list(adapt_delta =0.99 ,max_treedepth=15) ,
             warmup = 1000)
pp_check(a_pl1,nsamples=100)
#prior 5
a_pl2<- brm(Prey_E|vint(Prey_D)~ 1 + Treatment+(1|Block)
              ,data = A_pleb,
              family = beta_binomial2,  stanvars = stanvars,
              prior = c(set_prior ("normal(0,5)",class="b"),
                        set_prior ("normal(0,5)",class="Intercept")),
              chains = 4,cores=4,
              iter = 4000,control = list(adapt_delta =0.99 ,max_treedepth=15) ,
              warmup = 1000)
#prior 10
a_pl10<- brm(Prey_E|vint(Prey_D)~ 1 + Treatment+(1|Block)
              ,data = A_pleb,
              family = beta_binomial2,  stanvars = stanvars,
              prior = c(set_prior ("normal(0,10)",class="b"),
                        set_prior ("normal(0,10)",class="Intercept")),
              chains = 4,cores=4,
              iter = 4000,control = list(adapt_delta =0.99 ,max_treedepth=15) ,
              warmup = 1000)

##posterior distribution from the priors on the count scale 
n_apl_p <-data.frame(ABase=0,Aexp_1=0,Aexp_5=0,
                     BBase=0,Bexp_1=0,Bexp_5=0,
                     CBase=0,Cexp_1=0,Cexp_5=0,
                     DBase=0,Dexp_1=0,Dexp_5=0)

#DIFFERENT PRIOR MODEL OUTPUTS
n_apl_p[c(1:12000),1:3]<-   mod_out(a_pl1.4)
n_apl_p[c(1:12000),4:6]<-   mod_out(a_pl2)
n_apl_p[c(1:12000),7:9]<-   mod_out(a_pl10)
n_apl_p[c(1:12000),10:12]<- mod_out(a_pl1)

write.csv(n_apl_p,"A_plebR.csv")

#species predation estimates#################################################################################
Ap<-read.csv("A_plebR.csv",header=T)

Ad<-read.csv("A_dorsR.csv",header=T)

Bb<- read.csv("B_bullR.csv",header=T)

Hr<- read.csv ("H_rufiR.csv",header=T)

Ha<-read.csv("H_affR.csv",header=T)

Pc<-read.csv("P_cuprR.csv",header=T)

Pm<-read.csv("P_melaR.csv",header=T)

Pd<-read.csv("P_madiR.csv",header=T)

Nb<-read.csv("N_brevR.csv",header=T)

Cs<-read.csv("C_septR.csv",header=T)

Hx <- read.csv("H_axyR.csv",header=T) 

Pg<- read.csv("P_cognR.csv",header=T)

###SETUP OUTPUT DATAFRAME FOR FIGURE 1
Pred2.5<-data.frame(Species=0,Exposure=0,Mean_difference=0,Diff_lb=0,Diff_ub=0,Control=0,C_LB=0,U_LB=0)

#THIS FUNCTION EXTRACTS THE ESTIMATE BASED ON  mean = 0 and s.d. = 2.5 PRIOR
#x is a dataframe where Abase refers to control, exp1 = resistance and exp5 = recovery
#s refers to the species name

spec_est <-function (x,s){
  
  #posterior distributions are based on probability of consumptions and needs to be converted back to log odds
  #to retrieve model estimates
  # control
  XC2.5 <- logit(x$ABase/20)
  # resistance
  X12.5 <- logit(x$Aexp_1/20) - logit(x$ABase/20)
  # recovery
  X52.5 <- logit(x$Aexp_5/20) - logit(x$ABase/20)
  
  #this then recovers the mean and 95% credible intervals and puts them into a dataframe
  OUT<-data.frame(Species=0,Exposure=0,Mean_difference=0,Diff_lb=0,Diff_ub=0,Control=0,C_LB=0,U_LB=0)
  #species name
  OUT[1:2,1]<-s
  #treatment
  OUT[1,2]<-"Resistance"
  OUT[2,2]<-"Recovery"
  #resistance estimate
  OUT[1,3:5]<- round( c(mean(X12.5),posterior_interval(as.matrix(X12.5),prob=0.95)),2)
  #recovery estimate
  OUT[2,3:5]<- round( c(mean(X52.5),posterior_interval(as.matrix(X52.5),prob=0.95)),2)
  #control mean estimate populated for each row
  OUT[1,6:8]<- round( c(mean(XC2.5),posterior_interval(as.matrix(XC2.5),prob=0.95)),2)
  OUT[2,6:8]<- round( c(mean(XC2.5),posterior_interval(as.matrix(XC2.5),prob=0.95)),2)
  return(OUT)}
###################################################################################################
#POPULATE OUTPUT DATAFRAME WITH SPECIES ESTIMATES
Pred2.5[c(1:2),] <- spec_est(Ha,"Harpalus affinis")

Pred2.5[c(3:4),] <- spec_est(Hx,"Harmonia axyridis")

Pred2.5[c(5:6),] <- spec_est(Nb,"Nebria brevicollis")

#manually calculate badister bullatus as it does not have a recovery estimate
BbC2.5 <- logit(Bb$ABase/20)
Bb12.5 <- logit(Bb$Aexp_1/20) - logit(Bb$ABase/20)
Bb52.5 <- logit(Bb$Aexp_5/20) - logit(Bb$ABase/20)

Pred2.5[7,1]<-"Badister bullatus"
Pred2.5[7,2]<-"Resistance"
Pred2.5[7,3:5]<- round( c(mean(Bb12.5),posterior_interval(as.matrix(Bb12.5),prob=0.95)),2)
Pred2.5[7,6:8]<- round( c(mean(BbC2.5),posterior_interval(as.matrix(BbC2.5),prob=0.95)),2)
#
Pred2.5[c(8:9),] <- spec_est(Pg,"Philonthus cognatus")

Pred2.5[c(10:11),] <- spec_est(Pc,"Poecilus cupreus")

#manually calculate anchomenus dorsalis as it does not have a recovery estimate
AdC2.5 <- logit(Ad$ABase/20)
Ad12.5 <- logit(Ad$Aexp_1/20) - logit(Ad$ABase/20)
Pred2.5[12,1]<-"Anchomenus dorsalis"
Pred2.5[12,2]<-"Resistance"
Pred2.5[12,3:5]<- round( c(mean(Ad12.5),posterior_interval(as.matrix(Ad12.5),prob=0.95)),2)
Pred2.5[12,6:8]<- round( c(mean(AdC2.5),posterior_interval(as.matrix(AdC2.5),prob=0.95)),2)  

Pred2.5[c(13:14),] <- spec_est(Pd,"Pterostichus madidus")

Pred2.5[c(15:16),] <- spec_est(Pm,"Pterostichus melanarius")

Pred2.5[c(17:18),] <- spec_est(Hr,"Harpalus rufipes")

Pred2.5[c(19:20),] <- spec_est(Cs,"Coccinella septempunctata")

Pred2.5[c(21:22),] <- spec_est(Ap,"Amara plebja")

#for figures
dodge<-position_dodge(width=0.5)
#REMOVE UNIQUE ESTIMATES FOR EACH SPECIES CONTROL ESTIMATE
con1<-distinct(data.frame(Pred2.5$Species, Pred2.5$Control,Pred2.5$C_LB,Pred2.5$U_LB))
con1<- con1[order(con1[2]),]
con1


#FIGURE FOR DIFFERENCES FROM THE CONTROL
all1<-ggplot (Pred2.5,aes(x=reorder(Species,Control),y=Mean_difference,color=Exposure))+ 
  geom_point(aes(shape=Exposure),position= dodge)+
  geom_errorbar(aes(ymin=Diff_lb,ymax=Diff_ub),position=dodge)+ 
  geom_hline(yintercept = 0, linetype="dashed", 
             color = "black", size=0.5)+coord_flip()

all1


###plots of intercept for different prior distributions
p1<-ggplot (con1,aes(x=reorder(Pred2.5.Species,Pred2.5.Control),y=Pred2.5.Control))+ 
  geom_bar(stat="identity")+
  geom_errorbar(aes(ymin=Pred2.5.C_LB,ymax=Pred2.5.U_LB),position=dodge)+coord_flip()
p1

#graphs for paper main graph
p1<-p1+labs(y="Control log odds of predation (?95 CI)",x="")+ 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(),text=element_text(size=12) ,axis.line = element_line(colour = "black"))
all1<-all1+labs(y="Difference in log odds of predation from the control (?95 CI)",x="")+
  theme(panel.grid.major = element_blank(),text=element_text(size=12), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
library(ggpubr)
p1
all1
ggsave("Pred_plotR.png",plot=plots1, dpi=400, dev='png', height=6, width=8, units="in") 
####Alive and dead data
####Alive and dead data
ALL <-read.csv("ALL_FEEDING.csv",header=T)

sp <-table(ALL$species_name,ALL$Treatment, ALL$ID_2)

sp <-as.data.frame(sp)
sp0<-subset(sp,sp$Var3=="DEAD")
sp1<-subset(sp,sp$Var3=="ALIVE")
sp1$Freq0<-sp0$Freq
sp1$Tot<- sp1$Freq+sp1$Freq0

sp1$Freq0<-sp0$Freq

sp1$p<-round((sp1$Freq/sp1$Tot)*100,2)
D1<-distinct(Pred2.5[c(1,6,7,8)])
colnames(sp1)[1]<-"Species"
levels(sp1$Species)[1:12]<-c("Amara plebja","Anchomenus dorsalis","Badister bullatus",
                          "Coccinella septempunctata","Harmonia axyridis","Harpalus affinis",
                          "Harpalus rufipes","Nebria brevicollis","Philonthus cognatus",
                          "Poecilus cupreus", "Pterostichus madidus", "Pterostichus melanarius")
sp1 <-merge(sp1,D1,by="Species")




sp1$Var2 = factor(sp1$Var2,levels(sp1$Var2)[c(2,3,1)])
sp1<-sp1[-c(5,8),]

###plots of intercept for different prior distributions
p3<-ggplot (sp1,aes(x=reorder(sp1$Species,sp1$Control),y=sp1$p),group=sp1$Var2)+ 
  geom_bar(stat="identity",position=position_dodge(),
           aes(x=reorder(sp1$Species,sp1$Control),fill=sp1$Var2))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(),text=element_text(size=10) ,
        axis.text.x = element_text(angle = 90),axis.line = element_line(colour = "black"))+ 
  guides(fill = guide_legend(reverse = TRUE))+ylab("Survival (%)")+xlab("Species")
p3$labels$fill<-"Treatment"
p3


plots1<- ggarrange(p1,all1,p3,ncol=1,labels=c("2a)","2b)","2c"),font.label=list (face="plain")) 
plots1
