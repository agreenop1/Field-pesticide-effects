setwd("C:/Users/arrgre/iCloudDrive/PHD WORK/Field pest control/Paper_1/Analysis")
rm(list = ls(all.names = TRUE))

#packages
library(boot)
library(brms)
library(data.table)
library(ggplot2)
library(bayesplot)
library(ape)
library(broom) 
library(projpred)
library(vegan)
library(FD)
library(picante)
library(adephylo)

# ABUNDANCE OF ALL THE PREDATORS
red_pred<-read.csv("Abun1.csv",header=T)

# just abundances in a dataframe
red_pred1<-read.csv("Abun1.csv",header=T)

#DIVERSITY#####################################################################
J<-diversity(red_pred1)/log (specnumber(red_pred1))#this is community evenness
#abundance
red_pred$s_ab<- apply (red_pred1,1,sum)
#proportion
red_pred$prop <- round((red_pred$s_ab/red_pred$abun1)*100,2)
#species richness
red_pred$sr<-specnumber(red_pred1)

#PHYLOGENETIC DIVERSITY - the mean pairwise standardised PD based on sum branch length
#phylo file
P1<- read.csv("Phylo1.csv",header=T)
attach(P1)
phy2 <- as.phylo (~Order/Family/sub_family/tribe/Genus/Species, P1 ,collapse=F)#order for the phylo tree
plot(phy2)
phy3 <- as.matrix(  distTips (phy2,method = c("patristic"),tips="all"))

#phlo community standardused
sesmpd<-ses.mpd(as.matrix(red_pred1),phy3,null.model = c("taxa.labels"),
                runs = 999, iterations = 1000, abundance.weighted = TRUE) 

#FUNCTIONAL DIVERSITY
#traits
fds2  <- read.csv("TRAITS.csv",header=T,row.names=1)
#calculate functional diversity
allfd <- fdisp(gowdis(fds2),as.matrix(red_pred1))
#calculate the CWM
mass_g<-functcomp(fds2[1],as.matrix(red_pred1),CWM.type = "all")
wings_B<-functcomp(fds2[2],as.matrix(red_pred1),CWM.type = "all")
wings_M<-functcomp(fds2[3],as.matrix(red_pred1),CWM.type = "all")
wings_D<-functcomp(fds2[4],as.matrix(red_pred1),CWM.type = "all")

#############################################################################################
#FEEDING CAPACITY DATA ( at least two lists per species)
# amara
Ap<-read.csv("A_pleb.csv",header=T)[1:4]
Ap$ID <-"A_pleb"
A_pleb<- list ( c(Ap$ABase),c(Ap$Aexp_1))

# anchomenus
Ad<-read.csv("A_dors.csv",header=T)[1:3]
Ad$ID <-"A_dors"
A_dors<- list ( c(Ad$ABase),c(Ad$Aexp_1))

#badister
Bb<- read.csv("B_bull.csv",header=T)[1:3]
Bb$ID<-"B_bull"
B_bipu<- list ( c(Bb$ABase),c(Bb$Aexp_1))

#harpalus r
Hr<- read.csv ("H_rufi.csv",header=T)[1:4]
Hr$ID<-"H_rufi"
H_rufi<- list ( c(Hr$ABase),c(Hr$Aexp_1))

#harpalus a
Ha<-read.csv("H_aff.csv",header=T)[1:4]
Ha$ID<-"H_affi"
H_affi<- list ( c(Ha$ABase),c(Ha$Aexp_1))

#poecilius c
Pc<-read.csv("P_cupr.csv",header=T)[1:4]
Pc$ID<- "P_cupr"
P_cupr<- list ( c(Pc$ABase),c(Pc$Aexp_1))

#pterostichus m
Pm<-read.csv("P_mela.csv",header=T)[1:4]
Pm$ID<- "P_mela"
P_mela<- list ( c(Pm$ABase),c(Pm$Aexp_1))

#pterostichus d
Pd<-read.csv("P_madi.csv",header=T)[1:4]
Pd$ID<-"P_madi"
P_madi<- list ( c(Pd$ABase),c(Pd$Aexp_1))

#nebria brevicollis
Nb<-read.csv("N_brev.csv",header=T)[1:4]
Nb$ID<- "N_brev"
N_brev<- list ( c(Nb$ABase),c(Nb$Aexp_1))

#cocc sept
Cs<-read.csv("C_sept.csv",header=T)[1:4]
Cs$ID<-"C_sept"
C_sept<- list ( c(Cs$ABase),c(Cs$Aexp_1))

#########################################################################################################################################################
######################################FUNCTION TO CALCULATE COMMUNTIY DIFFERENCES
# for a selected predator this samples the functional capacity of predator by 
# the abundances across all sites 
# y = abundances of individual  | s = predation sampling distribution
pest_eff <-  function (y,s) {
  #if abundance is zero no sampling occurs and onto the next site
  if  ((y)>0){
      #sampe predation by abundance
    a0 <-  sample(s[[1]],(y),replace = T)
    a1 <-  sample(s[[2]],(y),replace= T)
  } else  {
      #absent so does not impact abundance
    a0<- 0
    a1<- 0
  }
  fin <- list (a0,a1)  
  return(fin)
}

#nested function one to amalgamate all pesticide responses
fin_counts <- function (x) {
  ##########################
  #all abundances for a named predator (calls abundance column based on name)
  op <- all_pred_d[x] 
  #get the functional capacity of predator from global environment
  samps <- get (noquote(x))
  #for this species calculate the number of aphids consumed across all sites
  all_abs <- apply (op,1,pest_eff,s=samps)
  #sum to total predation done by tht species at a0 and a1
  fin_pred0<-data.frame()
  fin_pred1<-data.frame()
  for (i in 1:nrow(all_pred_d)) {
    #this puts all the predator responses based on name into a single list - condition1
    fin_pred0[i,paste(noquote(x),0)]  <- as.numeric (sum (all_abs[[i]][[1]])) 
    #this puts all the predator responses based on name into a single list - condition2 
    fin_pred1[i,paste(noquote(x),1)]  <- as.numeric(sum (all_abs[[i]][[2]])) 
  } 
  fin_pred <-list (fin_pred0, fin_pred1)
  return(fin_pred)}
##########function for community differences defined below######## 

##function to calculate differences 
comm_diffs <- function (x){
  #returns list of all predation by all predators
  comps<- lapply (x,fin_counts)
  h0<-data.frame(matrix(nrow=nrow(all_pred_d)))
  h1<-data.frame(matrix(nrow=nrow(all_pred_d)))
  #extract total predation of all species for H0 and H1 !Species number needs to be changed based on subset!
  #this gives all responses at condition 1 and 2 in a single dataframe
  for (i in 1:length(x)) { 
    asd1<- as.data.frame( comps[[i]][1])
    h0[colnames(asd1)[1]] <- asd1[1]
    asd<- as.data.frame( comps[[i]][2]) 
    h1[colnames(asd)[1]] <- asd[1]
  } 
  #
  #calculate difference in predation between condition 0 and 1
  h0s <- data.frame(rowSums(h0[2:ncol(h0)]),rowSums(h1[2:ncol(h1)] ) )
  h0s$diff <- (h0s$rowSums.h1.2.ncol.h1.../h0s$rowSums.h0.2.ncol.h0...)
  
  return (list(h0s,h0,h1))}
############function end and community calculation

#abs is used to describe abundance of all species and must be included in 
# also in all_pred_d which is used within the comms diff function
abs<-round (red_pred1)
all_pred_d<-abs
all_names<-names(abs)
all_names
#data frames for the function output
comm_eff<-list()
#all_pred_d needs to be a dataframe with the abundances of all predators in 
#a list of the names of the predators then needs to be passed to the comm diff function
# each element of this list must have an associated functional predator variation
# data frame
for (j in 1:100){
  b<-comm_diffs(all_names)[[1]]
  b$ID<-paste0(j,"I")
  comm_eff[[j]] <-b
  
}

#all the raw values for each community (rep 100 times for every dataset)
comm_eff<-rbindlist(comm_eff)
#eve
comm_eff$eve<-rep(J,100)
#abs
comm_eff$abs<-rep(log(red_pred$s_ab),100) 
#species rich
comm_eff$sr <-rep(red_pred$sr,100)
#phy
comm_eff$PHY<-rep(sesmpd$mpd.obs.z,100)
#wing types
comm_eff$wings_M <- rep(wings_M$Wings_M_1,100)
comm_eff$wings_D <- rep(wings_D$WINGS_D_1,100)
comm_eff$wings_B <- rep(wings_B$Wings_B_1,100)
#mass
comm_eff$mas<-rep ((mass_g$Mass_g),100)
#FD
comm_eff$Fdisp<-rep(allfd$FDis,100)

##option below to use z scores of variables 
#!!!CHECK LOG ABUNDANCE OF PESTICIDE EXPOSURE!!!
##option below to use z scores of variables
#eve
comm_eff$eve_z<-   rep(as.numeric(scale(J,center = T,scale=T)),100)
#abs
comm_eff$abs_z<-   rep(as.numeric(scale(log(red_pred$s_ab),center = T,scale=T)),100)
#sr
comm_eff$sr_z<-    rep(as.numeric(scale(red_pred$sr,center = T,scale=T)),100)
#phy
comm_eff$PHY_Z<-   rep(as.numeric(scale(sesmpd$mpd.obs.z,center = T,scale=T)),100)
#mass
comm_eff$Mass_z<-  rep(as.numeric(scale(mass_g$Mass_g,scale=T,center=T)),100)
#FD
comm_eff$Fdisp_z<- rep(as.numeric(scale(allfd$FDis,scale=T,center=T)),100)
#site
comm_eff$site<- rep(paste0 (1:256,"A"),100)
comm_eff1 <-comm_eff



##REFERENCE MODEL###################################################
library(brms)
library(rstan)
library(rstanarm)
library(projpred)
library(doParallel)
library(parallel)
library(foreach)

######set up correct dataframes for output of comm models and species models
#remove sites with NA
comps1 <- droplevels(subset(comm_eff1,subset=comm_eff1$site!="250A") )
#set up correct row name
rownames(comps1)<- seq(length=nrow(comps1))
#split data into 100 even datasets
samps<- seq(1,25500,255)
comps1$run<- paste0(rep (1:100,each=255),"a")
all_dat <-split(comps1,comps1$run)

#set up output for both models
res_mod_NEW<-vector("list",100)
#set up list for R2 values
R2_Fd<-list()

##set up cores##### 
numCores <- 10
cl <- makeCluster(numCores)
registerDoParallel(cl)

#function to run through all data sets using the same basic models for community and species differences
#runs in parallel on computer cores !!can take a long time and use a lot of memory depending on specs!!

foreach (i = 1:100,.packages = c("rstan","rstanarm","projpred"),
         .combine = rbind) %dopar%{
           
res_mod_NEW <- stan_glm( log(diff )~ 1 +sr_z+abs_z+PHY_Z+Mass_z+Fdisp_z
                         ,data = as.data.frame(all_dat[[i]]),
                         family = gaussian(), prior=hs(global_scale = 1/(5-1) * 1/sqrt(255) , slab_scale = 1) 
                         ,prior_intercept = normal(0,10),
                         chains = 4,warmup = 1000,cores=1,
                         iter = 3000,control = list(adapt_delta =0.99 ,max_treedepth=13))

#this make sures the model selection (in this case k-fold) only runs on the current core 
#it may be faster or slower to parallelise this based on the computer or cluster used

options(mc.cores =1)
l<-cv_varsel(res_mod_NEW, cv_method='kfold',K=10, selection="forward")
b<- list(res_mod_NEW ,R2_Fd,l) # model object currently stored which uses a lot of space

#output for models 
nm<-paste0("OUT/","ITER",i,".rds")
saveRDS(b,nm)

}
stopCluster(cl)



#################################################################################################################################
###############################################MODEL COEFFICIENT OUTPUTS FROM 100 ITERATIONS#####################################
#################################################################################################################################
setwd("C:/Users/arrgre/iCloudDrive/PHD WORK/Field pest control/Paper_1/Analysis")
#####run funciton for base models across all simulated pest control values
library(data.table)
library(projpred)
library(brms)
library(dplyr)

#output files
lps<-list()
mods<-list()
dt<-list()

#cycle through each interation and pull out variable selection
for (i in 1:100){
  
#iteration !CHECK FILE!
nm<-paste0("OUT/","ITER",i,".rds")
fin<- readRDS(nm)

#suggest how many variables should be included in the subset
nv<-suggest_size(fin[[3]],alpha=0.1)

#project the posterior distribution based on that subset
lps[[i]] <-project(fin[[3]],nv=nv,ns=2000)

#iteration number
d<-i

#output for the model estimates for each variable subset for each iteration
mod <- data.frame (matrix (nrow=10,ncol=4))

  #complete models into indivudal data frame
  #length
  for (y in 1:length ( round(colMeans(as.matrix(lps[[d]]),2)))){
    
    #name
    mod[y,1]<-names(( round(colMeans(as.matrix(lps[[d]])),2)))[y]
    
    #mean estimate
    mod[y,2]<-round(colMeans(as.matrix(lps[[d]])),2)[y]
    
    #credible intervals
    mod[y,3:4]<- round(posterior_interval(as.matrix(lps[[d]])),2)[y,1:2]
  }
  mods[[i]]<-mod
  dt[[i]] <-fin[[1]]$data
}

setwd("C:/Users/arrgre/iCloudDrive/PHD WORK/Field pest control/Paper_1/Analysis/1exp_res")

##save output########
savo<-list (mods,lps,dt)
saveRDS(savo,"RESIL.rds")


###########################GET COEFFICENTS FOR MAIN TABLE IN THE PAPER###################
#these functions may only work on the output produced by the analysis in the paper
#check read in file!
dats <-readRDS("RESIL.rds")

#coef range function (actually used in the below function)
coef_ran<- function (v,x){
  #subset dataframe based on variable
  x1 <-subset(x,x$X1==v)
  
  #order by the coefficient
  x1<- x1[order(x1$X2,x1$X3),]
  
  #get min
  m_min<-x1[1,]
  
  #order by the coefficient
  x1<- x1[order(x1$X2,x1$X4),]
  
  #get max
  m_max<-x1[nrow(x1),]
  m_min$Ran<-"Min"
  m_max$Ran<-"Max"
  
  #both files put in rbind
  out<- as.data.frame(rbind(m_min,m_max))
  #
  out<-out %>%
    select_if(~ !any(is.na(.)))
  out$Inc<-nrow(x1)
  
  #reorder dataframe
  out<-out[,c("X1","Ran","X2","X3","X4","Inc")]
   return(out)
}
#################################################
#table output function for the diversity measure estimates 

proj_fun_t<- function(x){
  #all model coefficients
  mod_coef  <- rbindlist(x)
  
  #variables in subset
  coef_var<-levels(as.factor(mod_coef$X1))
  coef_var <-coef_var[coef_var!="sigma"]#remove sigma as it is not presented
  out  <-rbindlist (lapply(coef_var,coef_ran,x=mod_coef))
  
  #give proper names to the variables
  out$X1<-as.factor(out$X1)
  levels(out$X1)[1:6]<-c("Intercept","Log abundance","Functional diversity",
                         "Community weighted mean body mass","Phylogenetic diversity","Species richness")
  #order to match paper
  out$od<-c(1,1,6,6,4,4,3,3,2,2,5,5)
  out<-out[order(out$od),]
  out$od<-NULL
  #give column names
  colnames(out)[1:6]<-c("Parameter","Range","Mean","Lower 95% CI","Upper 95% CI","Inclusion")
  return(out)}


#output dataframe for all the coefficients###############################################################################
out <-proj_fun_t(dats[[1]])
out
write.csv(out,"RESIL_MORT_COEFS.csv")

##########################################################################################################################
###############################################FIGURES####################################################################
##########################################################################################################################
setwd("C:/Users/arrgre/iCloudDrive/PHD WORK/Field pest control/Paper_1/Analysis")
#MAIN FIGURE FUNCTION
# v = variable , x = data frame, s -> 1 - sublethal, 2 = mortality
coef_ran<- function (v,x,s){
  library(ggplot2)
  #subset dataframe based on variable
  x1 <-subset(x,x$X1==v)
  vn <- nrow(x1)
  #order by the coefficient and get model with minimum estimate
  x1<- x1[order(x1$X2,x1$X3),]
  mi<-as.numeric(x1[1,5])
  
  #order by the coefficient max
  x1<- x1[order(x1$X2,x1$X4),]
  mx<-as.numeric(x1[nrow(x1),5][1])
  
  #variable numbers
  #subset dataframe for correct models without NA
  n1<-as.data.frame(subset(x,x$ID==mi))
  n2<-as.data.frame(subset(x,x$ID==mx))
  #get variable number subtracting the intercept and sigma
  nv1<-nrow(subset(n1,!is.na(n1$X1)))-2
  nv2<-nrow(subset(n2,!is.na(n2$X2)))-2
  
  #get the correct data from the function
  if(s==1){
    #get original data set to hold other variables at their mean
    #min estimate
    s_m_min<-readRDS(paste0("OUT/","ITER",mi,".rds"))
    #max estimate
    s_m_max<-readRDS(paste0("OUT/","ITER",mx,".rds"))} else {
      
      #min estimate
      s_m_min<-readRDS(paste0("OUT/","ITER",mi,"_mort.rds"))
      #max estimate
      s_m_max<-readRDS(paste0("OUT/","ITER",mx,"_mort.rds")) 
      print("yes")  
      
    }
  
  #dataframe to do predictions on projections
  samps1<-as.data.frame(s_m_min[[1]]$data)
  samps2<-as.data.frame(s_m_max[[1]]$data)
  
  #PUT ALL OTHER FACTORS AT ZERO
  if (v=="PHY_Z"){
    samps1[c("abs_z","sr_z","Fdisp_z","Mass_z")] <-0
    samps2[c("abs_z","sr_z","Fdisp_z","Mass_z")] <-0} else if (v=="Mass_z"){
      
      samps1[c("abs_z","sr_z","Fdisp_z","PHY_Z")]<-0
      samps2[c("abs_z","sr_z","Fdisp_z","PHY_Z")]<-0  } else if (v=="Fdisp_z"){
        
        samps1[c("abs_z","sr_z","Mass_z","PHY_Z")]<-0
        samps2[c("abs_z","sr_z","Mass_z","PHY_Z")]<-0   } else if (v=="sr_z"){
          
          samps1[c("abs_z","Fdisp_z","Mass_z","PHY_Z")]<-0
          samps2[c("abs_z","Fdisp_z","Mass_z","PHY_Z")]<-0 }else if (v=="abs_z"){
            
            samps1[c("sr_z","Fdisp_z","Mass_z","PHY_Z")]<-0
            samps2[c("sr_z","Fdisp_z","Mass_z","PHY_Z")]<-0     
          }
  
  #form predictions MIN
  #get projected linear prediction
  M_min1 <- proj_linpred(s_m_min[[3]],xnew=samps1, integrated = F,nv=nv1)
  #get posterior intervals
  M_min_pr<-round(posterior_interval(as.matrix(M_min1)),15)
  #create dataframe
  M_min_pr<-data.frame(M_min_pr)
  #put mean effect in
  M_min_pr$mean <- colMeans( proj_linpred(s_m_min[[3]],xnew=samps1, integrated = F,nv=nv1)) 
  
  #form predictions MAX
  #get projected linear prediction
  M_max1 <- proj_linpred(s_m_max[[3]],xnew=samps2, integrated = F,nv=nv2)
  #get posterior intervals
  M_max_pr<-round(posterior_interval(as.matrix(M_max1)),15)
  #create dataframe
  M_max_pr<-data.frame(M_max_pr)
  #put mean effect in
  M_max_pr$mean <- colMeans( proj_linpred(s_m_max[[3]],xnew=samps2, integrated = F,nv=nv2))
  
  #Correct x variable - non standardised
  if (v=="PHY_Z"){xs="PHY" } else if (v=="Mass_z"){
    xs="mas"} else if (v=="Fdisp_z"){
      xs="Fdisp"  } else if (v=="sr_z"){
        xs="sr" }else if (v=="abs_z"){
          xs="abs"}
  
  
  #create a joint dataframe
  M_min_pr$Estimate<-"Min"
  M_min_pr$x<-samps1[,xs]
  
  M_max_pr$Estimate<-"Max"
  M_max_pr$x<-samps2[,xs]
  
  #ALL DATA USED FOR THE PLOT
  FDAT<-rbind.data.frame(M_min_pr,M_max_pr)
  
  #PLOT
  
  #Correct x variable 
  if (v=="PHY_Z"){xn="Phylogenetic diversity" } else if (v=="Mass_z"){
    xn="Community weighted mean body mass"} else if (v=="Fdisp_z"){
      xn="Functional diversity"  } else if (v=="sr_z"){
        xn="Species richness" }else if (v=="abs_z"){
          xn="Log abundance"}
  if (vn>1){
  #FINAL PLOT
  p1 <-ggplot(data=FDAT) +  geom_ribbon(aes(ymin=as.numeric(FDAT$X5.),
                                            ymax=as.numeric(FDAT$X95.),
                                            x=as.numeric(FDAT$x),fill=FDAT$Estimate),alpha=0.4)+
    geom_line(aes(x=as.numeric(FDAT$x),y=as.numeric(FDAT$mean),group=as.factor(FDAT$Estimate),color=FDAT$Estimate)) +
    labs(y = 'Resistance log response ratio', x = xn )+
    theme (axis.text.y = element_text(size=8),
           axis.text.x = element_text(size=8,angle=0),
           axis.title.x = element_text(size=9), 
           axis.title.y = element_text(size=9),
           panel.grid.major = element_blank(), 
           panel.grid.minor = element_blank(),
           panel.background = element_blank(),
           legend.title=element_text(size=9), 
           legend.text=element_text(size=8),
           plot.margin=grid::unit(c(3,0.5,0.5,3.0), "mm"),
           axis.line = element_line(colour = "black"),
           plot.title = element_text(face="bold",size=9),
           plot.subtitle=element_text(size=9))+ geom_hline(yintercept = 0,linetype="dashed")
  p1$labels$fill<-"Estimate"
  p1$labels$colour<-"Estimate"}else{
    p1 <-ggplot(data=FDAT) +  geom_ribbon(aes(ymin=as.numeric(FDAT$X5.),
                                              ymax=as.numeric(FDAT$X95.),
                                              x=as.numeric(FDAT$x)),alpha=0.4)+
      geom_line(aes(x=as.numeric(FDAT$x),y=as.numeric(FDAT$mean))) +
      labs(y = 'Resistance log response ratio', x = xn )+
      theme (axis.text.y = element_text(size=8),
             axis.text.x = element_text(size=8,angle=0),
             axis.title.x = element_text(size=9), 
             axis.title.y = element_text(size=9),
             panel.grid.major = element_blank(), 
             panel.grid.minor = element_blank(),
             panel.background = element_blank(),
             legend.title=element_text(size=9), 
             legend.text=element_text(size=8),
             plot.margin=grid::unit(c(3,0.5,0.5,3.0), "mm"),
             axis.line = element_line(colour = "black"),
             plot.title = element_text(face="bold",size=9),
             plot.subtitle=element_text(size=9))+ geom_hline(yintercept = 0,linetype="dashed")
  }
  
  return(p1)
}

fig_fun<- function(x,s){
  #all model coefficients
  mod_coef  <- rbindlist(x)
  mod_coef$ID <- rep (1:100,each=10)
  
  #variables in subset
  coef_var<-levels(as.factor(mod_coef$X1))
  coef_var <-coef_var[coef_var!="sigma"&coef_var!="(Intercept)"]#remove sigma as it is not presented + intercept

  #s = 1 sublethal , s = 2 for mortality | RETURNS ALL PLOTS FOR COEFFICIENTS
  all_p<-lapply(coef_var,coef_ran,x=mod_coef,s=s)
 
 return(all_p)
}

################ALL PLOTS#################
#fig fun 1 = sub lethal and 2 = mortalities
allp <- fig_fun(dats[[1]],s=1)
allp[[2]]
#create plot of responses for phylogenetic diversity
library(ggplot2)
library(ggpubr)

plots1<- ggarrange(allp[[4]],allp[[3]],
                   allp[[2]],allp[[5]],
                   allp[[1]],nrow =3,ncol=2,labels=c("3a)","3b)",
                                              "3c)","3d)",
                                              "3e)"),font.label=list (face="plain",size=10))
plots1
setwd("C:/Users/arrgre/iCloudDrive/PHD WORK/Field pest control/Paper_1/Analysis/1exp_res")
ggsave("RESIL_F.png",plot=plots1, dpi=400, dev='png', height=7, width=6.5, units="in")

