############################################################################### 
######  TUTORIAL TO REPRODUCE THE SIMULATIONS PUBLISHED IN MAZEL et al. #######
############################################################################### 

# Load neccessary functions and packages

source("PhyloSymbiosis_Functions.R") #functions developed 
library(parallel)
library(geiger)
library(abind)
library(vegan)
library(phytools)
library(GUniFrac)

#-------------------#
#-------------------#
#  Simple Example   #
#-------------------#
#-------------------#

#---------------------------#
#   Define parameter space  #
#---------------------------#

N_hosts=25
N_symbionts=150
V_delta_host=1
V_delta_symbiont=1
Ntraits=1
Breath=0.5
beta.comp = 0
years=20
K=150
beta.envT=.5
beta.abunT=.5
nrep=2

params=expand.grid(Simul="Example",rep=1:nrep,N_hosts=N_hosts,N_symbionts=N_symbionts,Ntraits=Ntraits,Breath=Breath,beta.comp=beta.comp,beta.abun=beta.abunT,beta.env=beta.envT,years=years,K=K,V_delta_host=V_delta_host,V_delta_symbiont=V_delta_symbiont)
Nparams=dim(params)[1]
ListParams=1:Nparams

#---------------------------#
#    Run the simulations    #
#---------------------------#

NoCores=1 #depends how many cores you can use
Res=mclapply(ListParams,Sim,params=params,mc.cores=NoCores)
Res

#--------------------------#
#--------------------------#
#  Results from the paper  #
#--------------------------#
#--------------------------#

#---------------------------#
#   Define parameter space  #
#---------------------------#

params=expand.grid(Simul=NA,rep=NA,N_hosts=NA,N_symbionts=NA,Ntraits=NA,Breath=NA,beta.comp=NA,beta.abun=NA,beta.env=NA,years=NA,K=NA,V_delta_host=NA,V_delta_symbiont=NA)

# default
N_hosts=25
N_symbionts=150
V_delta_host=c(.01,.1,1,10,1000)
V_delta_symbiont=1
Ntraits=1
Breath=0.5
beta.comp = 0
years=20
K=150
nrep=250

#S1 Varying filter Vs abundance effects
beta.envT=seq(from=0,to=1,by=.1)
beta.abunT=rev(seq(from=0,to=1,by=.1))
beta.envT=seq(from=0,to=1,by=.5)
beta.abunT=rev(seq(from=0,to=1,by=.5))
for (i in 1:length(beta.abunT))
{params=rbind(params,expand.grid(Simul="S1",rep=1:nrep,N_hosts=N_hosts,N_symbionts=N_symbionts,Ntraits=Ntraits,Breath=Breath,beta.comp=beta.comp,beta.abun=beta.abunT[i],beta.env=beta.envT[i],years=years,K=K,V_delta_host=V_delta_host,V_delta_symbiont=V_delta_symbiont))}
params=params[-1,]

#SM: we fix beta.env=0.5; beta.abun=0.5
beta.env=0.5
beta.abun=0.5

# S2. Varying years
years=c(20,40,80)
params=rbind(params,expand.grid(Simul="S6",rep=1:nrep,N_hosts=N_hosts,N_symbionts=N_symbionts,Ntraits=Ntraits,Breath=Breath,beta.comp=beta.comp,beta.abun=beta.abun,beta.env=beta.env,years=years,K=K,V_delta_host=V_delta_host,V_delta_symbiont=V_delta_symbiont))

# S3: Varying host numbers
years=20
N_hosts=10
params=rbind(params,expand.grid(Simul="S4",rep=1:nrep,N_hosts=N_hosts,N_symbionts=N_symbionts,Ntraits=Ntraits,Breath=Breath,beta.comp=beta.comp,beta.abun=beta.abun,beta.env=beta.env,years=years,K=K,V_delta_host=V_delta_host,V_delta_symbiont=V_delta_symbiont))

#S4: Varying symbiont numbers
N_hosts=25
N_symbionts=75
params=rbind(params,expand.grid(Simul="S5",rep=1:nrep,N_hosts=N_hosts,N_symbionts=N_symbionts,Ntraits=Ntraits,Breath=Breath,beta.comp=beta.comp,beta.abun=beta.abun,beta.env=beta.env,years=years,K=K,V_delta_host=V_delta_host,V_delta_symbiont=V_delta_symbiont))

# S5. varying individuals 
N_symbionts=150
K=300
params=rbind(params,expand.grid(Simul="S3",rep=1:nrep,N_hosts=N_hosts,N_symbionts=N_symbionts,Ntraits=Ntraits,Breath=Breath,beta.comp=beta.comp,beta.abun=beta.abun,beta.env=beta.env,years=years,K=K,V_delta_host=V_delta_host,V_delta_symbiont=V_delta_symbiont))

# S6. Varying delta symbiont
V_delta_symbiont=c(.01,.1,1,10,1000)
params=rbind(params,expand.grid(Simul="S7",rep=1:nrep,N_hosts=N_hosts,N_symbionts=N_symbionts,Ntraits=Ntraits,Breath=Breath,beta.comp=beta.comp,beta.abun=beta.abun,beta.env=beta.env,years=years,K=K,V_delta_host=V_delta_host,V_delta_symbiont=V_delta_symbiont))

# S7. Asymetric Competition 
V_delta_symbiont=1
beta.env=.5
beta.abun=.25
beta.comp=.25
params=rbind(params,expand.grid(Simul="2C",rep=1:nrep,N_hosts=N_hosts,N_symbionts=N_symbionts,Ntraits=Ntraits,Breath=Breath,beta.comp=beta.comp,beta.abun=beta.abun,beta.env=beta.env,years=years,K=K,V_delta_host=V_delta_host,V_delta_symbiont=V_delta_symbiont))

# S8. Multiple filter trtaits 
beta.env=.5
beta.abun=.5
beta.comp=0
Ntraits=c(1,2,4)
params=rbind(params,expand.grid(Simul="2B and 3B",rep=1:nrep,N_hosts=N_hosts,N_symbionts=N_symbionts,Ntraits=Ntraits,Breath=Breath,beta.comp=beta.comp,beta.abun=beta.abun,beta.env=beta.env,years=years,K=K,V_delta_host=V_delta_host,V_delta_symbiont=V_delta_symbiont))

# Give them an ID
params$IDparams=paste("No",1:dim(params)[1],sep="")
rownames(params)=params$IDparams


#---------------------------#
#    Run the simulations    #
#---------------------------#

NoCores=14 #depends how many cores you can use
Res=mclapply(ListParams,Sim,params=params,mc.cores=NoCores)


