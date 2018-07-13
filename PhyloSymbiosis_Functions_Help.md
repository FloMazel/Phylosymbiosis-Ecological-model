### Function _Sim_

**Description:** 

**Input:** The function take a set of input simulation parameters 

**Output:** 

**Parameters:**  

  - i: refers to the i-th line of the "params" (see below) matrix that is used to run simulations 

  - params: matrix that contains the simulation parameters; the matrix column names should contains:
     - rep: the number of simulation you want to run    
     - N_hosts: number of host species (only one microbiome by host)          
     - N_symbionts: number of microbial units (e.g. "species", OTUs, ESVs)
     - Ntraits: number of host traits used to fiulter microbiomes   

     - V_delta_host: Delta parameter (sensu Pagel 1999) used to simulate host traits along host phylogeny 
     - V_delta_symbiont: Delta parameter (sensu Pagel 1999) used to simulate microbial traits along microbial phylogeny 



###### Breath: value of standard deviation of the Gaussian distributions that describe the microbial niches (identical for all species)
###### beta.comp: value of the strength of the competition filter (see details)
###### beta.abun: value of the strength of the competition filter (see details)
###### beta.env: value of the strength of the environmental filter (see details)
###### years: number of simulated time-steps;
###### K: value of carrying capacity, i.e. number of microbe individuals in the microbial community                


**Details:**  

  -- *Community assembly*. Community assembly is simulated by asynchroneous updating of K individuals per year. For each update a random individual is removed and replaced by an individual from the species pool. The choice of that individual follows a multinomial distribution, with probabilities being driven by an environmental filter, a competition filter and a recruitment filter.

  -- *Model of trait evolution*. delta is a time-dependent model of trait evolution (Pagel 1999). The delta model is similar to ACDC insofar as the delta model fits the relative contributions of early versus late evolution in the tree to the covariance of species trait values. Where delta is greater than 1, recent evolution has been relatively fast; if delta is less than 1, recent evolution has been comparatively slow. Intrepreted as a tree transformation, the model raises all node depths to an estimated power (delta). 


Münkemüller T, Gallien L. 2015. VirtualCom: a simulation model for eco-evolutionary community assembly and invasion. Methods Ecol Evol 6:735–743.
Pagel M. 1999. Inferring the historical patterns of biological evolution. Nature 401:877-884


### Function _RFmeasures_

**Description:** 

**Input:** The function take a set of input simulation parameters 

**Output:** 

**Parameters:**  

### Function _Phylosymbiosis_

**Description:** 

**Input:** The function take a set of input simulation parameters 

**Output:** 

**Parameters:**  




