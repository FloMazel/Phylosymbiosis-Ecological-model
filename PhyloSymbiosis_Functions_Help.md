### Function _Sim_

**Description:** The function use simulation parameters input to (1) assemble microbiome across host species and (2) measure phylosymbiosis patterns using various approaches 

**Input parameters:**  

  - i: refers to the i-th line of the "params" (see below) matrix that is used to run simulations 

  - params: matrix that contains the simulation parameters; the matrix column names should contains:
     - rep: the number of simulation you want to run    
     - N_hosts: number of host species (only one microbiome by host)          
     - N_symbionts: number of microbial units (e.g. "species", OTUs, ESVs)
     - Ntraits: number of host traits used to fiulter microbiomes   
     - V_delta_host: Delta parameter (sensu Pagel 1999) used to simulate host traits along host phylogeny 
     - V_delta_symbiont: Delta parameter (sensu Pagel 1999) used to simulate microbial traits along microbial phylogeny 
     - Breath: value of standard deviation of the Gaussian distributions that describe the microbial niches (identical for all species)
     - beta.comp: value of the strength of the competition filter (see details)
     - beta.abun: value of the strength of the competition filter (see details)
     - beta.env: value of the strength of the environmental filter (see details)
     - years: number of simulated time-steps;
     - K: value of carrying capacity, i.e. number of microbe individuals in the microbial community                

**Output:** The function returns a list of results corresponding to each line of the "param" input paramter. For each of these independant simulations, the function returns a list with thwo elements  outputs:  

  - First element: an array summarizing phylosymbiosis statistics (statistic and associated [p-value]) for four diversity metrics (Unifrac, Weighted Unifrac, Bray Curtis and Jaccard) and four approaches: 
       - _RF_: Dendogram-nased method (RF metric compared to model null 1, where tip labels are shuffled)
       - _RFrd_: Dendogram-nased method (RF metric compared to model null 2. where tree are re-simulated)
       - _Mantel_: Mantel method (Mantel test between host phylogenetic distances and microbial beta-diversity; statistic is Speraman Rho)
       - _ MantelTraits_: Mantel method (Mantel test between host trait distances and microbial beta-diversity; statistic is Speraman Rho)

  - Second element: a vector depicting the microbiome richness (number of microbial units present in each host)

**Details:**  

  -- *Community assembly*. Community assembly is simulated by asynchroneous updating of K individuals per year. For each update a random individual is removed and replaced by an individual from the species pool. The choice of that individual follows a multinomial distribution, with probabilities being driven by an environmental filter, a competition filter and a recruitment filter.

  -- *Model of trait evolution*. delta is a time-dependent model of trait evolution (Pagel 1999). The delta model is similar to ACDC insofar as the delta model fits the relative contributions of early versus late evolution in the tree to the covariance of species trait values. Where delta is greater than 1, recent evolution has been relatively fast; if delta is less than 1, recent evolution has been comparatively slow. Intrepreted as a tree transformation, the model raises all node depths to an estimated power (delta). 


Münkemüller T, Gallien L. 2015. VirtualCom: a simulation model for eco-evolutionary community assembly and invasion. Methods Ecol Evol 6:735–743.

Pagel M. 1999. Inferring the historical patterns of biological evolution. Nature 401:877-884


### Function _RFmeasures_

**Description:** The functions compares the observed host tree and abserved microbial community dendrogram with the Robinson Fould (RF) metric, and contrast the value of the pbserved RF metric to null expectations

**Input parameters:** 

  - HostTree: Host Phylogenetic tree

  - SymbiontDendro: microbial community dendrogram
  
  - nRandom: number of randomizations (default=100)

**Output:**  A table RF statistics (statistic and associated [p-value]) for two approaches: 

  - RF: Dendogram-nased method (RF metric compared to model null 1, where tip labels are shuffled)
  
  - RFrd: Dendogram-nased method (RF metric compared to model null 2. where tree are re-simulated)

### Function _Phylosymbiosis_

**Description:** The function compute phylosymbiosis statistics from an OTU table, an host phylogenetic tree and a set of hos traits (optional)

**Input parameters:**

  - OTUtable: an host species * microbial unit table 
  
  - HostTree: an host phylogenetic tree 
  
  - SymbiontTree: a symbiont phylogenetic tree (for phylogenetic beta-diversity mnetrics such as Unifrac)
  
  - TraitDist: ahost trait distance matrix 


**Output:** An array summarizing phylosymbiosis statistics (statistic and associated [p-value]) for four diversity metrics (Unifrac, Weighted Unifrac, Bray Curtis and Jaccard) and four approaches: 

       - _RF_: Dendogram-nased method (RF metric compared to model null 1, where tip labels are shuffled)
       
       - _RFrd_: Dendogram-nased method (RF metric compared to model null 2. where tree are re-simulated)
       
       - _Mantel_: Mantel method (Mantel test between host phylogenetic distances and microbial beta-diversity; statistic is Speraman Rho)
       
       - _ MantelTraits_: Mantel method (Mantel test between host trait distances and microbial beta-diversity; statistic is Speraman Rho)