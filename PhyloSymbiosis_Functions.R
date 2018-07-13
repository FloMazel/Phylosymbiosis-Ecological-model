


Sim=function(i,params)
{
  print(paste("Progression:",round(i/dim(params)[1],2)*100,"%",sep=" "))
  l=try(Simuls(paramVirt=params[i,],NoCores=1))
  if (!class(l)=="try-error"){Res= l[[1]];SR=l[[2]]} 
  if (class(l)=="try-error"){Res=SR="Failed Simulation"}
  return(list(Res,SR))
}

range02 <- function(x,newMax,newMin){ (x - min(x))/(max(x)-min(x)) * (newMax - newMin) + newMin } # to make comaprable host trait and microbe preferences

shuffleTips=function(x,tree)
{
  tree$tip.label=sample(tree$tip.label)
  return(tree)
}

RFmeasures=function(HostTree,SymbiontDendro,nRandom=100)
{
  ##obs=RF.dist(HostTree,as.phylo(SymbiontDendro),normalize = T)#phangorn
  obs=multiRF(HostTree,as.phylo(SymbiontDendro))#phangorn
  
  #shuffling tips ("null model")
  nullSymbiontDendro=lapply(1:nRandom,shuffleTips,as.phylo(SymbiontDendro))
  class(nullSymbiontDendro)="multiPhylo"
  null=sapply(1:nRandom,FUN=function(x){multiRF(HostTree,nullSymbiontDendro[[x]])})
  
  #new (random) trees ("random model")
  tips=as.phylo(SymbiontDendro)$tip.label
  randomSymbiontDendro=rmtree(N=nRandom,n=length(tips),tip.label = tips)
  rd=sapply(1:nRandom,FUN=function(x){multiRF(HostTree,randomSymbiontDendro[[x]])})
  
  #null=sapply(1:nRandom,FUN=function(x){RF.dist(HostTree,nullSymbiontDendro[[x]],normalize = T)})

  pval=sum(obs>null)/(nRandom+1)
  pvalRd=sum(obs>rd)/(nRandom+1)
  
  res=matrix(NA,ncol=2,nrow=2,dimnames = list(c("stat","pval"),c("RF","RFrd")))
  res["stat",]=rep(obs,2)
  res["pval","RF"]=pval
  res["pval","RFrd"]=pvalRd
  
  return(res)
}

Phylosymbiosis=function(OTUtable,HostTree,SymbiontTree,TraitDist=NA)
{
  conds=list(c('stat','pval','BlomK_Host','BlomK_Microbes'),c("RF","RFrd","Mantel","MantelTraits"),c("Unifrac","W_Unifrac","Bray","Jaccard"))
  Res=array(NA,dim=c(4,4,4),dimnames = conds)
  
  #Compute Beta
  unifracs <- aperm(GUniFrac(OTUtable, SymbiontTree, alpha=c(0, 0.5, 1))$unifracs,c(3,1,2)) #Unifrac
  bray <-vegdist(x=OTUtable,method="bray") #vegan
  jac <-vegdist(x=OTUtable,method="jaccard") #vegan
  
  
  #Compute correlation
    # Original method
      # Hierchachical clustering of microbiomes
      Dendro_W_Uni=hclust(as.dist(unifracs['d_1',,]),method = "average")
      Dendro_UW_Uni=hclust(as.dist(unifracs['d_UW',,]),method = "average")
      Dendro_Bray=hclust(bray,method = "average")
      Dendro_Jac=hclust(jac,method = "average")
      
      # Tree match
      Res[c("stat","pval"),c("RF","RFrd"),"W_Unifrac"] =RFmeasures(HostTree,SymbiontDendro=Dendro_W_Uni,nRandom=1000)
      Res[c("stat","pval"),c("RF","RFrd"),"Unifrac"] =RFmeasures(HostTree,SymbiontDendro=Dendro_UW_Uni,nRandom=1000)
      Res[c("stat","pval"),c("RF","RFrd"),"Bray"] =RFmeasures(HostTree,SymbiontDendro=Dendro_Bray,nRandom=1000)
      Res[c("stat","pval"),c("RF","RFrd"),"Jaccard"] =RFmeasures(HostTree,SymbiontDendro=Dendro_Jac,nRandom=1000)
      
      # Mantel
      Mantel_W_Uni=mantel(unifracs['d_1',,],cophenetic(HostTree))
      Res[c("stat","pval"),"Mantel","W_Unifrac"]=c(Mantel=Mantel_W_Uni$statistic,pval=Mantel_W_Uni$signif)
      Mantel_UW_Uni=mantel(unifracs['d_UW',,],cophenetic(HostTree))
      Res[c("stat","pval"),"Mantel","Unifrac"]=c(Mantel=Mantel_UW_Uni$statistic,pval=Mantel_UW_Uni$signif)
      Mantel_Bray=mantel(bray,cophenetic(HostTree))
      Res[c("stat","pval"),"Mantel","Bray"]=c(Mantel=Mantel_Bray$statistic,pval=Mantel_Bray$signif)
      Mantel_Jac=mantel(jac,cophenetic(HostTree))
      Res[c("stat","pval"),"Mantel","Jaccard"]=c(Mantel=Mantel_Jac$statistic,pval=Mantel_Jac$signif)
      
      
      
      # Mantel traits
      Mantel_W_Uni=mantel(unifracs['d_1',,],TraitDist[dimnames(unifracs)[[2]],dimnames(unifracs)[[2]]])
      Res[c("stat","pval"),"MantelTraits","W_Unifrac"]=c(Mantel=Mantel_W_Uni$statistic,pval=Mantel_W_Uni$signif)
      Mantel_UW_Uni=mantel(unifracs['d_UW',,],TraitDist[dimnames(unifracs)[[2]],dimnames(unifracs)[[2]]])
      Res[c("stat","pval"),"MantelTraits","Unifrac"]=c(Mantel=Mantel_UW_Uni$statistic,pval=Mantel_UW_Uni$signif)
      Mantel_Bray=mantel(bray,TraitDist[dimnames(unifracs)[[2]],dimnames(unifracs)[[2]]])
      Res[c("stat","pval"),"MantelTraits","Bray"]=c(Mantel=Mantel_Bray$statistic,pval=Mantel_Bray$signif)
      Mantel_Jac=mantel(jac,TraitDist[dimnames(unifracs)[[2]],dimnames(unifracs)[[2]]])
      Res[c("stat","pval"),"MantelTraits","Jaccard"]=c(Mantel=Mantel_Jac$statistic,pval=Mantel_Jac$signif)
      
      
      
return(Res)   
}


rescaleTree=function(Tree)
{
  H=max(nodeHeights(Tree))
  Tree$edge.length=Tree$edge.length/H
  return(Tree)
}

phylosigM=function(tree, x)
{
  if (class(x)=="matrix")
  {
    K=c()
    for (i in 1:dim(x)[2]){K=c(K,phylosig(tree, x[,i]))}
  }
  else K=phylosig(tree, x)
  K=mean(K)
  return(K)
}

Simuls=function(paramVirt,NoCores=3)
{
  
  #print("---------------------------------------")
  #print(paste("PARAMETERS --","  V_delta_host=",V_delta_host,"  / V_delta_symbiont=",V_delta_symbiont,sep=""))
  #print(paste("Rep#",rep,"--Simulating species",sep=""))
  #print(paramVirt)
  #print("Create Pools")
  
  # Create Hosts
  HostTree=rescaleTree(sim.bdtree(n=paramVirt[['N_hosts']]))  # Simulate a host tree
  HostTrait=sim.char(nsim=paramVirt[['Ntraits']],par=1,phy=rescale(HostTree,model="delta",paramVirt[['V_delta_host']]))[,,1]  # Evolve Host Traits 
  if (paramVirt[['Ntraits']]>1){HostTrait=sim.char(nsim=paramVirt[['Ntraits']],par=1,phy=rescale(HostTree,model="delta",paramVirt[['V_delta_host']]))[,1,]}  # Evolve Host Traits 
  
  if (paramVirt[['Ntraits']]==1){hostsNames=names(HostTrait)}
  if (paramVirt[['Ntraits']]>1){hostsNames=rownames(HostTrait)}
  
  HostBlomK=phylosigM(HostTree,HostTrait)   #Measure Phylosignal
  
  # Create Symbionts
  SymbiontTree=rescaleTree(sim.bdtree(n=paramVirt[['N_symbionts']]))   # Simulate a symbiont tree
  SymbiontTrait=sim.char(nsim=paramVirt[['Ntraits']],par=1,phy=rescale(SymbiontTree,model="delta",paramVirt[['V_delta_symbiont']]))[,,1]  # Evolve Symbiont Traits 
  if (paramVirt[['Ntraits']]>1){SymbiontTrait=sim.char(nsim=paramVirt[['Ntraits']],par=1,phy=rescale(SymbiontTree,model="delta",paramVirt[['V_delta_symbiont']]))[,1,]}  # Evolve Symbiont Traits 
  
  if (paramVirt[['Ntraits']]==1){SymbiontTrait=range02(SymbiontTrait,newMax = max(HostTrait),newMin = min(HostTrait))   } # RESCALING to make sure that Host Traits and Symbionts matches?
  if (paramVirt[['Ntraits']]>1) {for (i in 1:paramVirt[['Ntraits']]){ SymbiontTrait[,i]=range02(SymbiontTrait[,i],newMax = max(HostTrait[,i]),newMin = min(HostTrait[,i]))  }}
  
  SymbiontBlomK=phylosigM(SymbiontTree,SymbiontTrait)   #Measure Phylosignal
  
  #hist( SymbiontTrait)
  #hist(  range02(SymbiontTrait,newMax = max(HostTrait),newMin = min(HostTrait)),col="blue",add=T)
  #hist( HostTrait,col="red",add=T)

  # Assemble Microbial communities

 # print(paste("Rep#",paramVirt[['rep']],"--Assembling microbiomes"))
  #print("Assembling microbiomes")
  
  if (paramVirt[['Ntraits']]==1){  lres <- mclapply(1:paramVirt[['N_hosts']], function(x) {tamaure(niche.breadth=paramVirt[['Breath']], niche.optima=SymbiontTrait, env=HostTrait[x],  beta.env=paramVirt[['beta.env']], beta.comp=paramVirt[['beta.comp']],beta.abun=paramVirt[['beta.abun']], years=paramVirt[['years']], K=paramVirt[['K']], plot=F)},mc.cores=NoCores)} # 
  if (paramVirt[['Ntraits']]>1){  lres <- mclapply(1:paramVirt[['N_hosts']], function(x) {tamaureM(Ntraits=paramVirt[['Ntraits']],niche.breadth=paramVirt[['Breath']], niche.optima=SymbiontTrait, env=HostTrait[x,],  beta.env=paramVirt[['beta.env']], beta.abun=paramVirt[['beta.abun']], years=paramVirt[['years']], K=paramVirt[['K']], plot=F)},mc.cores=NoCores)} # 
  
  com <- t(sapply(lres, function(x) x$abundances)) # Get site x abundance matrix
  rownames(com)=hostsNames
  
  TraitDist=as.matrix(dist(HostTrait))
  
  # Measure Phylosymbiosis signal
 # print("Measuring Phylosymbiosis")
  
  PP=Phylosymbiosis(OTUtable=com,HostTree=HostTree,SymbiontTree=SymbiontTree,TraitDist=TraitDist)
  PhySig=c(HostBlomK=HostBlomK,SymbiontBlomK=SymbiontBlomK)
  PP[3:4,,]=PhySig
  
  SR=apply(com>1,1,sum)
  
  return(list(PP,SR))    
  #print("---------------------------------------")
}





tamaure=function (niche.breadth = 5, niche.optima, env, beta.env = 0, 
          beta.comp = 0, beta.abun = 0, years = 20, K = 20, community.in = NA, 
          species.pool.abundance = NA, plot = FALSE, competition = "symmetric", 
          intra.sp.com = 1) 
{
  species.count <- length(niche.optima)
  if (is.na(species.pool.abundance[1])) 
    species.pool.abundance <- rep(1, species.count)
  species.niche.dist <- as.matrix(dist(niche.optima, diag = TRUE, 
                                       upper = TRUE))
  species.niche.overlap.sym <- outer(niche.optima, niche.optima, 
                                     function(x, y) {
                                       2 * pnorm(-abs((x - y))/2, mean = 0, sd = niche.breadth)
                                     })
  species.niche.overlap.asym <- outer(niche.optima, niche.optima, 
                                      function(x, y) {
                                        sign <- ifelse(x > y, 1, 0)
                                        overlap <- 2 * pnorm(-abs((x - y))/2, mean = 0, sd = niche.breadth)
                                        sign * overlap
                                      })
  species.comp.hierarchy <- outer(niche.optima, niche.optima, 
                                  function(x, y) {
                                    ifelse(x > y, 1, 0)
                                  })
  diag(species.niche.overlap.sym) <- intra.sp.com
  diag(species.niche.overlap.asym) <- intra.sp.com
  diag(species.comp.hierarchy) <- intra.sp.com
  out.community <- as.data.frame(matrix(NA, nrow = years + 
                                          1, ncol = K))
  if (plot == TRUE) {
    community.null <- sample(names(niche.optima), K, replace = TRUE)
    mpd.now <- mean(species.niche.dist[community.null, community.null])
  }
  if (is.na(community.in[1])) {
    community <- out.community[1, ] <- sample(names(niche.optima), 
                                              K, replace = TRUE)
  }
  else {
    community <- out.community[1, ] <- community.in
  }
  log.p.env <- dnorm(x = niche.optima, mean = env, sd = niche.breadth, 
                     log = TRUE) - dnorm(x = env, mean = env, sd = niche.breadth, 
                                         log = TRUE)
  for (year in 1:years) {
    for (i in 1:K) {
      abundance <- as.numeric(table(community)[names(niche.optima)])
      abundance <- ifelse(is.na(abundance), 0, abundance)
      if (competition == "symmetric") 
        p.comp <- 1 - colSums(species.niche.overlap.sym[community, 
                                                        ])/K
      if (competition == "asymmetric") 
        p.comp <- 1 - colSums(species.niche.overlap.asym[community, 
                                                         ])/K
      if (competition == "hierarchical") 
        p.comp <- 1 - colSums(species.comp.hierarchy[community, 
                                                     ])/K
      p.all <- exp(beta.env * log.p.env + beta.comp * log(p.comp) + 
                     log(species.pool.abundance + beta.abun * abundance))
      p.all <- ifelse(is.na(p.all), min(p.all, na.rm = TRUE), 
                      p.all)
      if (sd(p.all, na.rm = TRUE) == 0) 
        p.all = NULL
      out <- sample(seq(community), 1)
      community[out] <- sample(names(niche.optima), 1, 
                               prob = p.all)
    }
    out.community[year + 1, ] <- community
    if (plot == TRUE) 
      mpd.now <- c(mpd.now, mean(species.niche.dist[community, 
                                                    community]))
    if (plot == TRUE & year == years) {
      n.rand <- round(years/4) + 1
      if (is.null(p.all)) 
        p.all <- rep(0, species.count)
      filters = cbind(comp = round(p.comp, 2), env = round(exp(log.p.env), 
                                                           2), abun = round(abundance, 2), all = round(p.all, 
                                                                                                       2))
      mpd.now <- c(mpd.now, mean(species.niche.dist[community, 
                                                    community]))
      mypar <- par(mfcol = c(3, 2))
      index <- order(niche.optima)
      barplot(filters[, 1][index], xlab = "Species ID", 
              ylab = "p.comp")
      barplot(filters[, 2][index], xlab = "Species ID", 
              ylab = "p.env")
      barplot(filters[, 3][index], xlab = "Species ID", 
              ylab = "abundance")
      barplot(filters[, 4][index], xlab = "Species ID", 
              ylab = "p.all")
      hist(table(community), xlab = "Abundance", ylab = "Frequency", 
           main = "")
      mpd.null <- sapply(1:1000, function(x) {
        community.null <- sample(names(niche.optima), 
                                 K, replace = TRUE)
        mean(species.niche.dist[community.null, community.null])
      })
      low.border <- sort(mpd.null)[50]
      high.border <- sort(mpd.null)[950]
      plot((1:length(mpd.now)), mpd.now, type = "l", xlab = "Time", 
           ylab = "mpd", ylim = c(min(c(low.border, mpd.now)), 
                                  max(c(high.border, mpd.now))))
      abline(h = c(low.border, high.border), col = 3)
      par(mypar)
    }
  }
  abundance.vector <- as.numeric(table(community)[names(niche.optima)])
  abundance.vector <- ifelse(is.na(abundance.vector), 0, abundance.vector)
  names(abundance.vector) <- names(niche.optima)
  return(list(community = community, abundances = abundance.vector, 
              communities.over.time = out.community))
}

# Laure Gallien
# 1 March 2016

# Community assembly from a given species pool


tamaureM <- function(niche.breadth = 5, 
                    niche.optima, 
                    env,
                    beta.env = 1, 
                    beta.abun = 0, 
                    years = 20, 
                    K = 20, 
                    Ntraits=2,
                    plot = FALSE ) {
  # niche.breadth = 5 ; niche.optima1 ; niche.optima2 ;  env1=50 ; env2=20 ; beta.env1 = 0.1 ; beta.env2 = 0.1 ; beta.abun = 1 ; years = 10 ; K = 100 ; plot = T			
  
  # ------- 1. getting input -------
  species.count <- dim(niche.optima)[1]  # nb of species in the species pool
  # if (is.na(species.pool.abundance[1])) 
  # species.pool.abundance <- rep(1, species.count)
  
  # ------- 2. The loop -------
  out.community <- as.data.frame(matrix(NA, nrow = years + 1, ncol = K))  # to store community composition over time
  if (plot == TRUE) {
    community.null <- sample(rownames(niche.optima), K, replace = TRUE)
    species.niche.dist <- as.matrix(dist(cbind(niche.optima, niche.optima), 
                                         diag = TRUE, upper = TRUE))  # to calculate mpd
    mpd.now <- mean(species.niche.dist[community.null, community.null])
  }
  
  # initialization: random
  community <- out.community[1, ] <- sample(rownames(niche.optima), K, replace = TRUE)
  
  
  # environmental filter (fix through time)
  
  log.p.env=list()
  for (i in 1:Ntraits)
  {  log.p.env[[i]] <- dnorm(x = niche.optima[,i], mean = env[i], sd = niche.breadth, log = TRUE) - dnorm(x = env[i], mean = env[i], sd = niche.breadth, log = TRUE)}
  log.p.env=do.call(cbind,log.p.env)
  
  # loop - years
  for (year in 1:years) {
    
    # loop - asynchroneous updating within years
    for (i in 1:K) {
      abundance <- as.numeric(table(community)[rownames(niche.optima)])
      # abundance <- ifelse(is.na(abundance), 0, abundance)
      abundance <- ifelse(is.na(abundance), 1, abundance)
      
      # p.all <- exp(beta.env * log.p.env  + log(beta.abun * abundance))
      #print(sum(p.all>1))
      #p.all <- exp(beta.env1*log.p.env1 + beta.env2*log.p.env2  + log(beta.abun*abundance))

      p.all <- exp(beta.env/Ntraits*apply(log.p.env,1,sum) + log(beta.abun*abundance))
      
      
      p.all <- ifelse(is.na(p.all), min(p.all, na.rm = TRUE), p.all)
      if (sd(p.all, na.rm = TRUE) == 0) p.all = NULL
      
      out <- sample(seq(community), 1)  # Choose randomly one individual to exclude
      community[out] <- sample(rownames(niche.optima), 1, prob = p.all)  # Replace this excluded species by lottery competition
    }
    out.community[year + 1, ] <- community
    if (plot == TRUE) mpd.now <- c(mpd.now, mean(species.niche.dist[community, community]))
    
    # ------- 3. Plot (if requested) -------
    if (plot == TRUE & year == years) {
      n.rand <- round(years/4) + 1
      if (is.null(p.all)) p.all <- rep(0, species.count)
      filters = cbind(env1 = round(exp(log.p.env1), 2), env2 = round(exp(log.p.env2), 2), 
                      abun = round(abundance, 2), all = round(p.all, 2))
      mpd.now <- c(mpd.now, mean(species.niche.dist[community, community]))
      mypar <- par(mfcol = c(3, 2))
      index <- order(niche.optima1)
      index2 <- order(niche.optima2)
      barplot(filters[, 1][index], xlab = "Species ID", ylab = "p.env1")
      barplot(filters[, 2][index2], xlab = "Species ID", ylab = "p.env2")
      barplot(filters[, 3][index], xlab = "Species ID", ylab = "abundance")
      barplot(filters[, 4][index], xlab = "Species ID", ylab = "p.all")
      hist(table(community), xlab = "Abundance", ylab = "Frequency", main = "")
      mpd.null <- sapply(1:1000, function(x) {
        community.null <- sample(names(niche.optima1), K, replace = TRUE)
        mean(species.niche.dist[community.null, community.null])
      })
      low.border <- sort(mpd.null)[50]
      high.border <- sort(mpd.null)[950]
      plot((1:length(mpd.now)), mpd.now, type = "l", xlab = "Time", ylab = "mpd", ylim = c(min(c(low.border, mpd.now)), max(c(high.border, mpd.now))))
      abline(h = c(low.border, high.border), col = 3)
      par(mypar)
    }
  }
  
  # ------- 5. Prepare the output -------
  abundance.vector <- as.numeric(table(community)[rownames(niche.optima)])
  abundance.vector <- ifelse(is.na(abundance.vector), 0, abundance.vector)
  names(abundance.vector) <- rownames(niche.optima)
  return(list(community = community, abundances = abundance.vector, communities.over.time = out.community))
}

nodeHeights=function (tree) 
{
  if (hasArg(root.edge)) 
    root.edge <- list(...)$root.edge
  else root.edge <- FALSE
  if (root.edge) 
    ROOT <- if (!is.null(tree$root.edge)) 
      tree$root.edge
  else 0
  else ROOT <- 0
  if (!inherits(tree, "phylo")) 
    stop("tree should be an")
  if (attr(tree, "order") != "cladewise" || is.null(attr(tree, 
                                                         "order"))) 
    t <- reorder(tree)
  else t <- tree
  root <- length(t$tip.label) + 1
  X <- matrix(NA, nrow(t$edge), 2)
  for (i in 1:nrow(t$edge)) {
    if (t$edge[i, 1] == root) {
      X[i, 1] <- 0
      X[i, 2] <- t$edge.length[i]
    }
    else {
      X[i, 1] <- X[match(t$edge[i, 1], t$edge[, 2]), 2]
      X[i, 2] <- X[i, 1] + t$edge.length[i]
    }
  }
  if (attr(tree, "order") != "cladewise" || is.null(attr(tree, 
                                                         "order"))) 
    o <- apply(matrix(tree$edge[, 2]), 1, function(x, y) which(x == 
                                                                 y), y = t$edge[, 2])
  else o <- 1:nrow(t$edge)
  return(X[o, ] + ROOT)
}



phylosig=function (tree, x, method = "K", test = FALSE, nsim = 1000, se = NULL, 
          start = NULL, control = list()) 
{
  if (!inherits(tree, "phylo")) 
    stop("tree should be an object of cla")
  #x <- matchDatatoTree(tree, x, "x")
  #tree <- matchTreetoData(tree, x, "x")
  if (!is.null(se)) {
    se <- matchDatatoTree(tree, se, "se")
    tree <- matchTreetoData(tree, se, "se")
    me = TRUE
    M <- diag(se^2)
    rownames(M) <- colnames(M) <- names(se)
  }
  else me = FALSE
  if (!is.null(start) && !is.null(se)) {
    if (start[1] <= 0 || start[2] < 0 || start[2] > maxLambda(tree)) {
      message("some of the elements of 'start' are invalid, resetting to random")
      start <- NULL
    }
  }
  if (method == "K") {
    C <- vcv.phylo(tree)
    x <- x[rownames(C)]
    n <- nrow(C)
    if (!me) {
      invC <- solve(C)
      a <- sum(invC %*% x)/sum(invC)
      K <- (t(x - a) %*% (x - a)/(t(x - a) %*% invC %*% 
                                    (x - a)))/((sum(diag(C)) - n/sum(invC))/(n - 
                                                                               1))
      if (!test) 
        return(as.numeric(K))
      else {
        P = 0
        simX <- x
        for (i in 1:nsim) {
          a <- sum(invC %*% simX)/sum(invC)
          simK <- (t(simX - a) %*% (simX - a)/(t(simX - 
                                                   a) %*% invC %*% (simX - a)))/((sum(diag(C)) - 
                                                                                    n/sum(invC))/(n - 1))
          if (simK >= K) 
            P <- P + 1/nsim
          simX <- sample(simX)
        }
        return(list(K = as.numeric(K), P = P))
      }
    }
    else {
      likelihoodK <- function(theta, C, M, y) {
        Ce <- theta * C + M
        invCe <- solve(Ce)
        a <- as.numeric(sum(invCe %*% y)/sum(invCe))
        logL <- -t(y - a) %*% invCe %*% (y - a)/2 - n * 
          log(2 * pi)/2 - determinant(Ce, logarithm = TRUE)$modulus/2
        return(logL)
      }
      M <- M[rownames(C), colnames(C)]
      invC <- solve(C)
      maxSig2 <- as.numeric(t(x - as.numeric(sum(invC %*% 
                                                   x)/sum(invC))) %*% invC %*% (x - as.numeric(sum(invC %*% 
                                                                                                     x)/sum(invC)))/n)
      res <- optimize(f = likelihoodK, interval = c(0, 
                                                    maxSig2), y = x, C = C, M = M, maximum = TRUE)
      sig2 <- res$maximum * n/(n - 1)
      Ce <- sig2 * C + M
      invCe <- solve(Ce)
      a <- as.numeric(sum(invCe %*% x)/sum(invCe))
      K <- (t(x - a) %*% (x - a)/(t(x - a) %*% invCe %*% 
                                    (x - a)))/((sum(diag(Ce)) - n/sum(invCe))/(n - 
                                                                                 1))
      if (!test) 
        return(list(K = as.numeric(K), sig2 = as.numeric(sig2), 
                    logL = res$objective[1, 1]))
      else {
        P = 0
        simX <- x
        for (i in 1:nsim) {
          maxSig2 <- as.numeric(t(simX - as.numeric(sum(invC %*% 
                                                          simX)/sum(invC))) %*% invC %*% (simX - as.numeric(sum(invC %*% 
                                                                                                                  simX)/sum(invC)))/n)
          simRes <- optimize(f = likelihoodK, interval = c(0, 
                                                           maxSig2), y = simX, C = C, M = M, maximum = TRUE)
          simSig2 <- simRes$maximum * n/(n - 1)
          Ce <- simSig2 * C + M
          invCe <- solve(Ce)
          a <- as.numeric(sum(invCe %*% simX)/sum(invCe))
          simK <- (t(simX - a) %*% (simX - a)/(t(simX - 
                                                   a) %*% invCe %*% (simX - a)))/((sum(diag(Ce)) - 
                                                                                     n/sum(invCe))/(n - 1))
          if (simK >= K) 
            P <- P + 1/nsim
          o <- sample(1:n)
          simX <- x[o]
          M <- diag(se[o]^2)
        }
        return(list(K = as.numeric(K), P = P, sig2 = as.numeric(sig2), 
                    logL = res$objective[1, 1]))
      }
    }
  }
  else if (method == "lambda") {
    lambda.transform <- function(C, lambda) {
      dC <- diag(diag(C))
      C <- lambda * (C - dC) + dC
      return(C)
    }
    likelihoodLambda <- function(theta, C, y) {
      Cl <- lambda.transform(C, theta)
      invCl <- solve(Cl)
      n <- nrow(Cl)
      y <- y[rownames(Cl)]
      a <- as.numeric(sum(invCl %*% y)/sum(invCl))
      sig2 <- as.numeric(t(y - a) %*% invCl %*% (y - a)/n)
      logL <- -t(y - a) %*% (1/sig2 * invCl) %*% (y - a)/2 - 
        n * log(2 * pi)/2 - determinant(sig2 * Cl, logarithm = TRUE)$modulus/2
      return(logL)
    }
    likelihoodLambda.me <- function(theta, C, y, M) {
      Cl <- theta[1] * lambda.transform(C, theta[2])
      V <- Cl + M
      invV <- solve(V)
      n <- nrow(Cl)
      y <- y[rownames(Cl)]
      a <- as.numeric(sum(invV %*% y)/sum(invV))
      logL <- -t(y - a) %*% invV %*% (y - a)/2 - n * log(2 * 
                                                           pi)/2 - determinant(V, logarithm = TRUE)$modulus/2
      return(-logL)
    }
    C <- vcv.phylo(tree)
    x <- x[rownames(C)]
    maxlam <- maxLambda(tree)
    if (!me) {
      res <- optimize(f = likelihoodLambda, interval = c(0, 
                                                         maxlam), y = x, C = C, maximum = TRUE)
      if (!test) 
        return(list(lambda = res$maximum, logL = res$objective[1, 
                                                               1]))
      else {
        logL0 <- likelihoodLambda(theta = 0, C = C, y = x)
        P <- as.numeric(pchisq(2 * (res$objective[1, 
                                                  1] - logL0), df = 1, lower.tail = FALSE))
        return(list(lambda = res$maximum, logL = res$objective[1, 
                                                               1], logL0 = logL0[1, 1], P = P))
      }
    }
    else {
      M <- M[rownames(C), colnames(C)]
      if (is.null(start)) 
        s <- c(0.02 * runif(n = 1) * mean(pic(x, multi2di(tree))^2), 
               runif(n = 1))
      else s <- start
      res <- optim(s, likelihoodLambda.me, C = C, y = x, 
                   M = M, method = "L-BFGS-B", lower = c(0, 0), 
                   upper = c(Inf, maxlam), control = control)
      if (!test) 
        return(list(lambda = res$par[2], sig2 = res$par[1], 
                    logL = -res$value, convergence = res$convergence, 
                    message = res$message))
      else {
        res0 <- optim(c(s[1], 0), likelihoodLambda.me, 
                      C = C, y = x, M = M, method = "L-BFGS-B", lower = c(0, 
                                                                          0), upper = c(Inf, 1e-10), control = control)
        P <- as.numeric(pchisq(2 * (res0$value - res$value), 
                               df = 1, lower.tail = FALSE))
        return(list(lambda = res$par[2], sig2 = res$par[1], 
                    logL = -res$value, convergence = res$convergence, 
                    message = res$message, logL0 = -res0$value, 
                    P = P))
      }
    }
  }
  else stop(paste("do not recognize method = ", 
                  sep = ""))
}

GUniFrac=function (otu.tab, tree, alpha = c(0, 0.5, 1)) 
{
  if (!is.rooted(tree)) 
    stop("Rooted phylogenetic tree required!")
  otu.tab <- as.matrix(otu.tab)
  row.sum <- rowSums(otu.tab)
  otu.tab <- otu.tab/row.sum
  n <- nrow(otu.tab)
  if (is.null(rownames(otu.tab))) {
    rownames(otu.tab) <- paste("comm", 1:n, sep = "_")
  }
  dimname3 <- c(paste("d", alpha, sep = "_"), "d_UW", "d_VAW")
  unifracs <- array(NA, c(n, n, length(alpha) + 2), dimnames = list(rownames(otu.tab), 
                                                                    rownames(otu.tab), dimname3))
  for (i in 1:(length(alpha) + 2)) {
    for (j in 1:n) {
      unifracs[j, j, i] <- 0
    }
  }
  if (sum(!(colnames(otu.tab) %in% tree$tip.label)) != 0) {
    stop("The OTU table contains unknown OTUs! OTU names\\n\\t\\t\\t\\t\\tin the OTU table and the tree should match!")
  }
  absent <- tree$tip.label[!(tree$tip.label %in% colnames(otu.tab))]
  if (length(absent) != 0) {
    tree <- drop.tip(tree, absent)
    warning("The tree has more OTU than the OTU table!")
  }
  tip.label <- tree$tip.label
  otu.tab <- otu.tab[, tip.label]
  ntip <- length(tip.label)
  nbr <- nrow(tree$edge)
  edge <- tree$edge
  edge2 <- edge[, 2]
  br.len <- tree$edge.length
  cum <- matrix(0, nbr, n)
  for (i in 1:ntip) {
    tip.loc <- which(edge2 == i)
    cum[tip.loc, ] <- cum[tip.loc, ] + otu.tab[, i]
    node <- edge[tip.loc, 1]
    node.loc <- which(edge2 == node)
    while (length(node.loc)) {
      cum[node.loc, ] <- cum[node.loc, ] + otu.tab[, i]
      node <- edge[node.loc, 1]
      node.loc <- which(edge2 == node)
    }
  }
  cum.ct <- round(t(t(cum) * row.sum))
  for (i in 2:n) {
    for (j in 1:(i - 1)) {
      cum1 <- cum[, i]
      cum2 <- cum[, j]
      ind <- (cum1 + cum2) != 0
      cum1 <- cum1[ind]
      cum2 <- cum2[ind]
      br.len2 <- br.len[ind]
      mi <- cum.ct[ind, i] + cum.ct[ind, j]
      mt <- row.sum[i] + row.sum[j]
      diff <- abs(cum1 - cum2)/(cum1 + cum2)
      for (k in 1:length(alpha)) {
        w <- br.len2 * (cum1 + cum2)^alpha[k]
        unifracs[i, j, k] <- unifracs[j, i, k] <- sum(diff * 
                                                        w)/sum(w)
      }
      ind2 <- (mt != mi)
      w <- br.len2 * (cum1 + cum2)/sqrt(mi * (mt - mi))
      unifracs[i, j, (k + 2)] <- unifracs[j, i, (k + 2)] <- sum(diff[ind2] * 
                                                                  w[ind2])/sum(w[ind2])
      cum1 <- (cum1 != 0)
      cum2 <- (cum2 != 0)
      unifracs[i, j, (k + 1)] <- unifracs[j, i, (k + 1)] <- sum(abs(cum1 - 
                                                                      cum2)/(cum1 + cum2) * br.len2)/sum(br.len2)
    }
  }
  return(list(unifracs = unifracs))
}


multiRF<-function(tree1,tree2){
  trees=list(tree1,tree2)
  class(trees)="multiPhylo"
  if(class(trees)!="multiPhylo")
    stop("trees should be an object of class \"multiPhylo\"")
  N<-length(trees)
  RF<-matrix(0,N,N)
  if(any(sapply(unclass(trees),is.rooted))){
    #cat("Some trees are rooted. Unrooting all trees.\n")
    trees<-lapply(unclass(trees),unroot)
  }
  foo<-function(pp) lapply(pp,function(x,pp)
    sort(attr(pp,"labels")[x]),pp=pp)
  xx<-lapply(unclass(trees),function(x) foo(prop.part(x)))
  for(i in 1:(N-1)) for(j in (i+1):N)
    RF[i,j]<-RF[j,i]<-2*sum(!xx[[i]]%in%xx[[j]])
  RF=RF[1,2]
  RF=RF / (Nnode(tree1) + Nnode(tree2) - 2)
  RF
}
