### Workflow -----
### for data analysis and prediction of Jared Zystro's se sweet corn breeding project

# Style guide -----
# ALLCAPS - global constants
# UpperCamel - functions
# lower_underscore - variables
# Space between operators
# No space after function call

# PART A - for all traits -----

# 1. Initialize -----
#     Clean up
#     Set global variables
InitializeProject <- function() {
  
  rm(list = ls(envir = .GlobalEnv), envir = .GlobalEnv) 
  
  WORKING_DIRECTORY <<- "/home/sweetcorn/sweetcornbreeding"
  # PHENO_DATA_FILE <<- "~/Dropbox/PhD/2015 Trials/Data.csv"
  PHENO_DATA_FILE <<- "./data/Data.csv"
  MARKER_DATA_FILE <<- "./data/data_parimp_37genos_dropminorhets_dropNAs_thinto250KB.hmp.txt"
  FUNCTIONS_SUBDIRECTORY <<- "./functions"
  NUM_GENOTYPES <<- 38
  TRAITS <<- c("PlantHt","EarHt","Flavor","Tenderness","Husk","Tip","Width","Length","Pollen")
  
  setwd(WORKING_DIRECTORY)
  
  if (!require("lme4")) install.packages("lme4")
  library (lme4)
  if (!require("lmerTest")) install.packages("lmerTest")
  library (lmerTest)
  if (!require("reshape")) install.packages("reshape")
  library (reshape)
  if (!require("Hmisc")) install.packages("Hmisc")
  library (Hmisc)
  if (!require("estimability")) install.packages("estimability")
  library (estimability)
  if (!require("parallel")) install.packages("parallel")
  library(parallel)
  if (!require("doSNOW")) install.packages("doSNOW")
  library(doSNOW)
  if (!require("bigmemory")) install.packages("bigmemory")
  library(bigmemory)
  if (!require("tidyr")) install.packages("tidyr")
  library(tidyr)
  
}
InitializeProject()

# 2. Import / create data tables -----
#     Load in 2015 and 2016 phenotype data
LoadPhenoData <- function() {
  
  return (read.csv(PHENO_DATA_FILE))
  
}
pheno_data <- LoadPhenoData()

#     Load in inbred marker data
LoadMarkerData <- function () {
  
  return (read.delim(MARKER_DATA_FILE))
  
}
marker_data <- LoadMarkerData()

# 3. Format data -----
#     Set factors in phenotype data
AddPhenoDataFactors <- function (pheno_data) {
  
pheno_data$Geno <- factor(pheno_data$Geno)
pheno_data$Block <- factor(pheno_data$Block)
pheno_data$Set <- factor(pheno_data$DIISet)
pheno_data$MaleGeno <- factor(pheno_data$MaleGeno)
pheno_data$FemaleGeno <- factor(pheno_data$FemaleGeno)
pheno_data$Geno <- factor(pheno_data$Geno)
pheno_data$Rep <- factor(pheno_data$Rep)
pheno_data$Row <- factor(pheno_data$Row)
pheno_data$Col <- factor(pheno_data$Col)
pheno_data$Year <- factor(pheno_data$Year)

pheno_data$Env <- factor(paste0(pheno_data$Location,pheno_data$Year))
pheno_data$EnvRep <- factor(paste0(pheno_data$Env,pheno_data$Rep))
pheno_data$EnvBlock <- factor(paste0(pheno_data$Env,pheno_data$Block))


# Trying another way to see if this has been messing up the model
# pheno_data$Env <- with (pheno_data,Location:Year)
# pheno_data$EnvRep <- with(pheno_data,Rep:Env)
# pheno_data$EnvBlock <-with(pheno_data,Block:Env)

return (pheno_data)

}
pheno_data <- AddPhenoDataFactors(pheno_data)

#     Create phenotypic data subsets
CreateDIISubset <- function (pheno_data) {
  
  return (droplevels(subset(pheno_data, Entry>40 & Entry<=140)))

}
hybrid_pheno <- CreateDIISubset(pheno_data)

CreateInbredSubset <- function (pheno_data) {
  
  return (droplevels((subset(pheno_data, Entry<=40))))
  
}
inbred_pheno <- CreateInbredSubset(pheno_data)

CreateCheckSubset <- function (pheno_data) {
  
  return (droplevels((subset(pheno_data, Entry>140))))
  
}
check_pheno <- CreateCheckSubset(pheno_data)

CreateSets <- function (hybrid_pheno) {
  
  dii_sets <- list()
  
  for (i in 1:4) {
    dii_sets[[i]] <- droplevels(subset(hybrid_pheno, Set==i))
  }
  
  return (dii_sets)  
}
dii_sets <- CreateSets(hybrid_pheno)

#     Format raw marker data
FormatMarkerData <- function (marker_data) {
  
  # Convert from HapMap format into 0,1,2 format
  # Code from the NAM package vignette: 
  # ftp://ftp.math.ethz.ch/sfs/CRAN/web/packages/NAM/vignettes/vignette.html
  G <- marker_data
  AA <- as.character(G[, 2]); 
  AA <- gsub('/', '', AA); 
  Gen <- t(G[, -c(1:11)]); 
  n <- nrow(Gen); 
  m <- ncol(Gen); 
  gen <- matrix(NA, n, m); 
  for (i in 1:m) { 
    A1 <- strsplit(AA[i], '') [[1]][1]; 
    A2 <- strsplit(AA[i], '')[[1]][2];
    BB <- paste(A1, A1, sep=''); 
    Bb <- paste(A1, A2, sep=''); 
    bb <- paste(A2, A2, sep=''); 
    M <- as.character(Gen[, i]); 
    gen[M == BB, i] <- 2; 
    gen[M == Bb, i] <- 1;
    gen[M == bb, i] <- 0}; 
  colnames(gen) <- paste(G[, 3],G[, 4],G[, 2],sep = '.')
  gen2 <- gen
  
  # Clean up names
  rownames(gen2) <- gsub("(.*?)\\.(.*)","\\1", rownames(Gen))
  rownames(gen2)[rownames(gen2) == "We1095047"] <- "We1095407-4"
  rownames(gen2)[rownames(gen2) == "We1095407"] <- "We1095407-2"
  rownames(gen2)[rownames(gen2) == "We09423"] <- "We09423a"
  rownames(gen2)[rownames(gen2) == "We07812"] <- "Wu07812"
  rownames(gen2)[rownames(gen2) == "We10802R"] <- "Wu10802R"
  rownames(gen2)[rownames(gen2) == "We10803R"] <- "Wu10803R"
  
  marker_data <- gen2
  
  return (marker_data)
  
}
marker_data <- FormatMarkerData(marker_data)

#     Create pedigree table
PedCombn <- function (parents) {
  
  return (c(paste0(parents[1], "x", parents[2]), paste0(parents[1]), paste0(parents[2])))
  
}
CreatePedTable <- function (parent_names) {
  
  inbredpeds <- data.frame(indiv=parent_names[,1],par1="NA",par2="NA")
  hybridpeds <- t(combn(as.vector(parent_names[1:NUM_GENOTYPES,1]), 2, FUN = PedCombn))
  hybridpeds <- data.frame(indiv=hybridpeds[,1],par1=hybridpeds[,2],par2=hybridpeds[,3])

  ped_table <- rbind(inbredpeds,hybridpeds)
  
  return (ped_table)
}
ped_table <- CreatePedTable(data.frame(rownames(marker_data)))

# Create marker data for hybrids from inbred marker data
CreateHybridMarkers <- function (marker_data, ped_table) {
  
  inbred_geno <- marker_data
  hybrid_peds <- ped_table[-(1:NUM_GENOTYPES), ]
  hybrid_geno <- matrix(NA, nrow = nrow(hybrid_peds), ncol = ncol(inbred_geno))
  
  # Average of parental alleles at each locus
  for (i in 1:nrow(hybrid_peds)) {
    hybrid_geno[i, ] <- (inbred_geno[rownames(inbred_geno) == hybrid_peds$par1[i]]+
                        inbred_geno[rownames(inbred_geno) == hybrid_peds$par2[i]])/2
  }
  
  rownames(hybrid_geno) <- hybrid_peds$indiv
  colnames(inbred_geno) <- 1:ncol(inbred_geno)
  colnames(hybrid_geno) <- 1:ncol(hybrid_geno)
  
  all_geno <- rbind(inbred_geno, hybrid_geno)
  
  return (all_geno)

}
marker_data <- CreateHybridMarkers(marker_data, ped_table)

# PART B - for each trait -----

ProtectedList <- function (list_name) {
  
  if (!exists(list_name)) {
    
    return_list <- list()
    
  } else {
    
    return_list <- get(list_name)
    
  }
  
  return (return_list)
  
}

SetTraitName <- function (trait_number) {

  return (TRAITS[trait_number])
  
}

for (trait_num in 1:length(TRAITS)) {

# trait_num <- 1

trait_name <- SetTraitName(trait_num) ### TO BE LOOPED?

# 4. Quality control ----- TO DO
#     Check for outliers
#     Test for normality
#     Test for equal variance

# ANOVA Lit ----
# From Suwarno, 2012 (Disseration): 

# Y = μ + Env + Rep(Env) + Block(Rep x Env) + Set + GCA1(Set) + GCA2(Set) + 
# SCA(Set) + Env x Set + Env x [GCA1(Set)] + Env x [GCA2(Set)] + Env x [SCA(Set)] + ε

# Where μ = grand mean, Env = environment, Rep = replication, GCA1 and GCA2 = 
# general combining ability of parent-1 and parent-2, respectively, SCA = specific
# combining ability, and ε = experimental error. Set, GCA, and SCA were considered 
# as fixed effects, whereas environment, replication, block, and all interactions 
# involving these factors were random. To understand the effects of hybrid and hybrid
# x environment interaction, the same model aggregating Set, GCA1(Set), GCA2(Set), 
# SCA(Set) as “hybrid” and their interaction with the environment as “hybrid x 
# environment” was fitted to the data.

# From lmer for SAS PROC MIXED Users:
# Multilocation$Grp <- with(Multilocation, Block:Location)
# (fm1Mult <- lmer(Adj ˜ Location * Trt + (1|Grp), Multilocation))

# 5. Test for experiment-wide GxE -----
#     ANOVA
CreateGxEModel<- function (trait_name, pheno_data) {

  formula_name <- as.formula(paste0(trait_name, " ~ Geno * Env + (1|EnvRep) + (1|EnvBlock)"))
  return (lmer(formula_name, data = pheno_data))  ### REQUIRES lme4
  
}
gxe_model <- ProtectedList("gxe_model")
gxe_model[[trait_name]] <- CreateGxEModel(trait_name, pheno_data)
gxe_anova <- ProtectedList("gxe_anova")
gxe_anova[[trait_name]] <- anova(gxe_model[[trait_name]]) ### REQUIRES lmerTest

#     Spearman rank correlation
GxESpearman <- function (trait_name, pheno_data) {
  
  trait_loc <- aggregate(as.formula(paste0(trait_name," ~ Geno:Env")), data = pheno_data, mean, na.rm = TRUE)
  trait_loc <- cast(trait_loc, Geno ~ Env, value = c(as.character(trait_name))) ### REQUIRES reshape
  
  return (rcorr(as.matrix(trait_loc), type = "spearman")) ### REQUIRES Hmisc
  
}
spear_corr <- ProtectedList ("spear_corr")
spear_corr[[trait_name]] <- GxESpearman(trait_name, pheno_data)

# 6. ANOVAs ----
#     Check ANOVAs
CreateCheckModel<- function (trait_name, check_pheno) {
  
  formula_name <- as.formula(paste0(trait_name, " ~ EnvBlock + (1|Geno)"))
  return (lmer(formula_name, data = check_pheno))
  
}
check_model <- ProtectedList("check_model")
check_model[[trait_name]] <- CreateCheckModel(trait_name, check_pheno)
check_anova <- ProtectedList("check_anova")
check_anova[[trait_name]] <- anova(check_model[[trait_name]])

#     Adjust means based on checks? TO DO

#     Inbred ANOVAs
CreateInbredModel<- function (trait_name, inbred_pheno) {
  
  formula_name <- as.formula(paste0(trait_name, " ~ Geno * Env + (1 | EnvRep)"))
  return (lmer(formula_name, data = inbred_pheno))
  
}
inbred_model <- ProtectedList("inbred_model")
inbred_model[[trait_name]] <- CreateInbredModel(trait_name, inbred_pheno)
inbred_anova <- ProtectedList("inbred_anova")
inbred_anova[[trait_name]] <- anova(inbred_model[[trait_name]])

#     Hybrid ANOVAs
CreateHybridModel<- function (trait_name, hybrid_pheno) {
  
  # Y = μ + Env + Rep(Env) + Block(Rep x Env) + Set + GCA1(Set) + GCA2(Set) + 
  # SCA(Set) + Env x Set + Env x [GCA1(Set)] + Env x [GCA2(Set)] + Env x [SCA(Set)] + ε
  
  
  formula_name <- as.formula(paste0(trait_name, " ~ Set + MaleGeno * FemaleGeno + (1 | EnvRep)"))
  return (lmer(formula_name, data = hybrid_pheno))
  
}
hybrid_model <- ProtectedList("hybrid_model")
hybrid_model[[trait_name]] <- CreateHybridModel(trait_name, hybrid_pheno)
hybrid_anova <- ProtectedList("hybrid_anova")
hybrid_anova[[trait_name]] <- anova(hybrid_model[[trait_name]])

# 7. Collect phenotype values -----
#     Inbred means

# NOTE: Inflated error b/c Env was dropped to address missing values
GetInbredMeans <- function (trait_name, inbred_pheno) {
  
  formula_name <- as.formula(paste0(trait_name, " ~ Geno + (1 | EnvRep)"))
  inbred_perse_model <- lmer(formula_name, data = inbred_pheno)
  inbred_lsmeans <- (lsmeans(inbred_perse_model, "Geno"))$lsmeans.table
  inbred_means <- data.frame (Genotype = inbred_lsmeans$Geno, 
                              Estimate = inbred_lsmeans$Estimate, 
                              StdErr = inbred_lsmeans$`Standard Error`)
  return (inbred_means)
  
}
inbred_means <- ProtectedList("inbred_means")
inbred_means[[trait_name]] <- GetInbredMeans(trait_name, inbred_pheno)

#     Tested hybrid means
#     Inbred GCAs

### Cannot seem to extract estimates from full model, splitting into parts
CreateSetModels <- function (trait_name, dii_sets) { ### DEBUG
  
  dii_model <- list()
  formula_name <- as.formula(paste0(trait_name, " ~ MaleGeno * FemaleGeno + (1 | EnvRep)"))
  
#  for (i in 1:3) { ### DEBUG
#    i <- 1 ### DEBUG
#    j <- as.numeric(i) ### DEBUG
#    print (j) ### DEBUG
#    ls() ### DEBUG
#    data_set_name <- get(paste0(as.character(paste0("dii_sets[", j, "]"))))
#    dii_model[[j]] <- lmer (formula_name, data = ) 
#  }
  
# Problem with with assigning list members using "i" in loop, dropped the loop
  dii_model[[1]] <- lmer (formula_name, data = dii_sets[[1]]) 
  dii_model[[2]] <- lmer (formula_name, data = dii_sets[[2]]) 
  dii_model[[3]] <- lmer (formula_name, data = dii_sets[[3]]) 
  dii_model[[4]] <- lmer (formula_name, data = dii_sets[[4]]) 
  
  return (dii_model) #DEBUG
  
} ### DEBUG
# debugonce(CreateSetModels) ### DEBUG
dii_model <- ProtectedList ("dii_model")
dii_model[[trait_name]] <- CreateSetModels(trait_name, dii_sets) ### DEBUG

GetGCAs <- function (dii_model) {
  
  gca <- list()
  
  for (i in 1:4) {
    
 #   i <- 1 ### DEBUG
    sca_lsmeans <- (lsmeans(dii_model[[i]], "MaleGeno:FemaleGeno"))$lsmeans.table ### REQUIRES estimability
    grandmean <- mean(sca_lsmeans$Estimate)
    
    fgca_lsmeans <- (lsmeans(dii_model[[i]], "FemaleGeno"))$lsmeans.table
    fgca <- data.frame (Genotype = fgca_lsmeans$FemaleGeno, 
                        GCA = fgca_lsmeans$Estimate - grandmean, 
                        StdErr = fgca_lsmeans$`Standard Error`)
    
    mgca_lsmeans <- (lsmeans(dii_model[[i]],"MaleGeno"))$lsmeans.table
    mgca <- data.frame (Genotype = mgca_lsmeans$MaleGeno, 
                        GCA = mgca_lsmeans$Estimate - grandmean, 
                        StdErr = mgca_lsmeans$`Standard Error`)
    
    gca[[i]] <- merge(fgca, mgca, all=TRUE)
  }
  
  return(do.call(rbind,gca))
  
}
gcas <- ProtectedList ("gcas") 
gcas[[trait_name]] <- GetGCAs(unlist(dii_model[[trait_name]]))

#     Tested hybrid SCAs
GetSCAs <- function (gcas, dii_model) {
  
  sca <- list()
  
  for (i in 1:4) {
#    i <- 1 ### DEBUG
    sca_lsmeans <- (lsmeans(dii_model[[i]], "MaleGeno:FemaleGeno"))$lsmeans.table
    grandmean <- mean(sca_lsmeans$Estimate)
    
    for (j in 1:nrow(sca_lsmeans)) {
#      j <- 1 ### DEBUG
      msca <- gcas[gcas$Genotype == paste0(sca_lsmeans[j, "MaleGeno"]), "GCA"]
      fsca <- gcas[gcas$Genotype == paste0(sca_lsmeans[j, "FemaleGeno"]), "GCA"]
      sca_lsmeans[j, "SCA"] <- sca_lsmeans[j, "Estimate"] - grandmean - msca - fsca
      sca_lsmeans[j, "Genotype"] <- paste0(sca_lsmeans[j, "FemaleGeno"], "x", sca_lsmeans[j, "MaleGeno"])
    }
    
    sca[[i]] <- data.frame (Genotype = sca_lsmeans$Genotype, 
                            Estimate = sca_lsmeans$Estimate, 
                            SCA = sca_lsmeans$SCA, 
                            StdErr = sca_lsmeans$`Standard Error`)
  }
  
  scas<-do.call(rbind,sca)
  
  return(scas)
}
# debugonce(GetSCAs) ### DEBUG
scas <- ProtectedList("scas")
scas[[trait_name]] <- GetSCAs(gcas[[trait_name]],unlist(dii_model[[trait_name]]))

# 8. Predict untested hybrid means -----
#     Classical (u + GCAf + GCAm)
##TODO
PredictMeansClassic <- function (ped_table, scas, gcas, inbred_means) {

  hybrid_mean <- mean(scas$Estimate)
  all_means <- data.frame(matrix(NA, nrow = nrow(ped_table), ncol = 2))
  
  for (i in 1:nrow(all_means)){
    
    ID <- paste0(ped_table$indiv[i])
    P1 <- paste0(ped_table$par1[i])
    P2 <- paste0(ped_table$par2[i])
    
    all_means[i, 1] <- ID
    
    if (!grepl("x", ID)) {
      if (length(inbred_means[inbred_means$Genotype == ID, "Estimate"]) == 0) {
        all_means[i, 2] <- NA
      }
      else {
        all_means[i, 2] <- inbred_means[inbred_means$Genotype == ID, "Estimate"]
      }
    }
    else if (paste0(P1, "x", P2) %in% scas$Genotype) {
      all_means[i, 2] <- scas[scas$Genotype == paste0(P1,"x",P2), "Estimate"]
    }
    else if (paste0(P2, "x", P1) %in% scas$Genotype) {
      all_means[i, 2] <- scas[scas$Genotype == paste0(P2,"x",P1), "Estimate"]
    }
    else {
      p1_value <- gcas[gcas$Genotype == P1, "GCA"]
      p2_value <- gcas[gcas$Genotype == P2, "GCA"]
      all_means[i, 2] <- hybrid_mean + p1_value + p2_value
    }
    
  }

  colnames(all_means) <- c("Genotype", "Estimate")
  return (all_means)
  
}

all_means_classic <- ProtectedList("all_means_classic")
# debugonce(PredictMeansClassic) ### DEBUG
all_means_classic[[trait_name]] <- PredictMeansClassic(ped_table, 
                                                          scas[[trait_name]], 
                                                          gcas[[trait_name]], 
                                                          inbred_means[[trait_name]])

#     GBLUP
#     Cross validiate GBLUP
#     Classical vs. GBLUP comparison

# 9. Predict synthetics -----
#     Classical (two step)

ConvertMeansToWide <- function (all_means) {
  
#  all_means <- all_means_classic[[trait_name]] ### DEBUG
  for (i in 1:NUM_GENOTYPES) {
    
    all_means[i,"Genotype"] <- paste0(all_means[i,"Genotype"], "x", all_means[i,"Genotype"])
  
  }
  sep_means <- separate(data = all_means, col = Genotype, into = c("Female", "Male"), sep = "x") ### REQUIRES tidyr
  sep_means2 <- separate(data = all_means, col = Genotype, into = c("Male", "Female"), sep = "x")
  all_means <- do.call(rbind,list(sep_means,sep_means2))
  wide_all <- reshape(data = all_means, idvar = "Female", timevar = "Male", direction = "wide", new.row.names = unique(all_means$Female))[,-1]
  names(wide_all) <- gsub("Estimate.", "", names(wide_all))
  
  return (wide_all)
}

wide_means_classic <- ProtectedList("wide_means_classic")
wide_means_classic[[trait_name]] <- ConvertMeansToWide(all_means_classic[[trait_name]]) 

# End Loop
}

Combinadic <- function (n, r, i) {
  
  # Allows slices of combinations in order to parallelize
  # http://msdn.microsoft.com/en-us/library/aa289166(VS.71).aspx
  # http://en.wikipedia.org/wiki/Combinadic
  
  if (i < 1 | i > choose(n, r)) stop ("'i' must be 0 < i <= n!/(n-r)!")
  
  largestV <- function (n, r, i) {
    #v <- n - 1
    v <- n                                  # Adjusted for one-based indexing
    #while (choose(v, r) > i) v <- v - 1
    while (choose(v, r) >= i) v <- v - 1        # Adjusted for one-based indexing
    return (v)
  }
  
  res <- rep(NA, r)
  for (j in 1:r) {
    res[j] <- largestV(n, r, i)
    i <- i - choose(res[j], r)
    n <- res[j]
    r <- r - 1
  }
  res <- res + 1
  return(res)
}

SubsetWright <- function (perse_hybrid, subset) {
  return (mean(perse_hybrid[subset, subset]))
}

TestSyn <- function (perse_hybrid, min, max, rows = "all", type = "p", cl) {
  
  if (min < 2) stop ("min must be at least 2")
  if (nrow(perse_hybrid) < 2) stop ("At least two inbreds must be included in perse_hybrid matrix")
  if (nrow(perse_hybrid) != ncol(perse_hybrid)) stop ("perse_hybrid matrix must be square")
  if (max > nrow(perse_hybrid)) stop ("max cannot be greater than the size of the matrix")
  
  cl <- makeCluster(detectCores() - 1) ### REQUIRES parallel
  clusterExport(cl, list("SubsetWright", "Combinadic"))
  
  total <- sum(apply(X = array(min:max), MARGIN = 1, FUN = choose, n = nrow(perse_hybrid)))
  n <- nrow(perse_hybrid)
  start <- 1
  stop <- choose(n, min)
  
  syn_data <- big.matrix(nrow = total, ncol = max + 1) ### REQUIRES big.matrix
  
  print.noquote (paste("Testing ", total, " combinations", sep = ""))
  
  for (i in min:max)
  {
    
    #add inbred numbers to syn_data
    inbreds <- t(parSapplyLB(cl = cl, X = 1 : ((stop - start) + 1), FUN = Combinadic, n = n, r = i))
    syn_data[start:stop, 1:i] <- inbreds
    
    #add sythetic values to syn_data
    syn_data[start:stop, max+1]  <- parApply(cl = cl, X = inbreds, MARGIN = 1,
                                             FUN = SubsetWright, perse_hybrid = perse_hybrid)
    
    start <- stop + 1
    stop <- start + choose(n, i + 1) - 1
    
  }
  
  stopCluster(cl)
  
  #remove NAs
  syn_data <- syn_data[!is.na(syn_data[, ncol(syn_data)]), ]
  
  mpermute(x = syn_data, cols = max + 1, decreasing = TRUE)
  
  if (rows == "all") rows <- nrow(syn_data)
  
  return (syn_data[c(1:rows, (nrow(syn_data) - rows):nrow(syn_data)), ])
}

for (trait_num in 1:length(TRAITS)) {
  
  # trait_num <- 1
  
  trait_name <- SetTraitName(trait_num) ### TO BE LOOPED?

syns <- ProtectedList("syns")
syns[[trait_name]] <- TestSyn(perse_hybrid = data.matrix(wide_means_classic[[trait_name]]),
                              min = 2, max = 7, rows = "all", cl = cl)
}

debugonce(TestSyn)
test <- TestSyn(perse_hybrid = data.matrix(wide_means_classic$PlantHt),
        min = 2, max = 5, rows = "all", cl = cl)

#     GBLUP (two step)
#     GLBUP (one step)

# PART C - for all traits -----

# 10. Merge trait data -----

# 11 . Filter and sort hybrids and synthetics based on index -----

# 12. Create tables and charts -----
#      QC info - QQplots, boxplots, etc
#      ANOVA tables
#      Spearkman correlations
#      Inbred means
#      Tested hybrid means
#      Inbred GCAs
#      Tested hybrid SCAs
#      Predicted hybrid means (Classical and GBLUP)
#      GBLUP cross-validation 
#      Predicted synthetic means (filtered)

printSyns <- function (syns, hybrids) {
  
  cat(paste0("Top ", (nrow(syns) - 1) / 2," and bottom ",(nrow(syns) - 1) / 2," predicted synthetics\n\n"))
  for (i in 1:nrow(syns)) {
    
    syn_row <- syns[i, 1:ncol(syns) - 1]
    syn_row <- syn_row[!is.na(syn_row)]
    print(hybrids[syn_row, syn_row])
    cat(paste0("Mean: ", round(mean(data.matrix(hybrids[syn_row, syn_row])), 4), "\n"))
    cat(paste0("Std Dev: ", round(sd(data.matrix(hybrids[syn_row, syn_row])), 4), "\n\n"))
    
  }
  
}
printSyns(syns[[trait_name]], wide_means_classic[[trait_name]])


# Saving and retriving environment -----
# save.image(file = "sweetcornobjects.rda")
# load("sweetcornobjects.rda")

# SCRATCH ------
# pred_names <- c("Location","Geno","LocRep","LocBlock")
# library(glmulti)
# 
# model_select <- glmulti(y=trait_name,xr=pred_names,data=pheno_data)
# model_select <- glmulti(y=trait_name,xr=pred_names,data=pheno_data,level=2,marginality=TRUE, method="d")




