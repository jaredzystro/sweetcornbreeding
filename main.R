
### Workflow for data analysis and prediction of 
### Jared Zystro's se sweet corn breeding project

# Style guide -----
# ALLCAPS - global variables
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
  
  WORKING_DIRECTORY <<- "~/Dropbox/PhD/SweetCornBreeding"
  PHENO_DATA_FILE <<- "~/Dropbox/PhD/2015 Trials/Data.csv"
  # PHENO_DATA_FILE <<- "./data/Data.csv"
  MARKER_DATA_FILE <<- "./data/data_parimp_37genos_dropminorhets_dropNAs_thinto250KB.hmp.txt"
  FUNCTIONS_SUBDIRECTORY <<- "./functions"
  NUM_GENOTYPES <<- 38
  TRAITS <<- c("PlantHt","EarHt","Flavor","Tenderness","Husk","Tip","Width","Length","Pollen")
  
  setwd(WORKING_DIRECTORY)
  
  library(lme4)
  library(lmerTest)
  library(reshape)
  library(Hmisc)
  library(estimability)
  
  
}
InitializeProject()

# 2. Import / create data tables -----
#     Load in 2015 and 2016 phenotype data
LoadPhenoData <- function() {
  
  return (read.csv(PHENO_DATA_FILE))
  
}
pheno_data <- LoadPhenoData()

#     Load in inbred marker data
LoadMarkerData <- function() {
  
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
CreateHybridMarkers <- function(marker_data,ped_table) {
  
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
marker_data <- CreateHybridMarkers(marker_data,ped_table)

# PART B - for each trait -----

SetTraitName <- function (trait_number) {

  return (TRAITS[trait_number])
  
}
trait_name <- SetTraitName(1) ### TO BE LOOPED?

#SCRATCH ------
# pred_names <- c("Location","Geno","LocRep","LocBlock")
# library(glmulti)
# 
# model_select <- glmulti(y=trait_name,xr=pred_names,data=pheno_data)
# model_select <- glmulti(y=trait_name,xr=pred_names,data=pheno_data,level=2,marginality=TRUE, method="d")


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
gxe_model <- list()
gxe_model[[trait_name]] <- CreateGxEModel(trait_name, pheno_data)
gxe_anova <- list()
gxe_anova[[trait_name]] <- anova(gxe_model[[trait_name]]) ### REQUIRES lmerTest

#     Spearman rank correlation
GxESpearman <- function (trait_name, pheno_data) {
  
  trait_loc <- aggregate(as.formula(paste0(trait_name," ~ Geno:Env")), data = pheno_data, mean, na.rm = TRUE)
  trait_loc <- cast(trait_loc, Geno ~ Env, value = c(as.character(trait_name))) ### REQUIRES reshape
  
  return (rcorr(as.matrix(trait_loc), type = "spearman")) ### REQUIRES Hmisc
  
}
spear_corr <- list ()
spear_corr[[trait_name]] <- GxESpearman(trait_name, pheno_data)

# 6. ANOVAs ----
#     Check ANOVAs
CreateCheckModel<- function (trait_name, check_pheno) {
  
  formula_name <- as.formula(paste0(trait_name, " ~ EnvBlock + (1|Geno)"))
  return (lmer(formula_name, data = check_pheno))
  
}
check_model <- list()
check_model[[trait_name]] <- CreateCheckModel(trait_name, check_pheno)
check_anova <- list()
check_anova[[trait_name]] <- anova(check_model[[trait_name]])

#     Adjust means based on checks? TO DO

#     Inbred ANOVAs
CreateInbredModel<- function (trait_name, inbred_pheno) {
  
  formula_name <- as.formula(paste0(trait_name, " ~ Geno * Env + (1 | EnvRep)"))
  return (lmer(formula_name, data = inbred_pheno))
  
}
inbred_model <- list()
inbred_model[[trait_name]] <- CreateInbredModel(trait_name, inbred_pheno)
inbred_anova <- list()
inbred_anova[[trait_name]] <- anova(inbred_model[[trait_name]])

#     Hybrid ANOVAs
CreateHybridModel<- function (trait_name, hybrid_pheno) {
  
  # Y = μ + Env + Rep(Env) + Block(Rep x Env) + Set + GCA1(Set) + GCA2(Set) + 
  # SCA(Set) + Env x Set + Env x [GCA1(Set)] + Env x [GCA2(Set)] + Env x [SCA(Set)] + ε
  
  
  formula_name <- as.formula(paste0(trait_name, " ~ Set + MaleGeno * FemaleGeno + (1 | EnvRep)"))
  return (lmer(formula_name, data = hybrid_pheno))
  
}
hybrid_model <- list()
hybrid_model[[trait_name]] <- CreateHybridModel(trait_name, hybrid_pheno)
hybrid_anova <- list()
hybrid_anova[[trait_name]] <- anova(hybrid_model[[trait_name]])

# 7. Collect phenotype values -----
#     Inbred means

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
  dii_model[[1]] <- lmer (formula_name, data = dii_sets[[1]]) 
  dii_model[[2]] <- lmer (formula_name, data = dii_sets[[2]]) 
  dii_model[[3]] <- lmer (formula_name, data = dii_sets[[3]]) 
  dii_model[[4]] <- lmer (formula_name, data = dii_sets[[4]]) 
  
  return (dii_model) #DEBUG
  
} ### DEBUG

# debugonce(CreateSetModels) ### DEBUG
dii_model <- CreateSetModels(trait_name, dii_sets) ### DEBUG

GetGCAs <- function (trait_name, dii_model) {
  
  gca <- list()
  
  for (i in 1:4) {
    
 #   i <- 1 ### DEBUG
    
    sca_lsmeans <- (lsmeans(dii_model[[i]], "MaleGeno:FemaleGeno"))$lsmeans.table ### REQUIRES estimability
    grandmean <- mean(sca_lsmeans$Estimate)
    
    fgca_lsmeans <- (lsmeans(dii_model[[i]], "FemaleGeno"))$lsmeans.table
    fgca <- data.frame (Genotype=fgca_lsmeans$FemaleGeno, GCA = fgca_lsmeans$Estimate - grandmean, StdErr = fgca_lsmeans$`Standard Error`)
    
    mgca_lsmeans <- (lsmeans(dii_model[[i]],"MaleGeno"))$lsmeans.table
    mgca <- data.frame (Genotype = mgca_lsmeans$MaleGeno, GCA = mgca_lsmeans$Estimate - grandmean, StdErr = mgca_lsmeans$`Standard Error`)
    
    gca[[i]] <- merge(fgca, mgca, all=TRUE)
  }
  
  return(do.call(rbind,gca))
  
}

gcas <- GetGCAs(trait_name, dii_model)

#     Tested hybrid means
#     Inbred GCAs
#     Tested hybrid SCAs

GetSCAs <- function (trait_name, dii_model) {
  
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
    
    sca[[i]] <- data.frame (Genotype=sca_lsmeans$Genotype, Estimate=sca_lsmeans$Estimate, SCA=sca_lsmeans$SCA, StdErr=sca_lsmeans$`Standard Error`)
  }
  
  scas<-do.call(rbind,sca)
  
  return(scas)
}

# debugonce(GetSCAs) ### DEBUG
scas <- GetSCAs(trait_name, dii_model)

# 8. Predict untested hybrid means -----
#     Classical (u + GCAf + GCAm)



#     GBLUP
#     Cross validiate GBLUP
#     Classical vs. GBLUP comparison

# 9. Predict synthetics -----
#     Classical (two step)
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






