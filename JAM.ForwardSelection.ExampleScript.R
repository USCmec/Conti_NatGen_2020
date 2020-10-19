################################################################################
# This is example code for the JAM analysis found in Conti et al. Nat. Genetics
# 
# Code is related to the following specific analysis: A marginal P-value from the multi-ethnic 
# meta-analysis GWAS less than 5.0x10-8 in either the population-specific or multiethnic analysis
# was used to define statistically significant genetic associations, with regions bounded within
# +/- 800 kb from the most significant variant. To determine if multiple independent associations
# exist within each region, we implemented a forward stepwise selection starting with the inclusion
# of the lowest multiethnic marginal P-value into a multivariate logistic regression model. 
# We used Joint Analysis of Marginal summary statistics (JAM) to obtain population-specific 
# conditional summary statistics from multivariate models. Conditional statistics were combined 
# with an inverse-variance weighted fixed effects meta-analysis to obtain multiethnic conditional 
# summary statistics (Supplementary Table 4). Variants with a conditional multiethnic 
# P-value < 5.0x10-8 were retained in the model. We excluded variants with a marginal multiethnic
# P-value > 5.0x10-4, MAF < 1% in all four populations, and correlation r2 ≥ 0.2 to any variants 
# included in the current model at each step.
#
#
# This code is intended to serve as an example and guide of how JAM was used in the forward selection analysis within a given region and ethnic group. 
# Input data is described but cannot be broadly shared due to IRB restrictions on individual level data. 
# Example code highlights the main steps for the JAM analysis and assumes the same SNPs are available within each ethnic group.
# Author: David Conti
# Date: 10/19/2020


# R libraries needed for this specific analysis
library(Rmpfr)
library(stringr)
library(mvtnorm)

# Data and Input:
# E.W: Reference data for genotypes (N X M) for ethnicity E: N number of individuals and M number of markers
# E.M: number of markers within the region
# E.beta: vector of length M with marginal effect estimates from the ethnic-specific meta-analysis in ethnicity E (note- the index.SNP is the first SNP in vector)
# E.se: vector of length M with standard errors corresponding to the marginal effect estimates from the ethnic-specific meta-analysis in ethnicity E
# E.pcases: proportion of cases from the ethnic-spcific meta analysis in ethnicity E
# E.N: number of individuals (cases and controls) from the ethnic-spcific meta analysis in ethnicity E

############### Run Binary transformation on each population
# function for binary transformation
Binary_Trans = function(Marg_LogORs,SE_LogORs,W,P_Cases,N){
  # 1) Infer z-scores from log-OR, SE and N
  Z_Scores = Marg_LogORs/(SE_LogORs*sqrt(N))
  # 2) Get SNP SDs from reference data and divide Z-scores to get allelic effects
  # Calculate SNP SDs in reference data
  SD_in_Ref = apply(W,2,function(v) sd(v) ) 
  # Divide standardised linear effects by SNP standard deviations
  Transformed_Linear_Betahats = Z_Scores/SD_in_Ref 
  # 3) Multiply by trait SD for effect on trait scale
  # Transformed_Linear_Betahats: NEW standardised "linear" effects to give to JAM (transformed to the linear scale)
  # Multiply by trait SD for effect on trait scale
  Transformed_Linear_Betahats = Transformed_Linear_Betahats*sqrt(P_Cases*(1-P_Cases)) 
  return(Transformed_Linear_Betahats)
}

# run for each ethnicity E={European, African, East Asian, Hispanic}
E.betahats <- Binary_Trans(E.beta, E.se, E.W, E.pcases, E.N)

############### create JAM data
# function for creating JAM data
getJAM.Variables <- function(betas, W, N, ridgeTerm=F) {
  M <- length(betas)  # transformation of betas from logistic model - Assume transformation before input
  p <- apply(W, 2, mean)/2
  n0 <- N*(1-p)^2
  n1 <- N*2*p*(1-p)
  n2 <- N*p^2
  y0 <- -(n1*betas+2*n2*betas)/(n0+n1+n2)
  y1 <- y0+betas
  y2 <- y0+2*betas
  z <- n1*y1 + 2*n2*y2
  W.scaled <- scale(W, center=T, scale=F)
  WW <- t(W.scaled)%*%W.scaled
  Dm <- 2*p*(1-p)*N
  D_sqrt = diag(sqrt(Dm))
  Dw_sqrt_inv = diag(1/sqrt(diag(WW)))
  scaled.WW = D_sqrt %*% Dw_sqrt_inv  %*% t(W.scaled) %*% W.scaled %*% Dw_sqrt_inv %*% D_sqrt
  ridgeValue <- ifelse(ridgeTerm, min(1,min(diag(scaled.WW)*.001)), 0)
  scaled.WW <- scaled.WW + ridgeValue*diag(M)
  L <- chol(scaled.WW) 
  zL <- solve(t(L))%*%z
  return(list(zL=zL, L=L, z=z, xTx=scaled.WW))
}

# run for each ethnicity E={European, African, East Asian, Hispanic}
E.jam <- getJAM.Variables(betas=E.betahats, W=E.W, N=E.N, ridgeTerm=T)

############### perform conditional selection for each ethnicity E={European, African, East Asian, Hispanic}
#Runs conditional JAM models for specified model size
conditional.Fun <- function(E.jam,group="E",M){
  z <- jam$z
  xTx <- jam$xTx
  n.prior <- 1
  p <- get(paste0("pcases.",group))
  b <- ((n.prior - 1) / 2 + 1) * p * (1 - p)
  a <- b / (p * (1 - p)) + 1
  traitVariance <- p*(1-p)
  #Start JAM
  indexSNP <- 1
  jam.res <- {}
  jam.betas <- {}
  for(m in 1:M) {
    if(m==indexSNP) { 
      jam.res <- rbind(jam.res,c(m, NA, NA))
      jam.betas <- rbind(jam.betas,c(NA,NA,NA,NA))
    }
    if(m!=indexSNP) {
      SNPs <- c(indexSNP, m)
      z.fit <- z[SNPs]
      xTx.fit <- as.matrix(xTx[SNPs, SNPs])
      r.jam <- fitJAM_v2(z=z.fit, xTx=xTx.fit, a_sigma=a, b_sigma=b, g=N, trait.variance=traitVariance, n.people=N)
      Z.score <- abs(r.jam$betas/r.jam$se)
      p.value <- 2*(pnorm(as(Z.score,"mpfr"),lower.tail = F))
      jam.res <- rbind(jam.res, c(m, as.numeric(p.value[1]),as.numeric(p.value[2])))
      jam.betas <- rbind(jam.betas,cbind(as.numeric(r.jam$betas[1]),as.numeric(r.jam$betas[2]),
                                         as.numeric(r.jam$se[1]),as.numeric(r.jam$se[2])))
    }
  }
  assign(paste0("jam.res.",group),jam.res)
  jam.betas <- as.data.frame(jam.betas)
  names(jam.betas) <- c("index.beta","snp.beta","index.se","snp.se")
  assign(paste0("jam.betas.",group),jam.betas)
  jam.results <- as.data.frame(cbind(jam.betas,jam.res))
  return(jam.results)
}

# run for each ethnicity E={European, African, East Asian, Hispanic}
E.cond <- as.data.frame(cbind(names(betahats.E),conditional.Fun(E.jam,"E",E.M)))


############### Conditional meta-analysis of the results
# Take each conditional analysis from every ethnic group
# Meta-analyze the results with the SNP list
# Return the original snp list with meta p-values and ethnic-specific p-values
meta.jam <- function(snp.list,pop1,pop2,pop3,pop4){
  snp.list <- as.data.frame(snp.list)
  names(snp.list) <- "query"
  snp.list$query <- as.character(snp.list$query)
  snp.list$index <- 1:nrow(snp.list)
  #Merge data sets
  Fdf0 <- merge(snp.list,pop1,by="query",all.x=T)
  colnames(Fdf0)[3:8] <- c("eur.beta.i","eur.beta.s","eur.se.i","eur.se.s","eur.p.index","eur.p.snp")
  Fdf1 <- merge(Fdf0,pop2,by="query",all.x=T)
  colnames(Fdf1)[9:14] <- c("afr.beta.i","afr.beta.s","afr.se.i","afr.se.s","afr.p.index","afr.p.snp")
  Fdf2 <- merge(Fdf1,pop3,by="query",all.x=T)
  colnames(Fdf2)[15:20] <- c("amr.beta.i","amr.beta.s","amr.se.i","amr.se.s","amr.p.index","amr.p.snp")
  Fdf3 <- merge(Fdf2,pop4,by="query",all.x=T)
  colnames(Fdf3)[21:26] <- c("eas.beta.i","eas.beta.s","eas.se.i","eas.se.s","eas.p.index","eas.p.snp")
  Fdf4 <- Fdf3[order(Fdf3$index),]
  
  #Create matrix of values to sum up in rows with na.rm=T
  #Index - Betas
  Fdf4$eur.sum.i <- (Fdf4$eur.beta.i/(Fdf4$eur.se.i)^2)
  Fdf4$afr.sum.i <- (Fdf4$afr.beta.i/(Fdf4$afr.se.i)^2)
  Fdf4$amr.sum.i <- (Fdf4$amr.beta.i/(Fdf4$amr.se.i)^2)
  Fdf4$eas.sum.i <- (Fdf4$eas.beta.i/(Fdf4$eas.se.i)^2)
  #SNP - Betas
  Fdf4$eur.sum.s <- (Fdf4$eur.beta.s/(Fdf4$eur.se.s)^2)
  Fdf4$afr.sum.s <- (Fdf4$afr.beta.s/(Fdf4$afr.se.s)^2)
  Fdf4$amr.sum.s <- (Fdf4$amr.beta.s/(Fdf4$amr.se.s)^2)
  Fdf4$eas.sum.s <- (Fdf4$eas.beta.s/(Fdf4$eas.se.s)^2)
  #Index - Variance
  Fdf4$eur.var.i <- (1/(Fdf4$eur.se.i)^2)
  Fdf4$afr.var.i <- (1/(Fdf4$afr.se.i)^2)
  Fdf4$amr.var.i <- (1/(Fdf4$amr.se.i)^2)
  Fdf4$eas.var.i <- (1/(Fdf4$eas.se.i)^2)
  #SNP - Variance
  Fdf4$eur.var.s <- (1/(Fdf4$eur.se.s)^2)
  Fdf4$afr.var.s <- (1/(Fdf4$afr.se.s)^2)
  Fdf4$amr.var.s <- (1/(Fdf4$amr.se.s)^2)
  Fdf4$eas.var.s <- (1/(Fdf4$eas.se.s)^2)
  
  #Combine the betas and se's from the ethnicities that have them
  Fdf4$sum.i <- apply(Fdf4[,c("eur.sum.i","afr.sum.i","amr.sum.i","eas.sum.i")],1,function(R){
    sum(R,na.rm = T)          
  })
  Fdf4$sum.s <- apply(Fdf4[,c("eur.sum.s","afr.sum.s","amr.sum.s","eas.sum.s")],1,function(R){
    sum(R,na.rm = T)
  })
  Fdf4$var.i <- apply(Fdf4[,c("eur.var.i","afr.var.i","amr.var.i","eas.var.i")],1,function(R){
    sqrt(sum(R,na.rm = T))
  })
  Fdf4$var.s <- apply(Fdf4[,c("eur.var.s","afr.var.s","amr.var.s","eas.var.s")],1,function(R){
    sqrt(sum(R,na.rm = T))
  })
  
  #Get meta values
  Fdf4$metaZ.i <- abs(Fdf4$sum.i/Fdf4$var.i)
  Fdf4$metaZ.s <- abs(Fdf4$sum.s/Fdf4$var.s)
  Fdf4$metaP.i <- 2*pnorm(Fdf4$metaZ.i,lower.tail=FALSE)
  Fdf4$metaP.s <- 2*pnorm(Fdf4$metaZ.s,lower.tail=FALSE)
  return(Fdf4)
}

# note snp.list is a vector of SNPs common to all populations
meta.cond <- meta.jam(snp.list,E1.cond,E2.cond,E3.cond,E4.cond)



