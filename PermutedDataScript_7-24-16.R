############################################################
########### Options, Packages, and Files to Load ###########  
############################################################

# Options #
options(stringsAsFactors = FALSE)

# Libraries #
require(FD)
require(picante)
require(cooccur)

# Files #
traits<-read.csv("TraitMatrix_Imputed_2016-06-27_pruned.csv", check.names=FALSE, row.names=1)
envs<-read.csv("isomatrix_pruned.csv", check.names=FALSE, row.names=1)
envt<-as.matrix(t(envs))
load("PermutedTraitxTrait_IsoList")

envs2 = envs
envs2$Species = rownames(envs2)

Isolations <- rownames(envt)

############## This is to generate permuted data - Permuted Data gets loaded above ############## 
# PermutedTraitxTrait_list= list()
# for(iter in 1:1000){
#   
#   PermutedTraitxTrait_IsoList[[iter]] = list()
#   tempMatrix <- randomizeMatrix(traits, null.model="independentswap", iterations=1000)
#   tempMatrix = data.frame(tempMatrix, check.names = F)
#   
#   #CWM<-functcomp(tempMatrix, envt, CWM.type="all")
#   #PermutedTraitxTrait_IsoList[[iter]]$CWM <- CWM 
#   
#   tempMatrix$Species = rownames(tempMatrix)
#   CompiledData = merge(envs2, tempMatrix, by = "Species")
#   PermutedTraitxTrait_IsoList[[iter]]$Permutation <- iter 
#   CompiledData= CompiledData[,-which(colnames(CompiledData) == "Species")]
#   PermutedTraitxTrait_IsoList[[iter]]$Permutation = iter
#   PermutedTraitxTrait_IsoList[[iter]]$CompiledData <- CompiledData 
#   print(iter)
# }
# 
# save(PermutedTraitxTrait_IsoList, file = paste("PermutedTraitxTrait_IsoList_", Sys.Date(), ".rda", sep = ""))

#############################################################################################################
########### Calculates counts of each type of trait presence/absence [(+/+), (-/-), (+/-), (-/+)] ########### 
########### for all trait pairs where PP is (+/+), NN is (-/-), PN is (+/-) and NP is (-/+)       ###########  
#############################################################################################################

# Calculates trait/trait counts for permuted data #
PermutedAllAssoc_list = list()
for(i in 1:1000){
  PermutedAllAssoc_list[[i]] = list()
  tempDF = PermutedTraitxTrait_IsoList[[i]]$CompiledData
  tempMatrix_PP = matrix(0, ncol = ncol(tempDF), nrow = ncol(tempDF))
  colnames(tempMatrix_PP) = colnames(tempDF)
  rownames(tempMatrix_PP) = colnames(tempDF)
  
  tempMatrix_NN = matrix(0, ncol = ncol(tempDF), nrow = ncol(tempDF))
  colnames(tempMatrix_NN) = colnames(tempDF)
  rownames(tempMatrix_NN) = colnames(tempDF)
  
  tempMatrix_PN = matrix(0, ncol = ncol(tempDF), nrow = ncol(tempDF))
  colnames(tempMatrix_PN) = colnames(tempDF)
  rownames(tempMatrix_PN) = colnames(tempDF)
  
  tempMatrix_NP = matrix(0, ncol = ncol(tempDF), nrow = ncol(tempDF))
  colnames(tempMatrix_NP) = colnames(tempDF)
  rownames(tempMatrix_NP) = colnames(tempDF)
  
  for(j in 1:(ncol(tempDF)-1)){
    for(l in (j+1):ncol(tempDF)){
      POS = length(which(tempDF[,j] == 1 & tempDF[,l] == 1))
      NEG = length(which(tempDF[,j] == 0 & tempDF[,l] == 0))
      DIFF_PN = length(which(tempDF[,j] == 1 & tempDF[,l] == 0))
      DIFF_NP = length(which(tempDF[,j] == 0 & tempDF[,l] == 1))
      tempMatrix_PP[j,l] = POS
      tempMatrix_NN[j,1] = NEG
      tempMatrix_PN[j,l] = DIFF_PN
      tempMatrix_NP[j,l] = DIFF_NP
    }
  }
  PermutedAllAssoc_list[[i]]$Both_P = tempMatrix_PP
  PermutedAllAssoc_list[[i]]$Both_N = tempMatrix_NN
  PermutedAllAssoc_list[[i]]$Diff_PN = tempMatrix_PN
  PermutedAllAssoc_list[[i]]$Diff_NP = tempMatrix_NP
  print(i)
}

# Calculates trait/trait counts for observed data #
TraitCounts_df = data.frame(TraitA = character(), 
                            TraitB = character(), 
                            Both_P = numeric(),
                            Both_N = numeric(),
                            Diff_PN = numeric(), 
                            Diff_NP = numeric())
k = 1
for(i in 1:(ncol(traits)-1)){
  for(j in (i+1):ncol(traits)){
    TraitCounts_df[k,1] = colnames(traits)[i]
    TraitCounts_df[k,2] = colnames(traits)[j]
    TraitCounts_df[k,3] = length(which(traits[,i] == 1 & traits[,j] == 1))
    TraitCounts_df[k,4] = length(which(traits[,i] == 0 & traits[,j] == 0))
    TraitCounts_df[k,5] = length(which(traits[,i] == 1 & traits[,j] == 0))
    TraitCounts_df[k,6] = length(which(traits[,i] == 0 & traits[,j] == 1))
    k = k+1
  }
}

write.csv(TraitCounts_df, file = paste("TraitCounts_df_", Sys.Date(), ".csv", sep = ""), row.names = F)

###### Determine positive associations between Trait pairs ###### 
# Calculate the pooled Positive associations between  Trait pairs [(+/+) and (-/-)] #
PooledAssoc_Pos_df = data.frame(Trait_A = character(),
                                Trait_B = character(),
                                Observed_PP = numeric(),
                                Percent_PP = numeric())
k = 1
for(i in 1:(length(Traits)-1)){
  for(j in (i +1):length(Traits)){
    tempObs1 = TraitCounts_df[which(TraitCounts_df[,1] == Traits[[i]] & TraitCounts_df[,2] == Traits[[j]]),3]
    tempObs2 = TraitCounts_df[which(TraitCounts_df[,1] == Traits[[i]] & TraitCounts_df[,2] == Traits[[j]]),4]
    tempList = list()
    for(perm in 1:1000){
      tempDF1 = PermutedAllAssoc_list[[perm]]$Both_P
      tempDF2 = PermutedAllAssoc_list[[perm]]$Both_N
      COLI= which(colnames(tempDF) == Traits[[i]])
      ROWJ= which(rownames(tempDF) == Traits[[j]])
      COLJ= which(colnames(tempDF) == Traits[[j]])
      ROWI= which(rownames(tempDF) == Traits[[i]])
      OBS1 = tempDF1[ROWJ,COLI]
      OBS2 = tempDF1[ROWI,COLJ]
      if(OBS1 > OBS2){
        permOBS = OBS1
      }else{
        permOBS = OBS2
      }
      OBS3 = tempDF2[ROWJ,COLI]
      OBS4 = tempDF2[ROWI,COLJ]
      if(OBS3 > OBS4){
        permOBS2 = OBS3
      }else{
        permOBS2 = OBS4
      }
      tempList[perm] = permOBS + permOBS2
    }
    temporary = unlist(tempList)
    PooledAssoc_Pos_df[k,1] = Traits[[i]]
    PooledAssoc_Pos_df[k,2] = Traits[[j]]
    tempObs = tempObs1 + tempObs2
    PooledAssoc_Pos_df[k,3] = tempObs
    PooledAssoc_Pos_df[k,4] = ecdf(temporary)(tempObs)
    k = k + 1
  }
}

write.csv(PooledAssoc_Pos_df, file = paste("PooledAssoc_Pos_df_", Sys.Date(), ".csv", sep = ""), row.names = F)

# Determine the number of observed (+/+)'s for Trait pairs for the observed data and compare it to the permuted data # 
Assoc_PP_Pos_df = data.frame(Trait_A = character(),
                             Trait_B = character(),
                             Observed_PP = numeric(),
                             Percent_PP = numeric())
k = 1
for(i in 1:(length(Traits)-1)){
  for(j in (i +1):length(Traits)){
    tempObs = TraitCounts_df[which(TraitCounts_df[,1] == Traits[[i]] & TraitCounts_df[,2] == Traits[[j]]),3]
    tempList = list()
    for(perm in 1:1000){
      tempDF1 = PermutedAllAssoc_list[[perm]]$Both_P
      COLI= which(colnames(tempDF) == Traits[[i]])
      ROWJ= which(rownames(tempDF) == Traits[[j]])
      COLJ= which(colnames(tempDF) == Traits[[j]])
      ROWI= which(rownames(tempDF) == Traits[[i]])
      OBS1 = tempDF1[ROWJ,COLI]
      OBS2 = tempDF1[ROWI,COLJ]
      if(OBS1 > OBS2){
        permOBS = OBS1
      }else{
        permOBS = OBS2
      }
      tempList[perm] = permOBS
    }
    temporary = unlist(tempList)
    Assoc_PP_Pos_df[k,1] = Traits[[i]]
    Assoc_PP_Pos_df[k,2] = Traits[[j]]
    Assoc_PP_Pos_df[k,3] = tempObs
    Assoc_PP_Pos_df[k,4] = ecdf(temporary)(tempObs)
    k = k + 1
  }
}

write.csv(Assoc_PP_Pos_df, file = paste("Assoc_PP_Pos_df_", Sys.Date(), ".csv", sep = ""), row.names = F)

# Determine the number of observed (-/-)'s for Trait pairs for the observed data and compare it to the permuted data # 
Assoc_NN_Pos_df = data.frame(Trait_A = character(),
                             Trait_B = character(),
                             Observed_NN = numeric(),
                             Percent_NN = numeric())
k = 1
for(i in 1:(length(Traits)-1)){
  for(j in (i +1):length(Traits)){
    tempObs = TraitCounts_df[which(TraitCounts_df[,1] == Traits[[i]] & TraitCounts_df[,2] == Traits[[j]]),4]
    tempList = list()
    for(perm in 1:1000){
      tempDF1 = PermutedAllAssoc_list[[perm]]$Both_N
      COLI= which(colnames(tempDF) == Traits[[i]])
      ROWJ= which(rownames(tempDF) == Traits[[j]])
      COLJ= which(colnames(tempDF) == Traits[[j]])
      ROWI= which(rownames(tempDF) == Traits[[i]])
      OBS1 = tempDF1[ROWJ,COLI]
      OBS2 = tempDF1[ROWI,COLJ]
      if(OBS1 > OBS2){
        permOBS = OBS1
      }else{
        permOBS = OBS2
      }
      tempList[perm] = permOBS
    }
    temporary = unlist(tempList)
    Assoc_NN_Pos_df[k,1] = Traits[[i]]
    Assoc_NN_Pos_df[k,2] = Traits[[j]]
    Assoc_NN_Pos_df[k,3] = tempObs
    Assoc_NN_Pos_df[k,4] = ecdf(temporary)(tempObs)
    k = k + 1
  }
}

write.csv(Assoc_NN_Pos_df, file = paste("Assoc_NN_Pos_df_", Sys.Date(), ".csv", sep = ""), row.names = F)

###### Determine negative associations between Trait pairs ###### 
# Calculate the pooled Negative associations between  Trait pairs [(+/-) and (-/+)] #
PooledAssoc_Neg_df = data.frame(Trait_A = character(),
                                Trait_B = character(),
                                Observed_PP = numeric(),
                                Percent_PP = numeric())
k = 1
for(i in 1:(length(Traits)-1)){
  for(j in (i +1):length(Traits)){
    tempObs1 = TraitCounts_df[which(TraitCounts_df[,1] == Traits[[i]] & TraitCounts_df[,2] == Traits[[j]]),5]
    tempObs2 = TraitCounts_df[which(TraitCounts_df[,1] == Traits[[i]] & TraitCounts_df[,2] == Traits[[j]]),6]
    tempList = list()
    for(perm in 1:1000){
      tempDF1 = PermutedAllAssoc_list[[perm]]$Diff_PN
      tempDF2 = PermutedAllAssoc_list[[perm]]$Diff_NP
      COLI= which(colnames(tempDF) == Traits[[i]])
      ROWJ= which(rownames(tempDF) == Traits[[j]])
      COLJ= which(colnames(tempDF) == Traits[[j]])
      ROWI= which(rownames(tempDF) == Traits[[i]])
      OBS1 = tempDF1[ROWJ,COLI]
      OBS2 = tempDF1[ROWI,COLJ]
      if(OBS1 > OBS2){
        permOBS = OBS1
      }else{
        permOBS = OBS2
      }
      OBS3 = tempDF2[ROWJ,COLI]
      OBS4 = tempDF2[ROWI,COLJ]
      if(OBS3 > OBS4){
        permOBS2 = OBS3
      }else{
        permOBS2 = OBS4
      }
      tempList[perm] = permOBS + permOBS2
    }
    temporary = unlist(tempList)
    PooledAssoc_Neg_df[k,1] = Traits[[i]]
    PooledAssoc_Neg_df[k,2] = Traits[[j]]
    tempObs = tempObs1 + tempObs2
    PooledAssoc_Neg_df[k,3] = tempObs
    PooledAssoc_Neg_df[k,4] = ecdf(temporary)(tempObs)
    k = k + 1
  }
}

write.csv(PooledAssoc_Neg_df, file = paste("PooledAssoc_Neg_df_", Sys.Date(), ".csv", sep = ""), row.names = F)

# Determine the number of observed (+/-)'s for Trait pairs for the observed data and compare it to the permuted data # 
Assoc_PN_Neg_df = data.frame(Trait_A = character(),
                             Trait_B = character(),
                             Observed_PN = numeric(),
                             Percent_PN = numeric())
k = 1
for(i in 1:(length(Traits)-1)){
  for(j in (i +1):length(Traits)){
    tempObs = TraitCounts_df[which(TraitCounts_df[,1] == Traits[[i]] & TraitCounts_df[,2] == Traits[[j]]),5]
    tempList = list()
    for(perm in 1:1000){
      tempDF1 = PermutedAllAssoc_list[[perm]]$Diff_PN
      COLI= which(colnames(tempDF) == Traits[[i]])
      ROWJ= which(rownames(tempDF) == Traits[[j]])
      COLJ= which(colnames(tempDF) == Traits[[j]])
      ROWI= which(rownames(tempDF) == Traits[[i]])
      OBS1 = tempDF1[ROWJ,COLI]
      OBS2 = tempDF1[ROWI,COLJ]
      if(OBS1 > OBS2){
        permOBS = OBS1
      }else{
        permOBS = OBS2
      }
      tempList[perm] = permOBS
    }
    temporary = unlist(tempList)
    Assoc_PN_Neg_df[k,1] = Traits[[i]]
    Assoc_PN_Neg_df[k,2] = Traits[[j]]
    Assoc_PN_Neg_df[k,3] = tempObs
    Assoc_PN_Neg_df[k,4] = ecdf(temporary)(tempObs)
    k = k + 1
  }
}

write.csv(Assoc_PN_Neg_df, file = paste("Assoc_PN_Neg_df_", Sys.Date(), ".csv", sep = ""), row.names = F)

# Determine the number of observed (-/+)'s for Trait pairs for the observed data and compare it to the permuted data # 
Assoc_NP_Neg_df = data.frame(Trait_A = character(),
                             Trait_B = character(),
                             Observed_NP = numeric(),
                             Percent_NP = numeric())
k = 1
for(i in 1:(length(Traits)-1)){
  for(j in (i +1):length(Traits)){
    tempObs = TraitCounts_df[which(TraitCounts_df[,1] == Traits[[i]] & TraitCounts_df[,2] == Traits[[j]]),6]
    tempList = list()
    for(perm in 1:1000){
      tempDF1 = PermutedAllAssoc_list[[perm]]$Diff_NP
      COLI= which(colnames(tempDF) == Traits[[i]])
      ROWJ= which(rownames(tempDF) == Traits[[j]])
      COLJ= which(colnames(tempDF) == Traits[[j]])
      ROWI= which(rownames(tempDF) == Traits[[i]])
      OBS1 = tempDF1[ROWJ,COLI]
      OBS2 = tempDF1[ROWI,COLJ]
      if(OBS1 > OBS2){
        permOBS = OBS1
      }else{
        permOBS = OBS2
      }
      tempList[perm] = permOBS
    }
    temporary = unlist(tempList)
    Assoc_NP_Neg_df[k,1] = Traits[[i]]
    Assoc_NP_Neg_df[k,2] = Traits[[j]]
    Assoc_NP_Neg_df[k,3] = tempObs
    Assoc_NP_Neg_df[k,4] = ecdf(temporary)(tempObs)
    k = k + 1
  }
}

write.csv(Assoc_NP_Neg_df, file = paste("Assoc_NP_Neg_df_", Sys.Date(), ".csv", sep = ""), row.names = F)

envs3 = envs2[,-which(colnames(envs2) == "Species")]
traits2 = traits[,-which(colnames(traits) == "Species")]
CombData = merge(envs2, traits, by = "Species")
CombData = CombData[,-which(colnames(CombData) == "Species")]

############################################################################
########### Analysis of Isolation Environment/Trait Associations ###########  
############################################################################

# Calculate the number of observed [(+/+), (-,-), (+,-), and (-/+)] for each Isolation/Trait pair # 
IsolationXTraitCounts_df = data.frame(Isolation = character(), 
                            Trait = character(), 
                            Both_P = numeric(),
                            Both_N = numeric(),
                            Diff_PN = numeric(), 
                            Diff_NP = numeric())
k = 1
for(i in 1:length(Isolations)){
  for(j in 1:length(Traits)){
    IsolationXTraitCounts_df[k,1] = Isolations[[i]]
    IsolationXTraitCounts_df[k,2] = Traits[[j]]
	tempISO = data.frame(CombData[,which(colnames(CombData) == Isolations[i])])
	tempTRAIT = data.frame(CombData[,which(colnames(CombData) == Traits[j])])
	tempData = cbind(tempISO, tempTRAIT)
    IsolationXTraitCounts_df[k,3] = length(which(tempData[,1] == 1 & tempData[,2] == 1))
    IsolationXTraitCounts_df[k,4] = length(which(tempData[,1] == 0 & tempData[,2] == 0))
    IsolationXTraitCounts_df[k,5] = length(which(tempData[,1] == 1 & tempData[,2] == 0))
    IsolationXTraitCounts_df[k,6] = length(which(tempData[,1] == 0 & tempData[,2] == 1))
    k = k+1
  }
}

write.csv(IsolationXTraitCounts_df, file = paste("IsolationXTraitCounts_df_", Sys.Date(), ".csv", sep = ""), row.names = F)

###### Determine positive associations between Isolation environments and Traits ###### 
# Calculate the pooled Positive associations between Isolation Environment and Trait [(+/+) and (-/-)] #
PooledIsoTraitAssoc_Pos_df = data.frame(Isolation = character(),
                                Trait_B = character(),
                                Observed_PP = numeric(),
                                Percent_PP = numeric())
k = 1
for(i in 1:length(Isolations)){
  for(j in 1:length(Traits)){
    tempObs1 = IsolationXTraitCounts_df[which(IsolationXTraitCounts_df[,1] == Isolations[[i]] & IsolationXTraitCounts_df[,2] == Traits[[j]]),3]
    tempObs2 = IsolationXTraitCounts_df[which(IsolationXTraitCounts_df[,1] == Isolations[[i]] & IsolationXTraitCounts_df[,2] == Traits[[j]]),4]
    tempList = list()
    for(perm in 1:1000){
      tempDF1 = PermutedAllAssoc_list[[perm]]$Both_P
      tempDF2 = PermutedAllAssoc_list[[perm]]$Both_N
      COLI= which(colnames(tempDF) == Isolations[[i]])
      ROWJ= which(rownames(tempDF) == Traits[[j]])
      COLJ= which(colnames(tempDF) == Traits[[j]])
      ROWI= which(rownames(tempDF) == Isolations[[i]])
      OBS1 = tempDF1[ROWJ,COLI]
      OBS2 = tempDF1[ROWI,COLJ]
      if(OBS1 > OBS2){
        permOBS = OBS1
      }else{
        permOBS = OBS2
      }
      OBS3 = tempDF2[ROWJ,COLI]
      OBS4 = tempDF2[ROWI,COLJ]
      if(OBS3 > OBS4){
        permOBS2 = OBS3
      }else{
        permOBS2 = OBS4
      }
      tempList[perm] = permOBS + permOBS2
    }
    temporary = unlist(tempList)
    PooledIsoTraitAssoc_Pos_df[k,1] = Isolations[[i]]
    PooledIsoTraitAssoc_Pos_df[k,2] = Traits[[j]]
    tempObs = tempObs1 + tempObs2
    PooledIsoTraitAssoc_Pos_df[k,3] = tempObs
    PooledIsoTraitAssoc_Pos_df[k,4] = ecdf(temporary)(tempObs)
    k = k + 1
  }
}

write.csv(PooledIsoTraitAssoc_Pos_df, file = paste("PooledIsoTraitAssoc_Pos_df_", Sys.Date(), ".csv", sep = ""), row.names = F)

# Determine the number of observed (+/+)'s for Isolation and Trait for the observed data and compare it to the permuted data # 
AssocIsoTrait_PP_Pos_df = data.frame(Isolation = character(),
                             Trait = character(),
                             Observed_PP = numeric(),
                             Percent_PP = numeric())
k = 1
for(i in 1:length(Isolations)){
  for(j in 1:length(Traits)){
    tempObs = IsolationXTraitCounts_df[which(IsolationXTraitCounts_df[,1] == Isolations[[i]] & IsolationXTraitCounts_df[,2] == Traits[[j]]),3]
    tempList = list()
    for(perm in 1:1000){
      tempDF1 = PermutedAllAssoc_list[[perm]]$Both_P
      COLI= which(colnames(tempDF) == Isolations[[i]])
      ROWJ= which(rownames(tempDF) == Traits[[j]])
      COLJ= which(colnames(tempDF) == Traits[[j]])
      ROWI= which(rownames(tempDF) == Isolations[[i]])
      OBS1 = tempDF1[ROWJ,COLI]
      OBS2 = tempDF1[ROWI,COLJ]
      if(OBS1 > OBS2){
        permOBS = OBS1
      }else{
        permOBS = OBS2
      }
      tempList[perm] = permOBS
    }
    temporary = unlist(tempList)
    AssocIsoTrait_PP_Pos_df[k,1] = Isolations[[i]]
    AssocIsoTrait_PP_Pos_df[k,2] = Traits[[j]]
    AssocIsoTrait_PP_Pos_df[k,3] = tempObs
    AssocIsoTrait_PP_Pos_df[k,4] = ecdf(temporary)(tempObs)
    k = k + 1
  }
}

write.csv(AssocIsoTrait_PP_Pos_df, file = paste("AssocIsoTrait_PP_Pos_df_", Sys.Date(), ".csv", sep = ""), row.names = F)

# Determine the cut-off for signficant postive (+/+) traits in an environment by calculating the values of the data found in the 97.5th% of the distribution #
PosAssocIso_PP_split = splitBy("Isolation", AssocIsoTrait_PP_Pos_df)

PosIsoCutOffs_PP_df = data.frame(Isolation = character(), CutOff_97.5 = numeric(), CutOff_95 = numeric())

for(i in 1:length(PosAssocIso_PP_split)){
	PosIsoCutOffs_PP_df[i,1] = names(PosAssocIso_PP_split)[i]
	PosIsoCutOffs_PP_df[i,2] = quantile(PosAssocIso_PP_split[[i]][,4], 0.975)
	PosIsoCutOffs_PP_df[i,3] = quantile(PosAssocIso_PP_split[[i]][,4], 0.95)
}

write.csv(PosIsoCutOffs_PP_df, file = paste("PosIsoCutOffs_PP_df_", Sys.Date(), ".csv", sep = ""), row.names = F)

# Determine which traits are positive (+,+) significant in an environment #
Pos_SigTraitsIso_PP_list = list()
for(i in 1:nrow(PosIsoCutOffs_PP_df)){
  Pos_SigTraitsIso_PP_list[[i]] = list()
  temp_df = AssocIsoTrait_PP_Pos_df[which(AssocIsoTrait_PP_Pos_df[,1] == PosIsoCutOffs_PP_df[i,1]),]
  Pos_SigTraitsIso_PP_list[[i]]$Isolation = PosIsoCutOffs_PP_df[i,1]
  Pos_SigTraitsIso_PP_list[[i]]$Traits = temp_df[which(temp_df[,4] >= PosIsoCutOffs_PP_df[i,2]),2]
}

# Determine the number of observed (-/-)'s for Isolation and Trait for the observed data and compare it to the permuted data # 
AssocIsoTrait_NN_Pos_df = data.frame(Isolation = character(),
                             Trait = character(),
                             Observed_NN = numeric(),
                             Percent_NN = numeric())
k = 1
for(i in 1:length(Isolations)){
  for(j in 1:length(Traits)){
    tempObs = IsolationXTraitCounts_df[which(IsolationXTraitCounts_df[,1] == Isolations[[i]] & IsolationXTraitCounts_df[,2] == Traits[[j]]),4]
    tempList = list()
    for(perm in 1:1000){
      tempDF1 = PermutedAllAssoc_list[[perm]]$Both_N
      COLI= which(colnames(tempDF) == Isolations[[i]])
      ROWJ= which(rownames(tempDF) == Traits[[j]])
      COLJ= which(colnames(tempDF) == Traits[[j]])
      ROWI= which(rownames(tempDF) == Isolations[[i]])
      OBS1 = tempDF1[ROWJ,COLI]
      OBS2 = tempDF1[ROWI,COLJ]
      if(OBS1 > OBS2){
        permOBS = OBS1
      }else{
        permOBS = OBS2
      }
      tempList[perm] = permOBS
    }
    temporary = unlist(tempList)
    AssocIsoTrait_NN_Pos_df[k,1] = Isolations[[i]]
    AssocIsoTrait_NN_Pos_df[k,2] = Traits[[j]]
    AssocIsoTrait_NN_Pos_df[k,3] = tempObs
    AssocIsoTrait_NN_Pos_df[k,4] = ecdf(temporary)(tempObs)
    k = k + 1
  }
}

write.csv(AssocIsoTrait_NN_Pos_df, file = paste("AssocIsoTrait_NN_Pos_df_", Sys.Date(), ".csv", sep = ""), row.names = F)

# Determine the cut-off for signficant traits in an environment by calculating the values of the data found in the 97.5th% of the distribution #
PosAssocIso_NN_split = splitBy("Isolation", AssocIsoTrait_NN_Pos_df)

PosIsoCutOffs_NN_df = data.frame(Isolation = character(), CutOff_97.5 = numeric(), CutOff_95 = numeric())

for(i in 1:length(PosAssocIso_NN_split)){
	PosIsoCutOffs_NN_df[i,1] = names(PosAssocIso_NN_split)[i]
	PosIsoCutOffs_NN_df[i,2] = quantile(PosAssocIso_NN_split[[i]][,4], 0.975)
	PosIsoCutOffs_NN_df[i,3] = quantile(PosAssocIso_NN_split[[i]][,4], 0.95)
}

write.csv(PosIsoCutOffs_NN_df, file = paste("PosIsoCutOffs_NN_df_", Sys.Date(), ".csv", sep = ""), row.names = F)

# Determine which traits are significant in an environment #
Pos_SigTraitsIso_NN_list = list()
for(i in 1:nrow(PosIsoCutOffs_NN_df)){
  Pos_SigTraitsIso_NN_list[[i]] = list()
  temp_df = AssocIsoTrait_NN_Pos_df[which(AssocIsoTrait_NN_Pos_df[,1] == PosIsoCutOffs_NN_df[i,1]),]
  Pos_SigTraitsIso_NN_list[[i]]$Isolation = PosIsoCutOffs_NN_df[i,1]
  Pos_SigTraitsIso_NN_list[[i]]$Traits = temp_df[which(temp_df[,4] >= PosIsoCutOffs_NN_df[i,2]),2]
}

###### Determine negative associations between Isolation environments and Traits ###### 
# Calculate the pooled Negative associations between Isolation Environment and Trait [(+/-) and (-/+)] #
PooledIsoTraitAssoc_Neg_df = data.frame(Isolation = character(),
                                Trait = character(),
                                Observed_PP = numeric(),
                                Percent_PP = numeric())
k = 1
for(i in 1:length(Isolations)){
  for(j in 1:length(Traits)){
    tempObs1 = IsolationXTraitCounts_df[which(IsolationXTraitCounts_df[,1] == Isolations[[i]] & IsolationXTraitCounts_df[,2] == Traits[[j]]),5]
    tempObs2 = IsolationXTraitCounts_df[which(IsolationXTraitCounts_df[,1] == Isolations[[i]] & IsolationXTraitCounts_df[,2] == Traits[[j]]),6]
    tempList = list()
    for(perm in 1:1000){
      tempDF1 = PermutedAllAssoc_list[[perm]]$Diff_PN
      tempDF2 = PermutedAllAssoc_list[[perm]]$Diff_NP
      COLI= which(colnames(tempDF) == Isolations[[i]])
      ROWJ= which(rownames(tempDF) == Traits[[j]])
      COLJ= which(colnames(tempDF) == Traits[[j]])
      ROWI= which(rownames(tempDF) == Isolations[[i]])
      OBS1 = tempDF1[ROWJ,COLI]
      OBS2 = tempDF1[ROWI,COLJ]
      if(OBS1 > OBS2){
        permOBS = OBS1
      }else{
        permOBS = OBS2
      }
      OBS3 = tempDF2[ROWJ,COLI]
      OBS4 = tempDF2[ROWI,COLJ]
      if(OBS3 > OBS4){
        permOBS2 = OBS3
      }else{
        permOBS2 = OBS4
      }
      tempList[perm] = permOBS + permOBS2
    }
    temporary = unlist(tempList)
    PooledIsoTraitAssoc_Neg_df[k,1] = Isolations[[i]]
    PooledIsoTraitAssoc_Neg_df[k,2] = Traits[[j]]
    tempObs = tempObs1 + tempObs2
    PooledIsoTraitAssoc_Neg_df[k,3] = tempObs
    PooledIsoTraitAssoc_Neg_df[k,4] = ecdf(temporary)(tempObs)
    k = k + 1
  }
}

write.csv(PooledIsoTraitAssoc_Neg_df, file = paste("PooledIsoTraitAssoc_Neg_df_", Sys.Date(), ".csv", sep = ""), row.names = F)

# Determine the number of observed (+/-)'s for Isolation and Trait for the observed data and compare it to the permuted data # 
AssocIsoTrait_PN_Neg_df = data.frame(Isolation = character(),
                             Trait = character(),
                             Observed_PN = numeric(),
                             Percent_PN = numeric())
k = 1
for(i in 1:length(Isolations)){
  for(j in 1:length(Traits)){
    tempObs = IsolationXTraitCounts_df[which(IsolationXTraitCounts_df[,1] == Isolations[[i]] & IsolationXTraitCounts_df[,2] == Traits[[j]]),5]
    tempList = list()
    for(perm in 1:1000){
      tempDF1 = PermutedAllAssoc_list[[perm]]$Diff_PN
      COLI= which(colnames(tempDF) == Isolations[[i]])
      ROWJ= which(rownames(tempDF) == Traits[[j]])
      COLJ= which(colnames(tempDF) == Traits[[j]])
      ROWI= which(rownames(tempDF) == Isolations[[i]])
      OBS1 = tempDF1[ROWJ,COLI]
      OBS2 = tempDF1[ROWI,COLJ]
      if(OBS1 > OBS2){
        permOBS = OBS1
      }else{
        permOBS = OBS2
      }
      tempList[perm] = permOBS
    }
    temporary = unlist(tempList)
    AssocIsoTrait_PN_Neg_df[k,1] = Isolations[[i]]
    AssocIsoTrait_PN_Neg_df[k,2] = Traits[[j]]
    AssocIsoTrait_PN_Neg_df[k,3] = tempObs
    AssocIsoTrait_PN_Neg_df[k,4] = ecdf(temporary)(tempObs)
    k = k + 1
  }
}

write.csv(AssocIsoTrait_PN_Neg_df, file = paste("AssocIsoTrait_PN_Neg_df_", Sys.Date(), ".csv", sep = ""), row.names = F)

# Determine the cut-off for signficant traits in an environment by calculating the values of the data found in the 97.5th% of the distribution #
NegAssocIso_PN_split = splitBy("Isolation", AssocIsoTrait_PN_Neg_df)

NegIsoCutOffs_PN_df = data.frame(Isolation = character(), CutOff_97.5 = numeric(), CutOff_95 = numeric())

for(i in 1:length(NegAssocIso_PN_split)){
	NegIsoCutOffs_PN_df[i,1] = names(NegAssocIso_PN_split)[i]
	NegIsoCutOffs_PN_df[i,2] = quantile(NegAssocIso_PN_split[[i]][,4], 0.975)
	NegIsoCutOffs_PN_df[i,3] = quantile(NegAssocIso_PN_split[[i]][,4], 0.95)
}

write.csv(NegIsoCutOffs_PN_df, file = paste("NegIsoCutOffs_PN_df_", Sys.Date(), ".csv", sep = ""), row.names = F)

# Determine which traits are significant in an environment #
Neg_SigTraitsIso_PN_list = list()
for(i in 1:nrow(NegIsoCutOffs_PN_df)){
  Neg_SigTraitsIso_PN_list[[i]] = list()
  temp_df = AssocIsoTrait_PN_Neg_df[which(AssocIsoTrait_PN_Neg_df[,1] == NegIsoCutOffs_PN_df[i,1]),]
  Neg_SigTraitsIso_PN_list[[i]]$Isolation = NegIsoCutOffs_PN_df[i,1]
  Neg_SigTraitsIso_PN_list[[i]]$Traits = temp_df[which(temp_df[,4] >= NegIsoCutOffs_PN_df[i,2]),2]
}


# Determine the number of observed (-/+)'s for Isolation and Trait for the observed data and compare it to the permuted data # 
AssocIsoTrait_NP_Neg_df = data.frame(Isolation = character(),
                             Trait = character(),
                             Observed_NP = numeric(),
                             Percent_NP = numeric())
k = 1
for(i in 1:length(Isolations)){
  for(j in 1:length(Traits)){
    tempObs = IsolationXTraitCounts_df[which(IsolationXTraitCounts_df[,1] == Isolations[[i]] & IsolationXTraitCounts_df[,2] == Traits[[j]]),6]
    tempList = list()
    for(perm in 1:1000){
      tempDF1 = PermutedAllAssoc_list[[perm]]$Diff_NP
      COLI= which(colnames(tempDF) == Isolations[[i]])
      ROWJ= which(rownames(tempDF) == Traits[[j]])
      COLJ= which(colnames(tempDF) == Traits[[j]])
      ROWI= which(rownames(tempDF) == Isolations[[i]])
      OBS1 = tempDF1[ROWJ,COLI]
      OBS2 = tempDF1[ROWI,COLJ]
      if(OBS1 > OBS2){
        permOBS = OBS1
      }else{
        permOBS = OBS2
      }
      tempList[perm] = permOBS
    }
    temporary = unlist(tempList)
    AssocIsoTrait_NP_Neg_df[k,1] = Isolations[[i]]
    AssocIsoTrait_NP_Neg_df[k,2] = Traits[[j]]
    AssocIsoTrait_NP_Neg_df[k,3] = tempObs
    AssocIsoTrait_NP_Neg_df[k,4] = ecdf(temporary)(tempObs)
    k = k + 1
  }
}

write.csv(AssocIsoTrait_NP_Neg_df, file = paste("AssocIsoTrait_NP_Neg_df_", Sys.Date(), ".csv", sep = ""), row.names = F)

# Determine the cut-off for signficant traits in an environment by calculating the values of the data found in the 97.5th% of the distribution #
NegAssocIso_NP_split = splitBy("Isolation", AssocIsoTrait_NP_Neg_df)

NegIsoCutOffs_NP_df = data.frame(Isolation = character(), CutOff_97.5 = numeric(), CutOff_95 = numeric())

for(i in 1:length(NegAssocIso_NP_split)){
	NegIsoCutOffs_NP_df[i,1] = names(NegAssocIso_NP_split)[i]
	NegIsoCutOffs_NP_df[i,2] = quantile(NegAssocIso_NP_split[[i]][,4], 0.975)
	NegIsoCutOffs_NP_df[i,3] = quantile(NegAssocIso_NP_split[[i]][,4], 0.95)
}

write.csv(NegIsoCutOffs_NP_df, file = paste("NegIsoCutOffs_NP_df_", Sys.Date(), ".csv", sep = ""), row.names = F)

# Determine which traits are significant in an environment #
Neg_SigTraitsIso_NP_list = list()
for(i in 1:nrow(NegIsoCutOffs_NP_df)){
  Neg_SigTraitsIso_NP_list[[i]] = list()
  temp_df = AssocIsoTrait_NP_Neg_df[which(AssocIsoTrait_NP_Neg_df[,1] == NegIsoCutOffs_NP_df[i,1]),]
  Neg_SigTraitsIso_NP_list[[i]]$Isolation = NegIsoCutOffs_NP_df[i,1]
  Neg_SigTraitsIso_NP_list[[i]]$Traits = temp_df[which(temp_df[,4] >= NegIsoCutOffs_NP_df[i,2]),2]
}

#############################################################################################################
########### Calculates counts of each type of trait presence/absence [(+/+), (-/-), (+/-), (-/+)] ###########
########### within an isolation environment                                                       ###########
########### for all trait pairs where PP is (+/+), NN is (-/-), PN is (+/-) and NP is (-/+)       ###########  
#############################################################################################################

# Calculates trait pair counts for permuted data in an isolation environment #
PermutedIsoAllAssoc_list = list()
for(i in 1:1000){
  PermutedIsoAllAssoc_list[[i]] = list()
  CompiledData = PermutedTraitxTrait_IsoList[[i]]$CompiledData
  for(iso in 1:length(Isolations)){
    PermutedIsoAllAssoc_list[[i]][[iso]] = list()
    IsoCOL= which(colnames(CompiledData) == Isolations[[iso]])
    tempComp = CompiledData[which(CompiledData[,IsoCOL] == 1),c(56:ncol(CompiledData))] 
    
    tempMatrix_PP = matrix(0, ncol = ncol(tempComp), nrow = ncol(tempComp))
    colnames(tempMatrix_PP) = colnames(tempComp)
    rownames(tempMatrix_PP) = colnames(tempComp)
  
    tempMatrix_NN = matrix(0, ncol = ncol(tempComp), nrow = ncol(tempComp))
    colnames(tempMatrix_NN) = colnames(tempComp)
    rownames(tempMatrix_NN) = colnames(tempComp)
  
    tempMatrix_PN = matrix(0, ncol = ncol(tempComp), nrow = ncol(tempComp))
    colnames(tempMatrix_PN) = colnames(tempComp)
    rownames(tempMatrix_PN) = colnames(tempComp)
  
    tempMatrix_NP = matrix(0, ncol = ncol(tempComp), nrow = ncol(tempComp))
    colnames(tempMatrix_NP) = colnames(tempComp)
    rownames(tempMatrix_NP) = colnames(tempComp)
  
    for(j in 1:(ncol(tempComp)-1)){
      for(l in (j+1):ncol(tempComp)){
        POS = length(which(tempComp[,j] == 1 & tempComp[,l] == 1))
        NEG = length(which(tempComp[,j] == 0 & tempComp[,l] == 0))
        DIFF_PN = length(which(tempComp[,j] == 1 & tempComp[,l] == 0))
        DIFF_NP = length(which(tempComp[,j] == 0 & tempComp[,l] == 1))
        tempMatrix_PP[j,l] = POS
        tempMatrix_NN[j,1] = NEG
        tempMatrix_PN[j,l] = DIFF_PN
        tempMatrix_NP[j,l] = DIFF_NP
      }
    }
    PermutedIsoAllAssoc_list[[i]][[iso]]$Isolation = Isolations[[iso]]
    PermutedIsoAllAssoc_list[[i]][[iso]]$Both_P = tempMatrix_PP
    PermutedIsoAllAssoc_list[[i]][[iso]]$Both_N = tempMatrix_NN
    PermutedIsoAllAssoc_list[[i]][[iso]]$Diff_PN = tempMatrix_PN
    PermutedIsoAllAssoc_list[[i]][[iso]]$Diff_NP = tempMatrix_NP
  }
  print(i)
}

save(PermutedIsoAllAssoc_list, file = "PermutedIsoAllAssoc_list.rda")

traits$Species = rownames(traits)
envs$Species = rownames(envs)

CompiledData = merge(envs , traits, by = "Species")

# Calculates trait pair counts for observed data in an isolation environment #
IsoAllAssoc_list = list()
for(iso in 1:length(Isolations)){
  IsoAllAssoc_list[[iso]] = list()
  IsoCOL= which(colnames(CompiledData) == Isolations[[iso]])
  tempComp = CompiledData[which(CompiledData[,IsoCOL] == 1),c(57:ncol(CompiledData))] 
    
  tempMatrix_PP = matrix(0, ncol = ncol(tempComp), nrow = ncol(tempComp))
  colnames(tempMatrix_PP) = colnames(tempComp)
  rownames(tempMatrix_PP) = colnames(tempComp)
    
  tempMatrix_NN = matrix(0, ncol = ncol(tempComp), nrow = ncol(tempComp))
  colnames(tempMatrix_NN) = colnames(tempComp)
  rownames(tempMatrix_NN) = colnames(tempComp)
    
  tempMatrix_PN = matrix(0, ncol = ncol(tempComp), nrow = ncol(tempComp))
  colnames(tempMatrix_PN) = colnames(tempComp)
  rownames(tempMatrix_PN) = colnames(tempComp)
    
  tempMatrix_NP = matrix(0, ncol = ncol(tempComp), nrow = ncol(tempComp))
  colnames(tempMatrix_NP) = colnames(tempComp)
  rownames(tempMatrix_NP) = colnames(tempComp)
    
  for(j in 1:(ncol(tempComp)-1)){
    for(l in (j+1):ncol(tempComp)){
      POS = length(which(tempComp[,j] == 1 & tempComp[,l] == 1))
      NEG = length(which(tempComp[,j] == 0 & tempComp[,l] == 0))
      DIFF_PN = length(which(tempComp[,j] == 1 & tempComp[,l] == 0))
      DIFF_NP = length(which(tempComp[,j] == 0 & tempComp[,l] == 1))
      tempMatrix_PP[j,l] = POS
      tempMatrix_NN[j,1] = NEG
      tempMatrix_PN[j,l] = DIFF_PN
      tempMatrix_NP[j,l] = DIFF_NP
    }
  }
  IsoAllAssoc_list[[iso]]$Isolation = Isolations[[iso]]
  IsoAllAssoc_list[[iso]]$Both_P = tempMatrix_PP
  IsoAllAssoc_list[[iso]]$Both_N = tempMatrix_NN
  IsoAllAssoc_list[[iso]]$Diff_PN = tempMatrix_PN
  IsoAllAssoc_list[[iso]]$Diff_NP = tempMatrix_NP
  print(iso)
}

save(IsoAllAssoc_list, file = paste("IsoAllAssoc_list_", Sys.Date(), ".rda", sep = ""))

Traits = unique(colnames(PermutedIsoAllAssoc_list[[1]][[1]]$Both_P))

# Determine the number of observed (+/+)'s for Trait pairs in an isolation environment for the observed data and compare it to the permuted data # 
IsoAnalysis_PP_df = data.frame(Isolation = character(),
                                  Trait_A = character(),
                                  Trait_B = character(),
                                  Observed_PP = numeric(), 
                                  Percent_PP = numeric())

k = 1
for(iso in 1:length(Isolations)){
  tempObs = IsoAllAssoc_list[[which(IsoAllAssoc_list[[iso]]$Isolation == Isolations[[iso]])]]$Both_P
  for(i in 1:(length(Traits)-1)){
    for(j in (i+1):length(Traits)){
      temp_list = list()
      for(perm in 1:1000){
        Val1 = PermutedIsoAllAssoc_list[[perm]][[iso]]$Both_P[i,j]
        Val2 = PermutedIsoAllAssoc_list[[perm]][[iso]]$Both_P[j,i]
        if(Val1 != Val2 & Val1 > Val2){
          temp_list[perm] = Val1
        }else{ 
          temp_list[perm] = Val2
        }
      }
      templist = unlist(temp_list)
      ObsVal1 = tempObs[i,j]
      ObsVal2 = tempObs[j,1]
      if(ObsVal1 != ObsVal2 & ObsVal1 > ObsVal2){
        ObsValue = ObsVal1
      }else{
        ObsValue = ObsVal2
      }
      IsoAnalysis_PP_df[k,1] = Isolations[iso]
      IsoAnalysis_PP_df[k,2] = Traits[i]
      IsoAnalysis_PP_df[k,3] = Traits[j]
      IsoAnalysis_PP_df[k,4] = ObsValue
      IsoAnalysis_PP_df[k,5] = ecdf(templist)(ObsValue)
      k = k+1
    }
  } 
}

write.csv(IsoAnalysis_PP_df, file = paste("IsoAnalysis_PP_df_", Sys.Date(), ".csv", sep = ""), row.names = F)

# Limit Trait pairs in the isolation environment to trait pairs where one trait is significantly associated with the isolation environment - Positive (+/+) #
AllSigIsoAnalysis_PP_df = matrix(0, nrow = 1 ,ncol = ncol(IsoAnalysis_PP_df))
colnames(AllSigIsoAnalysis_PP_df) = colnames(IsoAnalysis_PP_df)

AllSigIsoAnalysis_PP_df = data.frame(AllSigIsoAnalysis_PP_df)

for(i in 1:length(Pos_SigTraitsIso_PP_list)){
  iso = Pos_SigTraitsIso_PP_list[[i]]$Isolation
  keepTraits = Pos_SigTraitsIso_PP_list[[i]]$Traits
  tempData = IsoAnalysis_PP_df[which(IsoAnalysis_PP_df$Isolation == iso),]
  dataKeep = tempData[which(tempData$Trait_A %in% keepTraits | tempData$Trait_B %in% keepTraits),]
  AllSigIsoAnalysis_PP_df = rbind(AllSigIsoAnalysis_PP_df, dataKeep)
}

write.csv(AllSigIsoAnalysis_PP_df, file = paste("AllSigIsoAnalysis_PP_df_", Sys.Date(), ".csv", sep = ""), row.names = F)

# Limit Trait pairs in the isolation environment to trait pairs where both traits are significantly associated with the isolation environment - Positive (+/+) #

BothSigIsoAnalysis_PP_df = matrix(0, nrow = 1 ,ncol = ncol(IsoAnalysis_PP_df))
colnames(BothSigIsoAnalysis_PP_df) = colnames(IsoAnalysis_PP_df)

BothSigIsoAnalysis_PP_df = data.frame(BothSigIsoAnalysis_PP_df)

for(i in 1:length(Pos_SigTraitsIso_PP_list)){
  iso = Pos_SigTraitsIso_PP_list[[i]]$Isolation
  keepTraits = Pos_SigTraitsIso_PP_list[[i]]$Traits
  tempData = IsoAnalysis_PP_df[which(IsoAnalysis_PP_df$Isolation == iso),]
  dataKeep = tempData[which(tempData$Trait_A %in% keepTraits),] 
  dataKeep = dataKeep[which(dataKeep$Trait_B %in% keepTraits),]
  BothSigIsoAnalysis_PP_df = rbind(BothSigIsoAnalysis_PP_df, dataKeep)
}

write.csv(BothSigIsoAnalysis_PP_df, file = paste("BothSigIsoAnalysis_PP_df_", Sys.Date(), ".csv", sep = ""), row.names = F)

# Determine the number of observed (-/-)'s for Trait pairs in an isolation environment for the observed data and compare it to the permuted data # 
IsoAnalysis_NN_df = data.frame(Isolation = character(),
                               Trait_A = character(),
                               Trait_B = character(),
                               Observed_NN = numeric(), 
                               Percent_NN = numeric())

k = 1
for(iso in 1:length(Isolations)){
  tempObs = IsoAllAssoc_list[[which(IsoAllAssoc_list[[iso]]$Isolation == Isolations[[iso]])]]$Both_N
  for(i in 1:(length(Traits)-1)){
    for(j in (i+1):length(Traits)){
      temp_list = list()
      for(perm in 1:1000){
        Val1 = PermutedIsoAllAssoc_list[[perm]][[iso]]$Both_P[i,j]
        Val2 = PermutedIsoAllAssoc_list[[perm]][[iso]]$Both_P[j,i]
        if(Val1 != Val2 & Val1 > Val2){
          temp_list[perm] = Val1
        }else{ 
          temp_list[perm] = Val2
        }
      }
      templist = unlist(temp_list)
      ObsVal1 = tempObs[i,j]
      ObsVal2 = tempObs[j,1]
      if(ObsVal1 != ObsVal2 & ObsVal1 > ObsVal2){
        ObsValue = ObsVal1
      }else{
        ObsValue = ObsVal2
      }
      IsoAnalysis_NN_df[k,1] = Isolations[iso]
      IsoAnalysis_NN_df[k,2] = Traits[i]
      IsoAnalysis_NN_df[k,3] = Traits[j]
      IsoAnalysis_NN_df[k,4] = ObsValue
      IsoAnalysis_NN_df[k,5] = ecdf(templist)(ObsValue)
      k = k+1
    }
  } 
}

write.csv(IsoAnalysis_NN_df, file = paste("IsoAnalysis_NN_df_", Sys.Date(), ".csv", sep = ""), row.names = F)

# Limit Trait pairs in the isolation environment to trait pairs where one trait is significantly associated with the isolation environment - Positive (-/-) #
AllSigIsoAnalysis_NN_df = matrix(0, nrow = 1 ,ncol = ncol(IsoAnalysis_NN_df))
colnames(AllSigIsoAnalysis_NN_df) = colnames(IsoAnalysis_NN_df)

AllSigIsoAnalysis_NN_df = data.frame(AllSigIsoAnalysis_NN_df)

for(i in 1:length(Pos_SigTraitsIso_NN_list)){
  iso = Pos_SigTraitsIso_NN_list[[i]]$Isolation
  keepTraits = Pos_SigTraitsIso_NN_list[[i]]$Traits
  tempData = IsoAnalysis_NN_df[which(IsoAnalysis_NN_df$Isolation == iso),]
  dataKeep = tempData[which(tempData$Trait_A %in% keepTraits | tempData$Trait_B %in% keepTraits),]
  AllSigIsoAnalysis_NN_df = rbind(AllSigIsoAnalysis_NN_df, dataKeep)
}

write.csv(AllSigIsoAnalysis_NN_df, file = paste("AllSigIsoAnalysis_NN_df_", Sys.Date(), ".csv", sep = ""), row.names = F)

# Limit Trait pairs in the isolation environment to trait pairs where both traits are significantly associated with the isolation environment - Positive (-/-) #
BothSigIsoAnalysis_NN_df = matrix(0, nrow = 1 ,ncol = ncol(IsoAnalysis_NN_df))
colnames(BothSigIsoAnalysis_NN_df) = colnames(IsoAnalysis_NN_df)

BothSigIsoAnalysis_NN_df = data.frame(BothSigIsoAnalysis_NN_df)

for(i in 1:length(Pos_SigTraitsIso_NN_list)){
  iso = Pos_SigTraitsIso_NN_list[[i]]$Isolation
  keepTraits = Pos_SigTraitsIso_NN_list[[i]]$Traits
  tempData = IsoAnalysis_NN_df[which(IsoAnalysis_NN_df$Isolation == iso),]
  dataKeep = tempData[which(tempData$Trait_A %in% keepTraits),] 
  dataKeep = dataKeep[which(dataKeep$Trait_B %in% keepTraits),]
  BothSigIsoAnalysis_NN_df = rbind(BothSigIsoAnalysis_NN_df, dataKeep)
}

write.csv(BothSigIsoAnalysis_NN_df, file = paste("BothSigIsoAnalysis_NN_df_", Sys.Date(), ".csv", sep = ""), row.names = F)

# Determine the number of observed (+/-)'s for Trait pairs in an isolation environment for the observed data and compare it to the permuted data # 
IsoAnalysis_PN_df = data.frame(Isolation = character(),
                               Trait_A = character(),
                               Trait_B = character(),
                               Observed_PN = numeric(), 
                               Percent_PN = numeric())

k = 1
for(iso in 1:length(Isolations)){
  tempObs = IsoAllAssoc_list[[which(IsoAllAssoc_list[[iso]]$Isolation == Isolations[[iso]])]]$Diff_PN
  for(i in 1:(length(Traits)-1)){
    for(j in (i+1):length(Traits)){
      temp_list = list()
      for(perm in 1:1000){
        Val1 = PermutedIsoAllAssoc_list[[perm]][[iso]]$Both_P[i,j]
        Val2 = PermutedIsoAllAssoc_list[[perm]][[iso]]$Both_P[j,i]
        if(Val1 != Val2 & Val1 > Val2){
          temp_list[perm] = Val1
        }else{ 
          temp_list[perm] = Val2
        }
      }
      templist = unlist(temp_list)
      ObsVal1 = tempObs[i,j]
      ObsVal2 = tempObs[j,1]
      if(ObsVal1 != ObsVal2 & ObsVal1 > ObsVal2){
        ObsValue = ObsVal1
      }else{
        ObsValue = ObsVal2
      }
      IsoAnalysis_PN_df[k,1] = Isolations[iso]
      IsoAnalysis_PN_df[k,2] = Traits[i]
      IsoAnalysis_PN_df[k,3] = Traits[j]
      IsoAnalysis_PN_df[k,4] = ObsValue
      IsoAnalysis_PN_df[k,5] = ecdf(templist)(ObsValue)
      k = k+1
    }
  } 
}

write.csv(IsoAnalysis_PN_df, file = paste("IsoAnalysis_PN_df_", Sys.Date(), ".csv", sep = ""), row.names = F)

# Limit Trait pairs in the isolation environment to trait pairs where one trait is significantly associated with the isolation environment - Negative (+/-) #
AllSigIsoAnalysis_PN_df = matrix(0, nrow = 1 ,ncol = ncol(IsoAnalysis_PN_df))
colnames(AllSigIsoAnalysis_PN_df) = colnames(IsoAnalysis_PN_df)

AllSigIsoAnalysis_PN_df = data.frame(AllSigIsoAnalysis_PN_df)

for(i in 1:length(Neg_SigTraitsIso_PN_list)){
  iso = Neg_SigTraitsIso_PN_list[[i]]$Isolation
  keepTraits = Neg_SigTraitsIso_PN_list[[i]]$Traits
  tempData = IsoAnalysis_PN_df[which(IsoAnalysis_PN_df$Isolation == iso),]
  dataKeep = tempData[which(tempData$Trait_A %in% keepTraits | tempData$Trait_B %in% keepTraits),]
  AllSigIsoAnalysis_PN_df = rbind(AllSigIsoAnalysis_PN_df, dataKeep)
}

write.csv(AllSigIsoAnalysis_PN_df, file = paste("AllSigIsoAnalysis_PN_df_", Sys.Date(), ".csv", sep = ""), row.names = F)

# Limit Trait pairs in the isolation environment to trait pairs where both traits are significantly associated with the isolation environment - Negative (+/-) #
BothSigIsoAnalysis_PN_df = matrix(0, nrow = 1 ,ncol = ncol(IsoAnalysis_PN_df))
colnames(BothSigIsoAnalysis_PN_df) = colnames(IsoAnalysis_PN_df)

BothSigIsoAnalysis_PN_df = data.frame(BothSigIsoAnalysis_PN_df)

for(i in 1:length(Neg_SigTraitsIso_PN_list)){
  iso = Neg_SigTraitsIso_PN_list[[i]]$Isolation
  keepTraits = Neg_SigTraitsIso_PN_list[[i]]$Traits
  tempData = IsoAnalysis_PN_df[which(IsoAnalysis_PN_df$Isolation == iso),]
  dataKeep = tempData[which(tempData$Trait_A %in% keepTraits),] 
  dataKeep = dataKeep[which(dataKeep$Trait_B %in% keepTraits),]
  BothSigIsoAnalysis_PN_df = rbind(BothSigIsoAnalysis_PN_df, dataKeep)
}

write.csv(BothSigIsoAnalysis_PN_df, file = paste("BothSigIsoAnalysis_PN_df_", Sys.Date(), ".csv", sep = ""), row.names = F)

# Determine the number of observed (-/+)'s for Trait pairs in an isolation environment for the observed data and compare it to the permuted data # 
IsoAnalysis_NP_df = data.frame(Isolation = character(),
                               Trait_A = character(),
                               Trait_B = character(),
                               Observed_NP = numeric(), 
                               Percent_NP = numeric())

k = 1
for(iso in 1:length(Isolations)){
  tempObs = IsoAllAssoc_list[[which(IsoAllAssoc_list[[iso]]$Isolation == Isolations[[iso]])]]$Diff_NP
  for(i in 1:(length(Traits)-1)){
    for(j in (i+1):length(Traits)){
      temp_list = list()
      for(perm in 1:1000){
        Val1 = PermutedIsoAllAssoc_list[[perm]][[iso]]$Both_P[i,j]
        Val2 = PermutedIsoAllAssoc_list[[perm]][[iso]]$Both_P[j,i]
        if(Val1 != Val2 & Val1 > Val2){
          temp_list[perm] = Val1
        }else{ 
          temp_list[perm] = Val2
        }
      }
      templist = unlist(temp_list)
      ObsVal1 = tempObs[i,j]
      ObsVal2 = tempObs[j,1]
      if(ObsVal1 != ObsVal2 & ObsVal1 > ObsVal2){
        ObsValue = ObsVal1
      }else{
        ObsValue = ObsVal2
      }
      IsoAnalysis_NP_df[k,1] = Isolations[iso]
      IsoAnalysis_NP_df[k,2] = Traits[i]
      IsoAnalysis_NP_df[k,3] = Traits[j]
      IsoAnalysis_NP_df[k,4] = ObsValue
      IsoAnalysis_NP_df[k,5] = ecdf(templist)(ObsValue)
      k = k+1
    }
  } 
}

write.csv(IsoAnalysis_NP_df, file = paste("IsoAnalysis_NP_df_", Sys.Date(), ".csv", sep = ""), row.names = F)

# Limit Trait pairs in the isolation environment to trait pairs where one trait is significantly associated with the isolation environment - Negative (-/+) #
AllSigIsoAnalysis_NP_df = matrix(0, nrow = 1 ,ncol = ncol(IsoAnalysis_NP_df))
colnames(AllSigIsoAnalysis_NP_df) = colnames(IsoAnalysis_NP_df)

AllSigIsoAnalysis_NP_df = data.frame(AllSigIsoAnalysis_NP_df)

for(i in 1:length(Neg_SigTraitsIso_NP_list)){
  iso = Neg_SigTraitsIso_NP_list[[i]]$Isolation
  keepTraits = Neg_SigTraitsIso_NP_list[[i]]$Traits
  tempData = IsoAnalysis_NP_df[which(IsoAnalysis_NP_df$Isolation == iso),]
  dataKeep = tempData[which(tempData$Trait_A %in% keepTraits | tempData$Trait_B %in% keepTraits),]
  AllSigIsoAnalysis_NP_df = rbind(AllSigIsoAnalysis_NP_df, dataKeep)
}

write.csv(AllSigIsoAnalysis_NP_df, file = paste("AllSigIsoAnalysis_NP_df_", Sys.Date(), ".csv", sep = ""), row.names = F)

# Limit Trait pairs in the isolation environment to trait pairs where both traits are significantly associated with the isolation environment - Negative (-/+) #
BothSigIsoAnalysis_NP_df = matrix(0, nrow = 1 ,ncol = ncol(IsoAnalysis_NP_df))
colnames(BothSigIsoAnalysis_NP_df) = colnames(IsoAnalysis_NP_df)

BothSigIsoAnalysis_NP_df = data.frame(BothSigIsoAnalysis_NP_df)

for(i in 1:length(Neg_SigTraitsIso_NP_list)){
  iso = Neg_SigTraitsIso_NP_list[[i]]$Isolation
  keepTraits = Neg_SigTraitsIso_NP_list[[i]]$Traits
  tempData = IsoAnalysis_NP_df[which(IsoAnalysis_NP_df$Isolation == iso),]
  dataKeep = tempData[which(tempData$Trait_A %in% keepTraits),] 
  dataKeep = dataKeep[which(dataKeep$Trait_B %in% keepTraits),]
  BothSigIsoAnalysis_NP_df = rbind(BothSigIsoAnalysis_NP_df, dataKeep)
}

write.csv(BothSigIsoAnalysis_NP_df, file = paste("BothSigIsoAnalysis_NP_df_", Sys.Date(), ".csv", sep = ""), row.names = F)
