options(stringsAsFactors = FALSE)

require(FD)
require(picante)
require(cooccur)

traits<-read.csv("TraitMatrix_Imputed_2016-06-27_pruned.csv", check.names=FALSE, row.names=1)
envs<-read.csv("isomatrix_pruned.csv", check.names=FALSE, row.names=1)
envt<-as.matrix(t(envs))
#txt_isodf <- read.csv("TraitxTrait_IsoDF.csv", header = TRUE, check.names = F)
#txi_df <- read.csv("IsolationTrait_CoOc.csv", header = TRUE, check.names = F)
load("PermutedTraitxTrait_IsoList")
load("PermIsoSubCoOccur.rda")

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
#
############## This is to generate permuted data - Permuted Data gets loaded above ############## 
# PermIsoSubCoOccur = list()
# for(iter in 1:1000){ 
#   CompiledData = PermutedTraitxTrait_IsoList[[iter]]$CompiledData
#   PermIsoSubCoOccur[[iter]] = list()
#   for(i in 1:length(Isolations)){	
#     TempCol = which(colnames(CompiledData) == Isolations[i])
#     TempData = CompiledData[which(CompiledData[,TempCol] == 1), 56:ncol(CompiledData)]
#     PermIsoSubCoOccur[[iter]][[i]] = list()
#     PermIsoSubCoOccur[[iter]][[i]]$Isolation = Isolations[i]   
#     TraitByTraitIso_co = cooccur(mat = t(TempData), type = "spp_site", thresh = T, spp_names = T, eff_matrix = T, only_effects = F)
#     PermIsoSubCoOccur[[iter]][[i]]$CoOcc = TraitByTraitIso_co[[2]]
#   }
#   print(iter)
# }
# 
# save(PermIsoSubCoOccur, file = "PermIsoSubCoOccur.rda")



# PermutedAssoc_df = data.frame(Permutation = numeric(), Isolation = character(), Positive = numeric(), Negative = numeric(), Random = numeric())
# k = 1
# for(i in 1:1000){
#   for(j in 1:length(PermIsoSubCoOccur[[i]])){
#     PermutedAssoc_df[k,1] = i
#     PermutedAssoc_df[k,2] = PermIsoSubCoOccur[[i]][[j]]$Isolation
#     tempDF = data.frame(PermIsoSubCoOccur[[i]][[j]]$CoOcc)
#     POS = length(which(tempDF$p_gt <= 0.05))
#     NEG = length(which(tempDF$p_lt <= 0.05))
#     PermutedAssoc_df[k,3] = POS
#     PermutedAssoc_df[k,4] = NEG
#     PermutedAssoc_df[k,5] = nrow(tempDF) - (POS + NEG)
#     k = k + 1
#   }
#   print(i)
# }
# write.csv(PermutedAssoc_df, file = paste("PermutedAssoc_df_", Sys.Date(), ".csv", sep = ""), row.names = F)


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

IsoAnalysis_NN_df = data.frame(Isolation = character(),
                               Trait_A = character(),
                               Trait_B = character(),
                               Observed_NN = numeric(), 
                               Percent_NN = numeric())

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

IsoAnalysis_PN_df = data.frame(Isolation = character(),
                               Trait_A = character(),
                               Trait_B = character(),
                               Observed_PN = numeric(), 
                               Percent_PN = numeric())

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

IsoAnalysis_NP_df = data.frame(Isolation = character(),
                               Trait_A = character(),
                               Trait_B = character(),
                               Observed_NP = numeric(), 
                               Percent_NP = numeric())

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
