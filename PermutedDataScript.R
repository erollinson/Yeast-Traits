options(stringsAsFactors = FALSE)

require(FD)
require(picante)
require(cooccur)

traits<-read.csv("TraitMatrix_Imputed_2016-06-27_pruned.csv", check.names=FALSE, row.names=1)
envs<-read.csv("isomatrix_pruned.csv", check.names=FALSE, row.names=1)
envt<-as.matrix(t(envs))
txt_isodf <- read.csv("TraitxTrait_IsoDF.csv", header = TRUE, check.names = F)
txi_df <- read.csv("IsolationTrait_CoOc.csv", header = TRUE, check.names = F)

envs2 = envs
envs2$Species = rownames(envs2)

options(stringsAsFactors = FALSE)

require(FD)
require(picante)
require(cooccur)

traits<-read.csv("TraitMatrix_Imputed_2016-06-27_pruned.csv", check.names=FALSE, row.names=1)
envs<-read.csv("isomatrix_pruned.csv", check.names=FALSE, row.names=1)
envt<-as.matrix(t(envs))
#txt_isodf <- read.csv("TraitxTrait_IsoDF.csv", header = TRUE, check.names = F)
#txi_df <- read.csv("IsolationTrait_CoOc.csv", header = TRUE, check.names = F)

envs2 = envs
envs2$Species = rownames(envs2)

Isolations <- rownames(envt)
PermutedTraitxTrait_list= list()
for(iter in 1:1000){
  
  PermutedTraitxTrait_IsoList[[iter]] = list()
  tempMatrix <- randomizeMatrix(traits, null.model="independentswap", iterations=1000)
  tempMatrix = data.frame(tempMatrix, check.names = F)
  
  #CWM<-functcomp(tempMatrix, envt, CWM.type="all")
  #PermutedTraitxTrait_IsoList[[iter]]$CWM <- CWM 
  
  tempMatrix$Species = rownames(tempMatrix)
  CompiledData = merge(envs2, tempMatrix, by = "Species")
  PermutedTraitxTrait_IsoList[[iter]]$Permutation <- iter 
  CompiledData= CompiledData[,-which(colnames(CompiledData) == "Species")]
  PermutedTraitxTrait_IsoList[[iter]]$Permutation = iter
  PermutedTraitxTrait_IsoList[[iter]]$CompiledData <- CompiledData 
  print(iter)
}

save(PermutedTraitxTrait_IsoList, file = paste("PermutedTraitxTrait_IsoList_", Sys.Date(), ".rda", sep = ""))

# for(iter in 1:1000){ 
#   for(i in 1:length(Isolations)){
#     TraitIsolation_m = matrix(0, ncol = length(CompiledData[56:ncol(CompiledData)]))
#     colnames(TraitIsolation_m) = colnames(CompiledData)[56:ncol(CompiledData)]
#     TraitIsolation_m = data.frame(TraitIsolation_m, check.names = F)
#     TempCol = which(colnames(CompiledData) == Isolations[i])
#     TempData = CompiledData[which(CompiledData[,TempCol] == 1), 56:ncol(CompiledData)]
#     PermutedTraitxTrait_IsoList[[iter]][[i]] = list()
#     PermutedTraitxTrait_IsoList[[iter]][[i]]$Isolation = Isolations[i]   
#     TraitByTraitIso_co = cooccur(mat = t(TempData), type = "spp_site", thresh = T, spp_names = T, eff_matrix = T, only_effects = F)
#     PermutedTraitxTrait_IsoList[[iter]][[i]]$CoOcc = TraitByTraitIso_co[[2]]
#   }
#   print(iter)
# }

#save(PermutedTraitxTrait_IsoList, file = "PermutedData.rda")

PermutedAssoc_df = data.frame(Permutation = numeric(), Isolation = character(), Positive = numeric(), Negative = numeric(), Random = numeric())
k = 1
for(i in 1:1000){
  for(j in 1:length(PermutedTraitxTrait_IsoList[[i]]$TraitxTrait_IsoList)){
    PermutedAssoc_df[k,1] = i
    PermutedAssoc_df[k,2] = PermutedTraitxTrait_IsoList[[i]]$TraitxTrait_IsoList[[j]]$Isolation
    tempDF = data.frame(PermutedTraitxTrait_IsoList[[i]]$TraitxTrait_IsoList[[j]]$CoOcc)
    POS = length(which(tempDF$p_gt <= 0.05))
    NEG = length(which(tempDF$p_lt <= 0.05))
    PermutedAssoc_df[k,3] = POS
    PermutedAssoc_df[k,4] = NEG
    PermutedAssoc_df[k,5] = nrow(tempDF) - (POS + NEG)
    k = k + 1
  }
  print(i)
}
write.csv(PermutedAssoc_df, file = paste("PermutedAssoc_df_", Sys.Date(), ".csv", sep = ""), row.names = F)


PermutedAllAssoc_list = list()
k = 1
for(i in 1:1000){
  PermutedAllAssoc_list[[i]] = list()
  tempDF = PermutedTraitxTrait_IsoList[[i]]$TraitxTrait_IsoList[[1]]$TraitMatrix
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
k = 1
for(i in 1:1000){
  PermutedIsoAllAssoc_list[[i]] = list()
  tempDF = data.frame(PermutedTraitxTrait_IsoList[[i]]$TraitxTrait_IsoList[[1]]$TraitMatrix)
  tempDF = data.frame(tempDF)
  rownames(tempDF) = rownames(traits)
  tempDF$Species = rownames(tempDF)
  CompiledData = merge(envs2, tempDF, by = "Species")
  CompiledData= CompiledData[,-which(colnames(CompiledData) == "Species")]
  for(iso in 1:length(Isolations))
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
        DIFF_NP = length(which(tempDF[,j] == 0 & tempDF[,l] == 1))
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

