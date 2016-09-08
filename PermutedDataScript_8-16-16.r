############################################################
########### Options, Packages, and Files to Load ###########  
############################################################

# Options #
options(stringsAsFactors = FALSE)

# Libraries #
require(FD)
require(picante)


require(cooccur)
require(ggplot2)

# Files #
traits<-read.csv("TraitMatrix_Imputed_2016-06-27_pruned.csv", check.names=FALSE, row.names=1)
envs<-read.csv("isomatrix_pruned.csv", check.names=FALSE, row.names=1)
envt<-as.matrix(t(envs))
load("PermutedTraitxTrait_IsoList")

envs2 = envs
envs2$Species = rownames(envs2)

Isolations <- rownames(envt)


############## Calculate the number of species isolated from an environment ##############

IsolationCounts_df = data.frame(Isolations = character(), SpeciesCount = numeric())

for(i in 1:length(Isolations)){
	tempData = envs2[,which(colnames(envs2) == Isolations[[i]])]
	IsolationCounts_df[i,1] = Isolations[i]
	IsolationCounts_df[i,2] = length(which(tempData > 0))
}
write.csv(IsolationCounts_df, file = paste("IsolationCounts_df_", Sys.Date(), ".csv", sep = ""), row.names = F)

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

# Determine the number of observed (+/+)'s for Trait pairs for the observed data and compare it to the permuted data # 
Assoc_PP_Pos_df = data.frame(Trait_A = character(),
                             Trait_B = character(),
                             Observed_PP = numeric(),
                             Expected_PP = numeric(),
                             Percent_PP = numeric())
k = 1
for(i in 1:(length(Traits)-1)){
  for(j in (i +1):length(Traits)){
    tempObs = TraitCounts_df[which(TraitCounts_df[,1] == Traits[[i]] & TraitCounts_df[,2] == Traits[[j]]),3]
    tempList = list()
    for(perm in 1:1000){
      tempDF1 = PermutedAllAssoc_list[[perm]]$Both_P
      COLI= which(colnames(tempDF1) == Traits[[i]])
      ROWJ= which(rownames(tempDF1) == Traits[[j]])
      COLJ= which(colnames(tempDF1) == Traits[[j]])
      ROWI= which(rownames(tempDF1) == Traits[[i]])
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
    Assoc_PP_Pos_df[k,4] = mean(temporary, na.rm = T)
    Assoc_PP_Pos_df[k,5] = ecdf(temporary)(tempObs)
    k = k + 1
  }
}

write.csv(Assoc_PP_Pos_df, file = paste("Assoc_PP_Pos_df_", Sys.Date(), ".csv", sep = ""), row.names = F)

SigCutoff_PP_97 = quantile(Assoc_PP_Pos_df[,5], 0.975)
SigCutoff_PP_95 = quantile(Assoc_PP_Pos_df[,5], 0.95)


Assoc_PP_Pos_df$TraitDiff = Assoc_PP_Pos_df[,3] - Assoc_PP_Pos_df[,4]
Assoc_PP_Pos_df$ID = c(1:nrow(Assoc_PP_Pos_df))
AvgDiff_PP = mean(Assoc_PP_Pos_df$TraitDiff, na.rm = T)

a = ggplot(Assoc_PP_Pos_df, aes(x = ID, y = TraitDiff))+geom_point(pch = 19, color = "grey70")
a = a + geom_hline(yintercept = AvgDiff_PP, color = "red", linetype  = 2)
a = a + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"),axis.text.x = element_text(angle = 90))  
a = a + xlab("Traits") + ylab("Observed - Expected")
ggsave(paste("Assoc_PP_Pos_", Sys.Date(), ".pdf"))


# Determine the number of observed (-/-)'s for Trait pairs for the observed data and compare it to the permuted data # 
Assoc_NN_Pos_df = data.frame(Trait_A = character(),
                             Trait_B = character(),
                             Observed_NN = numeric(),
                             Expected_NN = numeric(), 
                             Percent_NN = numeric())
k = 1
for(i in 1:(length(Traits)-1)){
  for(j in (i +1):length(Traits)){
    tempObs = TraitCounts_df[which(TraitCounts_df[,1] == Traits[[i]] & TraitCounts_df[,2] == Traits[[j]]),4]
    tempList = list()
    for(perm in 1:1000){
      tempDF1 = PermutedAllAssoc_list[[perm]]$Both_N
      COLI= which(colnames(tempDF1) == Traits[[i]])
      ROWJ= which(rownames(tempDF1) == Traits[[j]])
      COLJ= which(colnames(tempDF1) == Traits[[j]])
      ROWI= which(rownames(tempDF1) == Traits[[i]])
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
    Assoc_NN_Pos_df[k,4] = mean(temporary, na.rm = T)
    Assoc_NN_Pos_df[k,5] = ecdf(temporary)(tempObs)
    k = k + 1
  }
}

write.csv(Assoc_NN_Pos_df, file = paste("Assoc_NN_Pos_df_", Sys.Date(), ".csv", sep = ""), row.names = F)

SigCutoff_NN_97 = quantile(Assoc_NN_Pos_df[,5], 0.975)
SigCutoff_NN_95 = quantile(Assoc_NN_Pos_df[,5], 0.95)

Assoc_NN_Pos_df$TraitDiff = Assoc_NN_Pos_df[,3] - Assoc_NN_Pos_df[,4]
Assoc_NN_Pos_df$ID = c(1:nrow(Assoc_NN_Pos_df))
AvgDiff_NN = mean(Assoc_NN_Pos_df$TraitDiff, na.rm = T)

a = ggplot(Assoc_NN_Pos_df, aes(x = ID, y = TraitDiff))+geom_point(pch = 19, color = "grey70")
a = a + geom_hline(yintercept = AvgDiff_PP, color = "red", linetype  = 2)
a = a + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"),axis.text.x = element_text(angle = 90))  
a = a + xlab("Traits") + ylab("Observed - Expected")
ggsave(paste("Assoc_NN_Pos_", Sys.Date(), ".pdf"))


Traits = unique(Assoc_NN_Pos_df$Trait_A)

AsscMatrix_m = matrix(0, nrow = length(Traits), ncol = length(Traits))
colnames(AsscMatrix_m) = Traits
rownames(AsscMatrix_m) = Traits

AsscData_df = data.frame(Trait_A = character(), Trait_B = character(), Observed = numeric(), Expected = numeric(), Difference = numeric(), Percent = numeric(),ID = numeric())
k = 1
for(i in 1:(length(Traits)-1)){
  for(j in (i+1):length(Traits)){
    tempNN = Assoc_NN_Pos_df[which(Assoc_NN_Pos_df$Trait_A == Traits[i] | Assoc_NN_Pos_df$Trait_B == Traits[i]),] 
    tempNN = tempNN[which(tempNN$Trait_A == Traits[j] | tempNN$Trait_B == Traits[j]), ]
    
    tempPP = Assoc_PP_Pos_df[which(Assoc_PP_Pos_df$Trait_A == Traits[i] | Assoc_PP_Pos_df$Trait_B == Traits[i]),] 
    tempPP = tempPP[which(tempPP$Trait_A == Traits[j] | tempPP$Trait_B == Traits[j]), ]
    
    Observed = tempPP$Observed_PP + tempNN$Observed_NN
    Expected = tempPP$Expected_PP + tempNN$Expected_NN
  
    Difference = Observed - Expected
    
    PercentNN = tempNN$Percent_NN
    PercentPP = tempPP$Percent_PP
    
    AsscMatrix_m[i,j] = Difference
    
    AsscData_df[k,1] = Traits[i]
    AsscData_df[k,2] = Traits[j]
    AsscData_df[k,3] = Observed
    AsscData_df[k,4] = Expected
    AsscData_df[k,5] = Difference
    
    if(PercentPP == 1 | PercentNN == 1 & Difference > 0){
      AsscData_df[k,6] = 1  
    }else if(PercentPP == 1 | PercentNN == 1 & Difference < 0){
      AsscData_df[k,6] = -1
    }else{
      AsscData_df[k,6] = 0
    }
    
    AsscData_df[k,7] = k
    k = k+1
  }
}



COLORS = c("#766F54","#9F8351","white","#A53010" ,"#DE7E18", "#728653", "#6AAC91")
pdf(paste("AssocHeatMap.pdf", Sys.Date(),".pdf",sep = ""), height = 15, width = 15)
pheatmap(AsscMatrix_m, cluster_rows = F,cluster_cols = F,cellheight = 10, cellwidth = 10, border_color = "black",color = colorRampPalette(COLORS)(100))
dev.off()

MeanDiff = mean(AsscData_df$Difference, na.rm = T)

a = ggplot(AsscData_df, aes(x = ID, y = Difference))+geom_point(color = "grey70", pch = 16 ,fill = "black")
a = a + geom_hline(yintercept = MeanDiff, color = "red", lty = 2)
a = a + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"),axis.text.x = element_text(angle = 90, size = 0, color = "white"))  
a = a + xlab("Traits") + ylab("Observed - Expected")
ggsave(paste("AsscData_", Sys.Date(), ".pdf"))

###### Determine negative associations between Trait pairs ###### 

# Determine the number of observed (+/-)'s for Trait pairs for the observed data and compare it to the permuted data # 
Assoc_PN_Neg_df = data.frame(Trait_A = character(),
                             Trait_B = character(),
                             Observed_PN = numeric(),
                             Expected_PN = numeric(),
                             Percent_PN = numeric())
k = 1
for(i in 1:(length(Traits)-1)){
  for(j in (i +1):length(Traits)){
    tempObs = TraitCounts_df[which(TraitCounts_df[,1] == Traits[[i]] & TraitCounts_df[,2] == Traits[[j]]),5]
    tempList = list()
    for(perm in 1:1000){
      tempDF1 = PermutedAllAssoc_list[[perm]]$Diff_PN
      COLI= which(colnames(tempDF1) == Traits[[i]])
      ROWJ= which(rownames(tempDF1) == Traits[[j]])
      COLJ= which(colnames(tempDF1) == Traits[[j]])
      ROWI= which(rownames(tempDF1) == Traits[[i]])
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
    Assoc_PN_Neg_df[k,4] = mean(temporary, na.rm = T)
    Assoc_PN_Neg_df[k,5] = ecdf(temporary)(tempObs)
    k = k + 1
  }
}

write.csv(Assoc_PN_Neg_df, file = paste("Assoc_PN_Neg_df_", Sys.Date(), ".csv", sep = ""), row.names = F)

SigCutoff_PN_97 = quantile(Assoc_PN_Neg_df[,5], 0.975)
SigCutoff_PN_95 = quantile(Assoc_PN_Neg_df[,5], 0.95)

Assoc_PN_Neg_df$TraitDiff = Assoc_PN_Neg_df[,3] - Assoc_PN_Neg_df[,4]
Assoc_PN_Neg_df$ID = c(1:nrow(Assoc_PN_Neg_df))
AvgDiff_PN = mean(Assoc_PN_Neg_df$TraitDiff, na.rm = T)

a = ggplot(Assoc_PN_Neg_df, aes(x = ID, y = TraitDiff))+geom_point(pch = 19, color = "grey70")
a = a + geom_hline(yintercept = AvgDiff_PP, color = "red", linetype  = 2)
a = a + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"),axis.text.x = element_text(angle = 90))  
a = a + xlab("Traits") + ylab("Observed - Expected")
ggsave(paste("Assoc_PN_Neg_", Sys.Date(), ".pdf"))

# Determine the number of observed (-/+)'s for Trait pairs for the observed data and compare it to the permuted data # 
Assoc_NP_Neg_df = data.frame(Trait_A = character(),
                             Trait_B = character(),
                             Observed_NP = numeric(),
                             Expected_NP = numeric(),
                             Percent_NP = numeric())
k = 1
for(i in 1:(length(Traits)-1)){
  for(j in (i +1):length(Traits)){
    tempObs = TraitCounts_df[which(TraitCounts_df[,1] == Traits[[i]] & TraitCounts_df[,2] == Traits[[j]]),6]
    tempList = list()
    for(perm in 1:1000){
      tempDF1 = PermutedAllAssoc_list[[perm]]$Diff_NP
      COLI= which(colnames(tempDF1) == Traits[[i]])
      ROWJ= which(rownames(tempDF1) == Traits[[j]])
      COLJ= which(colnames(tempDF1) == Traits[[j]])
      ROWI= which(rownames(tempDF1) == Traits[[i]])
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
    Assoc_NP_Neg_df[k,4] = mean(temporary, na.rm = T)
    Assoc_NP_Neg_df[k,5] = ecdf(temporary)(tempObs)
    k = k + 1
  }
}

write.csv(Assoc_NP_Neg_df, file = paste("Assoc_NP_Neg_df_", Sys.Date(), ".csv", sep = ""), row.names = F)

SigCutoff_NP_97 = quantile(Assoc_NP_Neg_df[,5], 0.975)
SigCutoff_NP_95 = quantile(Assoc_NP_Neg_df[,5], 0.95)

Assoc_NP_Neg_df$TraitDiff = Assoc_NP_Neg_df[,3] - Assoc_NP_Neg_df[,4]
Assoc_NP_Neg_df$ID = c(1:nrow(Assoc_NP_Neg_df))
AvgDiff_NP = mean(Assoc_NP_Neg_df$TraitDiff, na.rm = T)

a = ggplot(Assoc_NP_Neg_df, aes(x = ID, y = TraitDiff))+geom_point(pch = 19, color = "grey70")
a = a + geom_hline(yintercept = AvgDiff_PP, color = "red", linetype  = 2)
a = a + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"),axis.text.x = element_text(angle = 90, size = 0))  
a = a + xlab("Traits") + ylab("Observed - Expected")
ggsave(paste("Assoc_NP_Neg_", Sys.Date(), ".pdf"))

Neg_AsscMatrix_m = matrix(0, nrow = length(Traits), ncol = length(Traits))
colnames(Neg_AsscMatrix_m) = Traits
rownames(Neg_AsscMatrix_m) = Traits

NegAsscData_df = data.frame(Trait_A = character(), Trait_B = character(), Observed = numeric(), Expected = numeric(), Difference = numeric(), Percent = numeric(),ID = numeric())
k = 1
for(i in 1:length(Traits)){
  for(j in 1:length(Traits)){
    tempNP = Assoc_NP_Neg_df[which(Assoc_NP_Neg_df$Trait_A == Traits[i]),] 
    tempNP = tempNP[which(tempNP$Trait_B == Traits[j]), ]
    
    tempPN = Assoc_PN_Neg_df[which(Assoc_PN_Neg_df$Trait_A == Traits[i]),] 
    tempPN = tempPN[which(tempPN$Trait_B == Traits[j]), ]
    
    Observed = tempPN$Observed_PN + tempNP$Observed_NP  
    Expected = tempPN$Expected_PN + tempNP$Expected_NP
    
    Difference = Observed - Expected
    
    PercentNP = tempNP$Percent_NP
    PercentPN = tempPN$Percent_PN
    
    Neg_AsscMatrix_m[i,j] = Difference
    
    NegAsscData_df[k,1] = Traits[i]
    NegAsscData_df[k,2] = Traits[j]
    NegAsscData_df[k,3] = Observed
    NegAsscData_df[k,4] = Expected
    NegAsscData_df[k,5] = Difference
    
    if(PercentPN == 1 | PercentNP == 1 & Difference > 0){
      NegAsscData_df[k,6] = -1  
    }else if(PercentPN == 1 | PercentNP == 1 & Difference < 0){
      NegAsscData_df[k,6] = -1
    }else{
      NegAsscData_df[k,6] = 0
    }
    
    NegAsscData_df[k,7] = k
    k = k+1
  }
}


#COLORS = c("#766F54","#9F8351","white","#A53010" ,"#DE7E18", "#728653", "#6AAC91")
#pdf(paste("AssocHeatMap.pdf", Sys.Date(),".pdf",sep = ""), height = 15, width = 15)
#pheatmap(AsscMatrix_m, cluster_rows = F,cluster_cols = F,cellheight = 10, cellwidth = 10, border_color = "black",color = colorRampPalette(COLORS)(100))
#dev.off()


PooledAsscData_df = rbind(AsscData_df, NegAsscData_df)
AllAsscData = matrix(0, nrow = 1, ncol = ncol(PooledAsscData_df))
colnames(AllAsscData) = colnames(PooledAsscData_df)
AllAsscData = data.frame(AllAsscData)

k = 1

for(i in 1:(length(Traits)-1)){
  for(j in (i+1):length(Traits)){
    tempData = PooledAsscData_df[which(PooledAsscData_df$Trait_A == Traits[i] | PooledAsscData_df$Trait_B == Traits[i]),] 
    tempData = tempData[which(tempData$Trait_A == Traits[j] | tempData$Trait_B == Traits[j]), ]
    if(length(which(tempData$Percent == 1)) == 1){
      AllAsscData[k,] = tempData[which(tempData$Percent == 1),]
      k = k + 1
    }else if(length(which(tempData$Percent == -1) == 1)){
      AllAsscData[k,] = tempData[which(tempData$Percent == -1),]
      k = k + 1
    }else if(length(which(tempData$Difference > 0)) == 1){
      AllAsscData[k, ] = tempData[which(tempData$Difference > 0),]
      k = k + 1
    }else{
      AllAsscData[k, ] = tempData[1,]
      k = k + 1
    }
  }
}

p = ggplot(AllAsscData, aes(Trait_A, Trait_B))
p = p + geom_tile(aes(fill = Percent, alpha = Difference), colour = "black")+scale_fill_gradient(high = "pink", low = "blue")
p = p + theme_bw()+theme(axis.title.x = element_text(size = 8), axis.text.x = element_text(size = 8, angle = 90),axis.text.y = element_text(size = 8))

MeanDiff = mean(AsscData_df$Difference, na.rm = T)

a = ggplot(AsscData_df, aes(x = ID, y = Difference))+geom_point(color = "grey70", pch = 16 ,fill = "black")
a = a + geom_hline(yintercept = MeanDiff, color = "red", lty = 2)
a = a + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"),axis.text.x = element_text(angle = 90, size = 0, color = "white"))  
a = a + xlab("Traits") + ylab("Observed - Expected")
ggsave(paste("AsscData_", Sys.Date(), ".pdf"))




a = ggplot(PooledData_df, aes(x = ID, y = Difference, color = factor(Association)))+geom_point(pch = 19)
a = a + geom_hline(yintercept = NegMeanDiff, color = "red", linetype  = 2)
a = a + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"),axis.text.x = element_text(angle = 90))  
a = a + xlab("Traits") + ylab("Observed - Expected")
ggsave(paste("PooledAssoc_", Sys.Date(), ".pdf"))

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

# Determine the number of observed (+/+)'s for Isolation and Trait for the observed data and compare it to the permuted data # 
AssocIsoTrait_PP_Pos_df = data.frame(Isolation = character(),
                             Trait = character(),
                             Observed_PP = numeric(),
							 Expected_PP = numeric(),
                             Percent_PP = numeric())
k = 1
for(i in 1:length(Isolations)){
  for(j in 1:length(Traits)){
    tempObs = IsolationXTraitCounts_df[which(IsolationXTraitCounts_df[,1] == Isolations[[i]] & IsolationXTraitCounts_df[,2] == Traits[[j]]),3]
    tempList = list()
    for(perm in 1:1000){
      tempDF1 = PermutedAllAssoc_list[[perm]]$Both_P
      COLI= which(colnames(tempDF1) == Isolations[[i]])
      ROWJ= which(rownames(tempDF1) == Traits[[j]])
      COLJ= which(colnames(tempDF1) == Traits[[j]])
      ROWI= which(rownames(tempDF1) == Isolations[[i]])
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
	AssocIsoTrait_PP_Pos_df[k,4] = mean(temporary, na.rm = T)
    AssocIsoTrait_PP_Pos_df[k,5] = ecdf(temporary)(tempObs)
    k = k + 1
  }
}

write.csv(AssocIsoTrait_PP_Pos_df, file = paste("AssocIsoTrait_PP_Pos_df_", Sys.Date(), ".csv", sep = ""), row.names = F)

# Determine the cut-off for signficant postive (+/+) traits in an environment by calculating the values of the data found in the 97.5th% of the distribution #
PosAssocIso_PP_split = splitBy("Isolation", AssocIsoTrait_PP_Pos_df)

PosIsoCutOffs_PP_df = data.frame(Isolation = character(), CutOff_97.5 = numeric(), CutOff_95 = numeric())

for(i in 1:length(PosAssocIso_PP_split)){
	PosIsoCutOffs_PP_df[i,1] = names(PosAssocIso_PP_split)[i]
	PosIsoCutOffs_PP_df[i,2] = quantile(PosAssocIso_PP_split[[i]][,5], 0.975)
	PosIsoCutOffs_PP_df[i,3] = quantile(PosAssocIso_PP_split[[i]][,5], 0.95)
}

write.csv(PosIsoCutOffs_PP_df, file = paste("PosIsoCutOffs_PP_df_", Sys.Date(), ".csv", sep = ""), row.names = F)

# Determine which traits are positive (+,+) significant in an environment #
Pos_SigTraitsIso_PP_list = list()
for(i in 1:nrow(PosIsoCutOffs_PP_df)){
  Pos_SigTraitsIso_PP_list[[i]] = list()
  temp_df = AssocIsoTrait_PP_Pos_df[which(AssocIsoTrait_PP_Pos_df[,1] == PosIsoCutOffs_PP_df[i,1]),]
  Pos_SigTraitsIso_PP_list[[i]]$Isolation = PosIsoCutOffs_PP_df[i,1]
  Pos_SigTraitsIso_PP_list[[i]]$Traits = temp_df[which(temp_df[,5] >= PosIsoCutOffs_PP_df[i,2]),2]
}

# Determine the number of observed (-/-)'s for Isolation and Trait for the observed data and compare it to the permuted data # 
AssocIsoTrait_NN_Pos_df = data.frame(Isolation = character(),
                             Trait = character(),
                             Observed_NN = numeric(),
							 Expected_NN = numeric(),
                             Percent_NN = numeric())
k = 1
for(i in 1:length(Isolations)){
  for(j in 1:length(Traits)){
    tempObs = IsolationXTraitCounts_df[which(IsolationXTraitCounts_df[,1] == Isolations[[i]] & IsolationXTraitCounts_df[,2] == Traits[[j]]),4]
    tempList = list()
    for(perm in 1:1000){
      tempDF1 = PermutedAllAssoc_list[[perm]]$Both_N
      COLI= which(colnames(tempDF1) == Isolations[[i]])
      ROWJ= which(rownames(tempDF1) == Traits[[j]])
      COLJ= which(colnames(tempDF1) == Traits[[j]])
      ROWI= which(rownames(tempDF1) == Isolations[[i]])
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
	AssocIsoTrait_NN_Pos_df[k,4] = mean(temporary, na.rm = T)
    AssocIsoTrait_NN_Pos_df[k,5] = ecdf(temporary)(tempObs)
    k = k + 1
  }
}

write.csv(AssocIsoTrait_NN_Pos_df, file = paste("AssocIsoTrait_NN_Pos_df_", Sys.Date(), ".csv", sep = ""), row.names = F)

# Determine the cut-off for signficant traits in an environment by calculating the values of the data found in the 97.5th% of the distribution #
PosAssocIso_NN_split = splitBy("Isolation", AssocIsoTrait_NN_Pos_df)

PosIsoCutOffs_NN_df = data.frame(Isolation = character(), CutOff_97.5 = numeric(), CutOff_95 = numeric())

for(i in 1:length(PosAssocIso_NN_split)){
	PosIsoCutOffs_NN_df[i,1] = names(PosAssocIso_NN_split)[i]
	PosIsoCutOffs_NN_df[i,2] = quantile(PosAssocIso_NN_split[[i]][,5], 0.975)
	PosIsoCutOffs_NN_df[i,3] = quantile(PosAssocIso_NN_split[[i]][,5], 0.95)
}

write.csv(PosIsoCutOffs_NN_df, file = paste("PosIsoCutOffs_NN_df_", Sys.Date(), ".csv", sep = ""), row.names = F)

# Determine which traits are significant in an environment #
Pos_SigTraitsIso_NN_list = list()
for(i in 1:nrow(PosIsoCutOffs_NN_df)){
  Pos_SigTraitsIso_NN_list[[i]] = list()
  temp_df = AssocIsoTrait_NN_Pos_df[which(AssocIsoTrait_NN_Pos_df[,1] == PosIsoCutOffs_NN_df[i,1]),]
  Pos_SigTraitsIso_NN_list[[i]]$Isolation = PosIsoCutOffs_NN_df[i,1]
  Pos_SigTraitsIso_NN_list[[i]]$Traits = temp_df[which(temp_df[,5] >= PosIsoCutOffs_NN_df[i,2]),2]
}

###### Determine negative associations between Isolation environments and Traits ###### 
# Determine the number of observed (+/-)'s for Isolation and Trait for the observed data and compare it to the permuted data # 
AssocIsoTrait_PN_Neg_df = data.frame(Isolation = character(),
                             Trait = character(),
                             Observed_PN = numeric(),
							 Expected_PN = numeric(),
                             Percent_PN = numeric())
k = 1
for(i in 1:length(Isolations)){
  for(j in 1:length(Traits)){
    tempObs = IsolationXTraitCounts_df[which(IsolationXTraitCounts_df[,1] == Isolations[[i]] & IsolationXTraitCounts_df[,2] == Traits[[j]]),5]
    tempList = list()
    for(perm in 1:1000){
      tempDF1 = PermutedAllAssoc_list[[perm]]$Diff_PN
      COLI= which(colnames(tempDF1) == Isolations[[i]])
      ROWJ= which(rownames(tempDF1) == Traits[[j]])
      COLJ= which(colnames(tempDF1) == Traits[[j]])
      ROWI= which(rownames(tempDF1) == Isolations[[i]])
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
	AssocIsoTrait_PN_Neg_df[k,4] = mean(temporary, na.rm = T)
    AssocIsoTrait_PN_Neg_df[k,5] = ecdf(temporary)(tempObs)
    k = k + 1
  }
}

write.csv(AssocIsoTrait_PN_Neg_df, file = paste("AssocIsoTrait_PN_Neg_df_", Sys.Date(), ".csv", sep = ""), row.names = F)

# Determine the cut-off for signficant traits in an environment by calculating the values of the data found in the 97.5th% of the distribution #
NegAssocIso_PN_split = splitBy("Isolation", AssocIsoTrait_PN_Neg_df)

NegIsoCutOffs_PN_df = data.frame(Isolation = character(), CutOff_97.5 = numeric(), CutOff_95 = numeric())

for(i in 1:length(NegAssocIso_PN_split)){
	NegIsoCutOffs_PN_df[i,1] = names(NegAssocIso_PN_split)[i]
	NegIsoCutOffs_PN_df[i,2] = quantile(NegAssocIso_PN_split[[i]][,5], 0.975)
	NegIsoCutOffs_PN_df[i,3] = quantile(NegAssocIso_PN_split[[i]][,5], 0.95)
}

write.csv(NegIsoCutOffs_PN_df, file = paste("NegIsoCutOffs_PN_df_", Sys.Date(), ".csv", sep = ""), row.names = F)

# Determine which traits are significant in an environment #
Neg_SigTraitsIso_PN_list = list()
for(i in 1:nrow(NegIsoCutOffs_PN_df)){
  Neg_SigTraitsIso_PN_list[[i]] = list()
  temp_df = AssocIsoTrait_PN_Neg_df[which(AssocIsoTrait_PN_Neg_df[,1] == NegIsoCutOffs_PN_df[i,1]),]
  Neg_SigTraitsIso_PN_list[[i]]$Isolation = NegIsoCutOffs_PN_df[i,1]
  Neg_SigTraitsIso_PN_list[[i]]$Traits = temp_df[which(temp_df[,5] >= NegIsoCutOffs_PN_df[i,2]),2]
}


# Determine the number of observed (-/+)'s for Isolation and Trait for the observed data and compare it to the permuted data # 
AssocIsoTrait_NP_Neg_df = data.frame(Isolation = character(),
                             Trait = character(),
                             Observed_NP = numeric(),
							 Expected_NP = numeric(),
                             Percent_NP = numeric())
k = 1
for(i in 1:length(Isolations)){
  for(j in 1:length(Traits)){
    tempObs = IsolationXTraitCounts_df[which(IsolationXTraitCounts_df[,1] == Isolations[[i]] & IsolationXTraitCounts_df[,2] == Traits[[j]]),6]
    tempList = list()
    for(perm in 1:1000){
      tempDF1 = PermutedAllAssoc_list[[perm]]$Diff_NP
      COLI= which(colnames(tempDF1) == Isolations[[i]])
      ROWJ= which(rownames(tempDF1) == Traits[[j]])
      COLJ= which(colnames(tempDF1) == Traits[[j]])
      ROWI= which(rownames(tempDF1) == Isolations[[i]])
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
	AssocIsoTrait_NP_Neg_df[k,4] = mean(temporary, na.rm = T)
    AssocIsoTrait_NP_Neg_df[k,5] = ecdf(temporary)(tempObs)
    k = k + 1
  }
}

write.csv(AssocIsoTrait_NP_Neg_df, file = paste("AssocIsoTrait_NP_Neg_df_", Sys.Date(), ".csv", sep = ""), row.names = F)

# Determine the cut-off for signficant traits in an environment by calculating the values of the data found in the 97.5th% of the distribution #
NegAssocIso_NP_split = splitBy("Isolation", AssocIsoTrait_NP_Neg_df)

NegIsoCutOffs_NP_df = data.frame(Isolation = character(), CutOff_97.5 = numeric(), CutOff_95 = numeric())

for(i in 1:length(NegAssocIso_NP_split)){
	NegIsoCutOffs_NP_df[i,1] = names(NegAssocIso_NP_split)[i]
	NegIsoCutOffs_NP_df[i,2] = quantile(NegAssocIso_NP_split[[i]][,5], 0.975)
	NegIsoCutOffs_NP_df[i,3] = quantile(NegAssocIso_NP_split[[i]][,5], 0.95)
}

write.csv(NegIsoCutOffs_NP_df, file = paste("NegIsoCutOffs_NP_df_", Sys.Date(), ".csv", sep = ""), row.names = F)

# Determine which traits are significant in an environment #
Neg_SigTraitsIso_NP_list = list()
for(i in 1:nrow(NegIsoCutOffs_NP_df)){
  Neg_SigTraitsIso_NP_list[[i]] = list()
  temp_df = AssocIsoTrait_NP_Neg_df[which(AssocIsoTrait_NP_Neg_df[,1] == NegIsoCutOffs_NP_df[i,1]),]
  Neg_SigTraitsIso_NP_list[[i]]$Isolation = NegIsoCutOffs_NP_df[i,1]
  Neg_SigTraitsIso_NP_list[[i]]$Traits = temp_df[which(temp_df[,5] >= NegIsoCutOffs_NP_df[i,2]),2]
}

#############################################################################################################
########### Calculates counts of each type of trait presence/absence [(+/+), (-/-), (+/-), (-/+)]  within an isolation environment  #########
########### for all trait pairs where PP is (+/+), NN is (-/-), PN is (+/-) and NP is (-/+)                                                                                  ###########  
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

###### Isolation Trait pair analysis for Positively Associated traits (+,+) ######
# Determine the number of observed (+/+)'s for Trait pairs in an isolation environment for the observed data and compare it to the permuted data # 
IsoAnalysis_PP_df = data.frame(Isolation = character(),
                                  Trait_A = character(),
                                  Trait_B = character(),
                                  Observed_PP = numeric(), 
								  Expected_PP = numeric(),
                                  Percent_PP = numeric())

k = 1
for(iso in 1:length(Isolations)){
  tempObs = IsoAllAssoc_list[[iso]]$Both_P
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
	  IsoAnalysis_PP_df[k,5] = mean(templist, na.rm = T)
      IsoAnalysis_PP_df[k,6] = ecdf(templist)(ObsValue)
      k = k+1
    }
  } 
}

write.csv(IsoAnalysis_PP_df, file = paste("IsoAnalysis_PP_df_", Sys.Date(), ".csv", sep = ""), row.names = F)


# Summary of permuted data trait pairs within an isolation environment - Postive (+,+) #
PermutedIsoSum_PP_df = data.frame(Isolation = character(),
                                  Trait_A = character(),
                                  Trait_B = character(),
                                  Mean_PP = numeric(), 
                                  SD_PP = numeric(), 
                                  Var_PP = numeric())

k = 1
for(iso in 1:length(Isolations)){
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
      PermutedIsoSum_PP_df[k,1] = Isolations[iso]
      PermutedIsoSum_PP_df[k,2] = Traits[i]
      PermutedIsoSum_PP_df[k,3] = Traits[j]
      PermutedIsoSum_PP_df[k,4] = mean(templist, na.rm = T)
      PermutedIsoSum_PP_df[k,5] =sd(templist, na.rm = T)
      PermutedIsoSum_PP_df[k,6] =var(templist, na.rm = T)
      k = k+1
    }
  } 
}

write.csv(PermutedIsoSum_PP_df, file = paste("PermutedIsoSum_PP_df_", Sys.Date(), ".csv", sep = ""), row.names = F)

# Limit Trait pairs in the isolation environment to trait pairs where one trait is significantly associated with the isolation environment - Positive (+/+) #
AllSigIsoAnalysis_PP_df = matrix(0, nrow = 1 ,ncol = ncol(IsoAnalysis_PP_df))
colnames(AllSigIsoAnalysis_PP_df) = colnames(IsoAnalysis_PP_df)

# Isolation Trait Pair Visualization - Positive (+,+) #
AllSigIsoAnalysis_PP_df = data.frame(AllSigIsoAnalysis_PP_df)

for(i in 1:length(Pos_SigTraitsIso_PP_list)){
  iso = Pos_SigTraitsIso_PP_list[[i]]$Isolation
  keepTraits = Pos_SigTraitsIso_PP_list[[i]]$Traits
  tempData = IsoAnalysis_PP_df[which(IsoAnalysis_PP_df$Isolation == iso),]
  dataKeep = tempData[which(tempData$Trait_A %in% keepTraits | tempData$Trait_B %in% keepTraits),]
  AllSigIsoAnalysis_PP_df = rbind(AllSigIsoAnalysis_PP_df, dataKeep)
}

write.csv(AllSigIsoAnalysis_PP_df, file = paste("AllSigIsoAnalysis_PP_df_", Sys.Date(), ".csv", sep = ""), row.names = F)

####### Isolation Trait Pair Visualization - Positive (+,+) ####### 
AllIsoSum_PP_df = merge(IsoAnalysis_PP_df, PermutedIsoSum_PP_df)
AllIsoSum_PP_df$Difference = AllIsoSum_PP_df$Observed_PP - AllIsoSum_PP_df$Mean_PP
AllIsoSum_PP_df$Color = rep(0, nrow(AllIsoSum_PP_df))

for(i in 1:nrow(AllIsoSum_PP_df)){
  SigValue = PosIsoCutOffs_PP_df[which(PosIsoCutOffs_PP_df$Isolation == AllIsoSum_PP_df[i,"Isolation"]), 2] 
  if(AllIsoSum_PP_df[i,"Difference"] > 0 & AllIsoSum_PP_df[i,"Percent_PP"] >= SigValue){
    AllIsoSum_PP_df[i, "Color"] = 1
  }else if(AllIsoSum_PP_df[i,"Difference"] < 0 & AllIsoSum_PP_df[i,"Percent_PP"] >= SigValue){
    AllIsoSum_PP_df[i, "Color"] = 2
  }else{
    AllIsoSum_PP_df[i, "Color"] = 3
  }
}
write.csv(AllIsoSum_PP_df, file = paste("AllIsoSum_PP_df_", Sys.Date(), ".csv", sep = ""), row.names = F)


pdf(paste("Isolation_TraitPair_PP_", Sys.Date(), ".pdf", sep = ""))
for(iso in 1:length(Isolations)){
  tempDF = AllIsoSum_PP_df[which(AllIsoSum_PP_df$Isolation == Isolations[iso]),]
  lengthCheck = length(unique(tempDF$Color))
  Values = unique(tempDF$Color)
  if(lengthCheck == 3){
    colorValues = c("#5C0023", "#034A58", "grey80")
  }else if(lengthCheck == 2 & length(which(Values == 1) > 0)){
    colorValues = c("#5C0023", "grey80")
  }else if(lengthCheck == 2 & length(which(Values == 2) > 0)){
    colorValues = c("#034A58", "grey80")
  }
  Mean_PP_Diff = mean(tempDF$Difference)
  Mean_PP_Diff = round(Mean_PP_Diff)
  
  a = ggplot(tempDF, aes(x = Trait_A, y = Difference))+geom_point()+scale_color_manual(values = colorValues)+geom_hline(yintercept = Mean_PP_Diff, color = "red")
  a = a + theme_bw()+ ggtitle(Isolations[iso]) + xlab("Trait Pairs")+ylab("Count")+scale_fill_grey()+ theme(axis.title.x = element_text(size = 12), axis.text.x = element_text(size = 0, angle = 90),axis.text.y = element_text(size = 8))
  print(a)
}
dev.off()

IsoPairStat_PP_df = data.frame(Isolation = character(), AvgDiff = numeric(), VarDiff = numeric(), AvgTraitCounts = numeric(), VarTraitCounts = numeric())
for(iso in 1:length(Isolations)){
  tempDF = AllIsoSum_PP_df[which(AllIsoSum_PP_df$Isolation == Isolations[iso]),]
  IsoPairStat_PP_df[iso, 1] = Isolations[iso]
  IsoPairStat_PP_df[iso, 2] = mean(tempDF$Difference, na.rm = TRUE)
  IsoPairStat_PP_df[iso, 3] = var(tempDF$Difference, na.rm = TRUE)
  IsoPairStat_PP_df[iso, 4] = mean(tempDF$Observed_PP, na.rm = TRUE)
  IsoPairStat_PP_df[iso, 5] = var(tempDF$Observed_PP, na.rm = TRUE)
}
write.csv(IsoPairStat_PP_df, file = paste("IsoPairStat_PP_df_", Sys.Date(), ".csv", sep = ""), row.names = F)

a = ggplot(IsoPairStat_PP_df, aes(x = Isolation, y = AvgDiff))+geom_point(pch = 19, fill = "grey60")
a = a + geom_errorbar(aes(ymin=AvgDiff-VarDiff, ymax=AvgDiff+VarDiff))
a = a + geom_hline(yintercept = 0, color = "red", linetype  = 2)
a = a + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"),axis.text.x = element_text(angle = 90))  
a = a + xlab("Isolation") + ylab("Average Difference")
ggsave(paste("IsoPairStat_PP_", Sys.Date(), ".pdf", sep = ""))

a = ggplot(IsoPairStat_PP_df, aes(x = Isolation, y = AvgTraitCounts))+geom_point(pch = 19, fill = "grey60")
a = a + geom_errorbar(aes(ymin=AvgTraitCounts-VarTraitCounts, ymax=AvgTraitCounts+VarTraitCounts))
a = a + geom_hline(yintercept = 0, color = "red", linetype  = 2)
a = a + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"),axis.text.x = element_text(angle = 90))  
a = a + xlab("Isolation") + ylab("Average Trait Counts")
ggsave(paste("IsoPair_TraitCount_PP_", Sys.Date(), ".pdf", sep = ""))

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

# Limit Trait pairs in the isolation environment to trait pairs where both traits are significantly associated with the isolation environment - Negative (-/+) #
Both_IsoSigSub_PP_df = matrix(0, nrow = 1 ,ncol = ncol(AllIsoSum_PP_df))
colnames(Both_IsoSigSub_PP_df) = colnames(AllIsoSum_PP_df)

Both_IsoSigSub_PP_df = data.frame(Both_IsoSigSub_PP_df)

for(i in 1:length(Pos_SigTraitsIso_PP_list)){
  iso = Pos_SigTraitsIso_PP_list[[i]]$Isolation
  keepTraits = Pos_SigTraitsIso_PP_list[[i]]$Traits
  tempData = AllIsoSum_PP_df[which(AllIsoSum_PP_df$Isolation == iso),]
  dataKeep = tempData[which(tempData$Trait_A %in% keepTraits),] 
  dataKeep = dataKeep[which(dataKeep$Trait_B %in% keepTraits),]
  Both_IsoSigSub_PP_df = rbind(Both_IsoSigSub_PP_df, dataKeep)
}

Both_IsoSigSub_PP_df = Both_IsoSigSub_PP_df[-1,]

write.csv(Both_IsoSigSub_PP_df, file = paste("Both_IsoSigSub_PP_df_", Sys.Date(), ".csv", sep = ""), row.names = F)

pdf(paste("BothSig_Isolation_TraitPair_PP_", Sys.Date(), ".pdf", sep = ""))
for(iso in 1:length(Isolations)){
  tempDF = Both_IsoSigSub_PP_df[which(Both_IsoSigSub_PP_df$Isolation == Isolations[iso]),]
  lengthCheck = length(unique(tempDF$Color))
  Values = unique(tempDF$Color)
  if(lengthCheck == 3){
    colorValues = c("#5C0023", "#034A58", "grey80")
  }else if(lengthCheck == 2 & length(which(Values == 1) > 0)){
    colorValues = c("#5C0023", "grey80")
  }else if(lengthCheck == 2 & length(which(Values == 2) > 0)){
    colorValues = c("#034A58", "grey80")
  }
  Both_Mean_PP_Diff = mean(tempDF$Difference)
  a = ggplot(tempDF, aes(x = Trait_A, y = Difference))+geom_point()+scale_color_manual(values = colorValues)+geom_hline(yintercept = Both_Mean_PP_Diff, color = "red")
  a = a + theme_bw()+ ggtitle(Isolations[iso]) + xlab("Trait Pairs")+ylab("Count")+scale_fill_grey()+ theme(axis.title.x = element_text(size = 12), axis.text.x = element_text(size = 0, angle = 90),axis.text.y = element_text(size = 8))
  print(a)
}
dev.off()

Both_IsoStat_PP_df = data.frame(Isolation = character(), AvgDiff = numeric(), VarDiff = numeric(), AvgTraitCounts = numeric(), VarTraitCounts = numeric())
for(iso in 1:length(Isolations)){
  tempDF = Both_IsoSigSub_PP_df[which(Both_IsoSigSub_PP_df$Isolation == Isolations[iso]),]
  Both_IsoStat_PP_df[iso, 1] = Isolations[iso]
  Both_IsoStat_PP_df[iso, 2] = mean(tempDF$Difference, na.rm = TRUE)
  Both_IsoStat_PP_df[iso, 3] = var(tempDF$Difference, na.rm = TRUE)
  Both_IsoStat_PP_df[iso, 4] = mean(tempDF$Observed_PP, na.rm = TRUE)
  Both_IsoStat_PP_df[iso, 5] = var(tempDF$Observed_PP, na.rm = TRUE)
}

write.csv(Both_IsoStat_PP_df, file = paste("Both_IsoStat_PP_df_", Sys.Date(), ".csv", sep = ""), row.names = F)

a = ggplot(Both_IsoStat_PP_df, aes(x = Isolation, y = AvgDiff))+geom_point(pch = 19, fill = "grey60")
a = a + geom_errorbar(aes(ymin=AvgDiff-VarDiff, ymax=AvgDiff+VarDiff))
a = a + geom_hline(yintercept = 0, color = "red", linetype  = 2)
a = a + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"),axis.text.x = element_text(angle = 90))  
a = a + xlab("Isolation") + ylab("Average Difference")
ggsave(file = paste("Both_IsoStat_PP_", Sys.Date(), ".pdf", sep = ""))

a = ggplot(Both_IsoStat_PP_df, aes(x = Isolation, y = AvgTraitCounts))+geom_point(pch = 19, fill = "grey60")
a = a + geom_errorbar(aes(ymin=AvgTraitCounts-VarTraitCounts, ymax=AvgTraitCounts+VarTraitCounts))
a = a + geom_hline(yintercept = 0, color = "red", linetype  = 2)
a = a + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"),axis.text.x = element_text(angle = 90))  
a = a + xlab("Isolation") + ylab("Average Trait Counts")
ggsave(paste("Both_TraitCount_PP_", Sys.Date(), ".pdf", sep = ""))

All_IsoSigSub_PP_df = matrix(0, nrow = 1 ,ncol = ncol(AllIsoSum_PP_df))
colnames(All_IsoSigSub_PP_df) = colnames(AllIsoSum_PP_df)

All_IsoSigSub_PP_df = data.frame(All_IsoSigSub_PP_df)

for(i in 1:length(Pos_SigTraitsIso_PP_list)){
  iso = Pos_SigTraitsIso_PP_list[[i]]$Isolation
  keepTraits = Pos_SigTraitsIso_PP_list[[i]]$Traits
  tempData = AllIsoSum_PP_df[which(AllIsoSum_PP_df$Isolation == iso),]
  dataKeep = tempData[which(tempData$Trait_A %in% keepTraits | tempData$Trait_B %in% keepTraits),]
  All_IsoSigSub_PP_df = rbind(All_IsoSigSub_PP_df, dataKeep)
}

All_IsoSigSub_PP_df = All_IsoSigSub_PP_df[-1,]

write.csv(All_IsoSigSub_PP_df, file = paste("All_IsoSigSub_PP_df", Sys.Date(), ".csv", sep = ""), row.names = F)

pdf(paste("AllSig_Isolation_TraitPair_PP_", Sys.Date(), ".pdf", sep = ""))
for(iso in 1:length(Isolations)){
  tempDF = All_IsoSigSub_PP_df[which(All_IsoSigSub_PP_df$Isolation == Isolations[iso]),]
  lengthCheck = length(unique(tempDF$Color))
  Values = unique(tempDF$Color)
  if(lengthCheck == 3){
    colorValues = c("#5C0023", "#034A58", "grey80")
  }else if(lengthCheck == 2 & length(which(Values == 1) > 0)){
    colorValues = c("#5C0023", "grey80")
  }else if(lengthCheck == 2 & length(which(Values == 2) > 0)){
    colorValues = c("#034A58", "grey80")
  }
  All_Mean_PP_Diff = mean(tempDF$Difference)
  a = ggplot(tempDF, aes(x = Trait_A, y = Difference))+geom_point()+scale_color_manual(values = colorValues)+geom_hline(yintercept = All_Mean_PP_Diff, color = "red")
  a = a + theme_bw()+ ggtitle(Isolations[iso]) + xlab("Trait Pairs")+ylab("Count")+scale_fill_grey()+ theme(axis.title.x = element_text(size = 12), axis.text.x = element_text(size = 0, angle = 90),axis.text.y = element_text(size = 8))
  print(a)
}
dev.off()

All_IsoStat_PP_df = data.frame(Isolation = character(), AvgDiff = numeric(), VarDiff = numeric(), AvgTraitCounts = numeric(), VarTraitCounts = numeric())
for(iso in 1:length(Isolations)){
  tempDF = All_IsoSigSub_PP_df[which(All_IsoSigSub_PP_df$Isolation == Isolations[iso]),]
  All_IsoStat_PP_df[iso, 1] = Isolations[iso]
  All_IsoStat_PP_df[iso, 2] = mean(tempDF$Difference, na.rm = TRUE)
  All_IsoStat_PP_df[iso, 3] = var(tempDF$Difference, na.rm = TRUE)
  All_IsoStat_PP_df[iso, 4] = mean(tempDF$Observed_PP, na.rm = TRUE)
  All_IsoStat_PP_df[iso, 5] = var(tempDF$Observed_PP, na.rm = TRUE)
}

write.csv(All_IsoStat_PP_df, file = paste("All_IsoStat_PP_df_", Sys.Date(), ".csv", sep = ""), row.names = F)

a = ggplot(All_IsoStat_PP_df, aes(x = Isolation, y = AvgDiff))+geom_point(pch = 19, fill = "grey60")
a = a + geom_errorbar(aes(ymin=AvgDiff-VarDiff, ymax=AvgDiff+VarDiff))
a = a + geom_hline(yintercept = 0, color = "red", linetype  = 2)
a = a + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"),axis.text.x = element_text(angle = 90))  
a = a + xlab("Isolation") + ylab("Average Difference")
ggsave(file = paste("All_IsoStat_PP_", Sys.Date(), ".pdf", sep = ""))

a = ggplot(All_IsoStat_PP_df, aes(x = Isolation, y = AvgTraitCounts))+geom_point(pch = 19, fill = "grey60")
a = a + geom_errorbar(aes(ymin=AvgTraitCounts-VarTraitCounts, ymax=AvgTraitCounts+VarTraitCounts))
a = a + geom_hline(yintercept = 0, color = "red", linetype  = 2)
a = a + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"),axis.text.x = element_text(angle = 90))  
a = a + xlab("Isolation") + ylab("Average Trait Counts")
ggsave(paste("All_TraitCount_PP_", Sys.Date(), ".pdf", sep = ""))

###### Isolation Trait pair analysis for Positively Associated traits (-,-) ######
# Determine the number of observed (-/-)'s for Trait pairs in an isolation environment for the observed data and compare it to the permuted data # 
IsoAnalysis_NN_df = data.frame(Isolation = character(),
                               Trait_A = character(),
                               Trait_B = character(),
                               Observed_NN = numeric(), 
                               Percent_NN = numeric())

k = 1
for(iso in 1:length(Isolations)){
  tempObs = IsoAllAssoc_list[[iso]]$Both_N
  for(i in 1:(length(Traits)-1)){
    for(j in (i+1):length(Traits)){
      temp_list = list()
      for(perm in 1:1000){
        Val1 = PermutedIsoAllAssoc_list[[perm]][[iso]]$Both_N[i,j]
        Val2 = PermutedIsoAllAssoc_list[[perm]][[iso]]$Both_N[j,i]
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

# Summary of permuted data trait pairs within an isolation environment - Postive (-,-) #
PermutedIsoSum_NN_df = data.frame(Isolation = character(),
                                  Trait_A = character(),
                                  Trait_B = character(),
                                  Mean_NN = numeric(), 
                                  SD_NN = numeric(), 
                                  Var_NN = numeric())

k = 1
for(iso in 1:length(Isolations)){
  for(i in 1:(length(Traits)-1)){
    for(j in (i+1):length(Traits)){
      temp_list = list()
      for(perm in 1:1000){
        Val1 = PermutedIsoAllAssoc_list[[perm]][[iso]]$Both_N[i,j]
        Val2 = PermutedIsoAllAssoc_list[[perm]][[iso]]$Both_N[j,i]
        if(Val1 != Val2 & Val1 > Val2){
          temp_list[perm] = Val1
        }else{ 
          temp_list[perm] = Val2
        }
      }
      templist = unlist(temp_list)
      PermutedIsoSum_NN_df[k,1] = Isolations[iso]
      PermutedIsoSum_NN_df[k,2] = Traits[i]
      PermutedIsoSum_NN_df[k,3] = Traits[j]
      PermutedIsoSum_NN_df[k,4] = mean(templist, na.rm = T)
      PermutedIsoSum_NN_df[k,5] =sd(templist, na.rm = T)
      PermutedIsoSum_NN_df[k,6] =var(templist, na.rm = T)
      k = k+1
    }
  } 
}

write.csv(PermutedIsoSum_NN_df, file = paste("PermutedIsoSum_NN_df_", Sys.Date(), ".csv", sep = ""), row.names = F)

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

####### Isolation Trait Pair Visualization - Positive (+,+) ####### 
AllIsoSum_NN_df = merge(IsoAnalysis_NN_df, PermutedIsoSum_NN_df)
AllIsoSum_NN_df$Difference = AllIsoSum_NN_df$Observed_NN - AllIsoSum_NN_df$Mean_NN
AllIsoSum_NN_df$Color = rep(0, nrow(AllIsoSum_NN_df))

for(i in 1:nrow(AllIsoSum_NN_df)){
  SigValue = PosIsoCutOffs_NN_df[which(PosIsoCutOffs_NN_df$Isolation == AllIsoSum_NN_df[i,"Isolation"]), 2] 
  if(AllIsoSum_NN_df[i,"Difference"] > 0 & AllIsoSum_NN_df[i,"Percent_NN"] >= SigValue){
    AllIsoSum_NN_df[i, "Color"] = 1
  }else if(AllIsoSum_NN_df[i,"Difference"] < 0 & AllIsoSum_NN_df[i,"Percent_NN"] >= SigValue){
    AllIsoSum_NN_df[i, "Color"] = 2
  }else{
    AllIsoSum_NN_df[i, "Color"] = 3
  }
}

write.csv(AllIsoSum_NN_df, file = paste("AllIsoSum_NN_df_", Sys.Date(), ".csv", sep = ""), row.names = F)


pdf(paste("Isolation_TraitPair_NN_", Sys.Date(), ".pdf", sep = ""))
for(iso in 1:length(Isolations)){
  tempDF = AllIsoSum_NN_df[which(AllIsoSum_NN_df$Isolation == Isolations[iso]),]
  lengthCheck = length(unique(tempDF$Color))
  Values = unique(tempDF$Color)
  if(lengthCheck == 3){
    colorValues = c("#5C0023", "#034A58", "grey80")
  }else if(lengthCheck == 2 & length(which(Values == 1) > 0)){
    colorValues = c("#5C0023", "grey80")
  }else if(lengthCheck == 2 & length(which(Values == 2) > 0)){
    colorValues = c("#034A58", "grey80")
  }
  
  Mean_NN_Diff = mean(tempDF$Difference)
  
  a = ggplot(tempDF, aes(x = Trait_A, y = Difference))+geom_point()+scale_color_manual(values = colorValues)+geom_hline(yintercept = Mean_NN_Diff, color = "red")
  a = a + theme_bw()+ ggtitle(Isolations[iso]) + xlab("Trait Pairs")+ylab("Count")+scale_fill_grey()+ theme(axis.title.x = element_text(size = 12), axis.text.x = element_text(size = 0, angle = 90),axis.text.y = element_text(size = 8))
  print(a)
}
dev.off()


IsoPairStat_NN_df = data.frame(Isolation = character(), AvgDiff = numeric(), VarDiff = numeric(), AvgTraitCounts = numeric(), VarTraitCounts = numeric())
for(iso in 1:length(Isolations)){
  tempDF = AllIsoSum_NN_df[which(AllIsoSum_NN_df$Isolation == Isolations[iso]),]
  IsoPairStat_NN_df[iso, 1] = Isolations[iso]
  IsoPairStat_NN_df[iso, 2] = mean(tempDF$Difference, na.rm = TRUE)
  IsoPairStat_NN_df[iso, 3] = var(tempDF$Difference, na.rm = TRUE)
  IsoPairStat_NN_df[iso, 4] = mean(tempDF$Observed_NN, na.rm = TRUE)
  IsoPairStat_NN_df[iso, 5] = var(tempDF$Observed_NN, na.rm = TRUE)
}

write.csv(IsoPairStat_NN_df, file = paste("IsoPairStat_NN_df_", Sys.Date(), ".csv", sep = ""), row.names = F)

a = ggplot(IsoPairStat_NN_df, aes(x = Isolation, y = AvgDiff))+geom_point(pch = 19, fill = "grey60")
a = a + geom_errorbar(aes(ymin=AvgDiff-VarDiff, ymax=AvgDiff+VarDiff))
a = a + geom_hline(yintercept = 0, color = "red", linetype  = 2)
a = a + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"),axis.text.x = element_text(angle = 90))  
a = a + xlab("Isolation") + ylab("Average Difference")
ggsave(paste("IsoPairStat_NN_", Sys.Date(), ".pdf", sep = ""))

a = ggplot(IsoPairStat_NN_df, aes(x = Isolation, y = AvgTraitCounts))+geom_point(pch = 19, fill = "grey60")
a = a + geom_errorbar(aes(ymin=AvgTraitCounts-VarTraitCounts, ymax=AvgTraitCounts+VarTraitCounts))
a = a + geom_hline(yintercept = 0, color = "red", linetype  = 2)
a = a + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"),axis.text.x = element_text(angle = 90))  
a = a + xlab("Isolation") + ylab("Average Trait Counts")
ggsave(paste("IsoPair_TraitCount_NN_", Sys.Date(), ".pdf", sep = ""))

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

# Limit Trait pairs in the isolation environment to trait pairs where both traits are significantly associated with the isolation environment - Negative (-/+) #
Both_IsoSigSub_NN_df = matrix(0, nrow = 1 ,ncol = ncol(AllIsoSum_NN_df))
colnames(Both_IsoSigSub_NN_df) = colnames(AllIsoSum_NN_df)

Both_IsoSigSub_NN_df = data.frame(Both_IsoSigSub_NN_df)

for(i in 1:length(Pos_SigTraitsIso_NN_list)){
  iso = Pos_SigTraitsIso_NN_list[[i]]$Isolation
  keepTraits = Pos_SigTraitsIso_NN_list[[i]]$Traits
  tempData = AllIsoSum_NN_df[which(AllIsoSum_NN_df$Isolation == iso),]
  dataKeep = tempData[which(tempData$Trait_A %in% keepTraits),] 
  dataKeep = dataKeep[which(dataKeep$Trait_B %in% keepTraits),]
  Both_IsoSigSub_NN_df = rbind(Both_IsoSigSub_NN_df, dataKeep)
}

Both_IsoSigSub_NN_df = Both_IsoSigSub_NN_df[-1,]
write.csv(Both_IsoSigSub_NN_df, file = paste("Both_IsoSigSub_NN_df_", Sys.Date(), ".csv", sep = ""), row.names = F)

pdf(paste("BothSig_Isolation_TraitPair_NN_", Sys.Date(), ".pdf", sep = ""))
for(iso in 1:length(Isolations)){
  tempDF = Both_IsoSigSub_NN_df[which(Both_IsoSigSub_NN_df$Isolation == Isolations[iso]),]
  lengthCheck = length(unique(tempDF$Color))
  Values = unique(tempDF$Color)
  if(lengthCheck == 3){
    colorValues = c("#5C0023", "#034A58", "grey80")
  }else if(lengthCheck == 2 & length(which(Values == 1) > 0)){
    colorValues = c("#5C0023", "grey80")
  }else if(lengthCheck == 2 & length(which(Values == 2) > 0)){
    colorValues = c("#034A58", "grey80")
  }
  Both_Mean_NN_Diff = mean(tempDF$Difference)
  
  a = ggplot(tempDF, aes(x = Trait_A, y = Difference))+geom_point()+scale_color_manual(values = colorValues)+geom_hline(yintercept = Both_Mean_NN_Diff, color = "red")
  a = a + theme_bw()+ ggtitle(Isolations[iso]) + xlab("Trait Pairs")+ylab("Count")+scale_fill_grey()+ theme(axis.title.x = element_text(size = 12), axis.text.x = element_text(size = 0, angle = 90),axis.text.y = element_text(size = 8))
  print(a)
}
dev.off()

Both_IsoStat_NN_df = data.frame(Isolation = character(), AvgDiff = numeric(), VarDiff = numeric(), AvgTraitCounts = numeric(), VarTraitCounts = numeric())
for(iso in 1:length(Isolations)){
  tempDF = Both_IsoSigSub_NN_df[which(Both_IsoSigSub_NN_df$Isolation == Isolations[iso]),]
  Both_IsoStat_NN_df[iso, 1] = Isolations[iso]
  Both_IsoStat_NN_df[iso, 2] = mean(tempDF$Difference, na.rm = TRUE)
  Both_IsoStat_NN_df[iso, 3] = var(tempDF$Difference, na.rm = TRUE)
  Both_IsoStat_NN_df[iso, 4] = mean(tempDF$Observed_NN, na.rm = TRUE)
  Both_IsoStat_NN_df[iso, 5] = var(tempDF$Observed_NN, na.rm = TRUE)
}

write.csv(Both_IsoStat_NN_df, file = paste("Both_IsoStat_NN_df_", Sys.Date(), ".csv", sep = ""), row.names = F)

a = ggplot(Both_IsoStat_NN_df, aes(x = Isolation, y = AvgDiff))+geom_point(pch = 19, fill = "grey60")
a = a + geom_errorbar(aes(ymin=AvgDiff-VarDiff, ymax=AvgDiff+VarDiff))
a = a + geom_hline(yintercept = 0, color = "red", linetype  = 2)
a = a + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"),axis.text.x = element_text(angle = 90))  
a = a + xlab("Isolation") + ylab("Average Difference")
ggsave(file = paste("Both_IsoStat_NN_", Sys.Date(), ".pdf", sep = ""))

a = ggplot(Both_IsoStat_NN_df, aes(x = Isolation, y = AvgTraitCounts))+geom_point(pch = 19, fill = "grey60")
a = a + geom_errorbar(aes(ymin=AvgTraitCounts-VarTraitCounts, ymax=AvgTraitCounts+VarTraitCounts))
a = a + geom_hline(yintercept = 0, color = "red", linetype  = 2)
a = a + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"),axis.text.x = element_text(angle = 90))  
a = a + xlab("Isolation") + ylab("Average Trait Counts")
ggsave(paste("Both_TraitCount_NN_", Sys.Date(), ".pdf", sep = ""))

All_IsoSigSub_NN_df = matrix(0, nrow = 1 ,ncol = ncol(AllIsoSum_NN_df))
colnames(All_IsoSigSub_NN_df) = colnames(AllIsoSum_NN_df)

All_IsoSigSub_NN_df = data.frame(All_IsoSigSub_NN_df)

for(i in 1:length(Pos_SigTraitsIso_NN_list)){
  iso = Pos_SigTraitsIso_NN_list[[i]]$Isolation
  keepTraits = Pos_SigTraitsIso_NN_list[[i]]$Traits
  tempData = AllIsoSum_NN_df[which(AllIsoSum_NN_df$Isolation == iso),]
  dataKeep = tempData[which(tempData$Trait_A %in% keepTraits | tempData$Trait_B %in% keepTraits),]
  All_IsoSigSub_NN_df = rbind(All_IsoSigSub_NN_df, dataKeep)
}

All_IsoSigSub_NN_df = All_IsoSigSub_NN_df[-1,]
write.csv(All_IsoSigSub_NN_df, file = paste("All_IsoSigSub_NN_df", Sys.Date(), ".csv", sep = ""), row.names = F)

pdf(paste("AllSig_Isolation_TraitPair_NN_", Sys.Date(), ".pdf", sep = ""))
for(iso in 1:length(Isolations)){
  tempDF = All_IsoSigSub_NN_df[which(All_IsoSigSub_NN_df$Isolation == Isolations[iso]),]
  lengthCheck = length(unique(tempDF$Color))
  Values = unique(tempDF$Color)
  if(lengthCheck == 3){
    colorValues = c("#5C0023", "#034A58", "grey80")
  }else if(lengthCheck == 2 & length(which(Values == 1) > 0)){
    colorValues = c("#5C0023", "grey80")
  }else if(lengthCheck == 2 & length(which(Values == 2) > 0)){
    colorValues = c("#034A58", "grey80")
  }
  All_Mean_NN_Diff = mean(tempDF$Difference)
  
  a = ggplot(tempDF, aes(x = Trait_A, y = Difference))+geom_point()+scale_color_manual(values = colorValues)+geom_hline(yintercept = All_Mean_NN_Diff, color = "red")
  a = a + theme_bw()+ ggtitle(Isolations[iso]) + xlab("Trait Pairs")+ylab("Count")+scale_fill_grey()+ theme(axis.title.x = element_text(size = 12), axis.text.x = element_text(size = 0, angle = 90),axis.text.y = element_text(size = 8))
  print(a)
}
dev.off()

All_IsoStat_NN_df = data.frame(Isolation = character(), AvgDiff = numeric(), VarDiff = numeric(), AvgTraitCounts = numeric(), VarTraitCounts = numeric())
for(iso in 1:length(Isolations)){
  tempDF = All_IsoSigSub_NN_df[which(All_IsoSigSub_NN_df$Isolation == Isolations[iso]),]
  All_IsoStat_NN_df[iso, 1] = Isolations[iso]
  All_IsoStat_NN_df[iso, 2] = mean(tempDF$Difference, na.rm = TRUE)
  All_IsoStat_NN_df[iso, 3] = var(tempDF$Difference, na.rm = TRUE)
  All_IsoStat_NN_df[iso, 4] = mean(tempDF$Observed_NN, na.rm = TRUE)
  All_IsoStat_NN_df[iso, 5] = var(tempDF$Observed_NN, na.rm = TRUE)
}

write.csv(All_IsoStat_NN_df, file = paste("All_IsoStat_NN_df_", Sys.Date(), ".csv", sep = ""), row.names = F)

a = ggplot(All_IsoStat_NN_df, aes(x = Isolation, y = AvgDiff))+geom_point(pch = 19, fill = "grey60")
a = a + geom_errorbar(aes(ymin=AvgDiff-VarDiff, ymax=AvgDiff+VarDiff))
a = a + geom_hline(yintercept = 0, color = "red", linetype  = 2)
a = a + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"),axis.text.x = element_text(angle = 90))  
a = a + xlab("Isolation") + ylab("Average Difference")
ggsave(file = paste("All_IsoStat_NN_", Sys.Date(), ".pdf", sep = ""))

a = ggplot(All_IsoStat_NN_df, aes(x = Isolation, y = AvgTraitCounts))+geom_point(pch = 19, fill = "grey60")
a = a + geom_errorbar(aes(ymin=AvgTraitCounts-VarTraitCounts, ymax=AvgTraitCounts+VarTraitCounts))
a = a + geom_hline(yintercept = 0, color = "red", linetype  = 2)
a = a + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"),axis.text.x = element_text(angle = 90))  
a = a + xlab("Isolation") + ylab("Average Trait Counts")
ggsave(paste("All_TraitCount_NN_", Sys.Date(), ".pdf", sep = ""))

###### Isolation Trait pair analysis for Negatively Associated traits (+,-) ######
# Determine the number of observed (+/-)'s for Trait pairs in an isolation environment for the observed data and compare it to the permuted data # 
IsoAnalysis_NP_df = data.frame(Isolation = character(),
                               Trait_A = character(),
                               Trait_B = character(),
                               Observed_PN = numeric(), 
                               Percent_PN = numeric())

k = 1
for(iso in 1:length(Isolations)){
  tempObs = IsoAllAssoc_list[[iso]]$Diff_PN
  for(i in 1:(length(Traits)-1)){
    for(j in (i+1):length(Traits)){
      temp_list = list()
      for(perm in 1:1000){
        Val1 = PermutedIsoAllAssoc_list[[perm]][[iso]]$Diff_PN[i,j]
        Val2 = PermutedIsoAllAssoc_list[[perm]][[iso]]$Diff_PN[j,i]
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

# Summary of permuted data trait pairs within an isolation environment - Negative (+,-) #
PermutedIsoSum_PN_df = data.frame(Isolation = character(),
                                  Trait_A = character(),
                                  Trait_B = character(),
                                  Mean_PN = numeric(), 
                                  SD_PN = numeric(), 
                                  Var_PN = numeric())

k = 1
for(iso in 1:length(Isolations)){
  for(i in 1:(length(Traits)-1)){
    for(j in (i+1):length(Traits)){
      temp_list = list()
      for(perm in 1:1000){
        Val1 = PermutedIsoAllAssoc_list[[perm]][[iso]]$Diff_PN[i,j]
        Val2 = PermutedIsoAllAssoc_list[[perm]][[iso]]$Diff_PN[j,i]
        if(Val1 != Val2 & Val1 > Val2){
          temp_list[perm] = Val1
        }else{ 
          temp_list[perm] = Val2
        }
      }
      templist = unlist(temp_list)
      PermutedIsoSum_PN_df[k,1] = Isolations[iso]
      PermutedIsoSum_PN_df[k,2] = Traits[i]
      PermutedIsoSum_PN_df[k,3] = Traits[j]
      PermutedIsoSum_PN_df[k,4] = mean(templist, na.rm = T)
      PermutedIsoSum_PN_df[k,5] =sd(templist, na.rm = T)
      PermutedIsoSum_PN_df[k,6] =var(templist, na.rm = T)
      k = k+1
    }
  } 
}

write.csv(PermutedIsoSum_PN_df, file = paste("PermutedIsoSum_PN_df_", Sys.Date(), ".csv", sep = ""), row.names = F)

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

####### Isolation Trait Pair Visualization - Positive (+,+) ####### 
AllIsoSum_PN_df = merge(IsoAnalysis_PN_df, PermutedIsoSum_PN_df)
AllIsoSum_PN_df$Difference = AllIsoSum_PN_df$Observed_PN - AllIsoSum_PN_df$Mean_PN
AllIsoSum_PN_df$Color = rep(0, nrow(AllIsoSum_PN_df))

for(i in 1:nrow(AllIsoSum_PN_df)){
  SigValue = NegIsoCutOffs_PN_df[which(NegIsoCutOffs_PN_df$Isolation == AllIsoSum_PN_df[i,"Isolation"]), 2] 
  if(AllIsoSum_PN_df[i,"Difference"] > 0 & AllIsoSum_PN_df[i,"Percent_PN"] >= SigValue){
    AllIsoSum_PN_df[i, "Color"] = 1
  }else if(AllIsoSum_PN_df[i,"Difference"] < 0 & AllIsoSum_PN_df[i,"Percent_PN"] >= SigValue){
    AllIsoSum_PN_df[i, "Color"] = 2
  }else{
    AllIsoSum_PN_df[i, "Color"] = 3
  }
}

write.csv(AllIsoSum_PN_df, file = paste("AllIsoSum_PN_df_", Sys.Date(), ".csv", sep = ""), row.names = F)

pdf(paste("Isolation_TraitPair_PN_", Sys.Date(), ".pdf", sep = ""))
for(iso in 1:length(Isolations)){
  tempDF = AllIsoSum_PN_df[which(AllIsoSum_PN_df$Isolation == Isolations[iso]),]
  lengthCheck = length(unique(tempDF$Color))
  Values = unique(tempDF$Color)
  if(lengthCheck == 3){
    colorValues = c("#5C0023", "#034A58", "grey80")
  }else if(lengthCheck == 2 & length(which(Values == 1) > 0)){
    colorValues = c("#5C0023", "grey80")
  }else if(lengthCheck == 2 & length(which(Values == 2) > 0)){
    colorValues = c("#034A58", "grey80")
  }
  Mean_PN_Diff = mean(tempDF$Difference)
  
  a = ggplot(tempDF, aes(x = Trait_A, y = Difference))+geom_point()+scale_color_manual(values = colorValues)+geom_hline(yintercept = Mean_PN_Diff, color = "red")
  a = a + theme_bw()+ ggtitle(Isolations[iso]) + xlab("Trait Pairs")+ylab("Count")+scale_fill_grey()+ theme(axis.title.x = element_text(size = 12), axis.text.x = element_text(size = 0, angle = 90),axis.text.y = element_text(size = 8))
  print(a)
}
dev.off()

IsoPairStat_PN_df = data.frame(Isolation = character(), AvgDiff = numeric(), VarDiff = numeric(), AvgTraitCounts = numeric(), VarTraitCounts = numeric())
for(iso in 1:length(Isolations)){
  tempDF = AllIsoSum_PN_df[which(AllIsoSum_PN_df$Isolation == Isolations[iso]),]
  IsoPairStat_PN_df[iso, 1] = Isolations[iso]
  IsoPairStat_PN_df[iso, 2] = mean(tempDF$Difference, na.rm = TRUE)
  IsoPairStat_PN_df[iso, 3] = var(tempDF$Difference, na.rm = TRUE)
  IsoPairStat_PN_df[iso, 4] = mean(tempDF$Observed_PN, na.rm = TRUE)
  IsoPairStat_PN_df[iso, 5] = var(tempDF$Observed_PN, na.rm = TRUE)
}

write.csv(IsoPairStat_PN_df, file = paste("IsoPairStat_PN_df_", Sys.Date(), ".csv", sep = ""), row.names = F)

a = ggplot(IsoPairStat_PN_df, aes(x = Isolation, y = AvgDiff))+geom_point(pch = 19, fill = "grey60")
a = a + geom_errorbar(aes(ymin=AvgDiff-VarDiff, ymax=AvgDiff+VarDiff))
a = a + geom_hline(yintercept = 0, color = "red", linetype  = 2)
a = a + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"),axis.text.x = element_text(angle = 90))  
a = a + xlab("Isolation") + ylab("Average Difference")
ggsave(paste("IsoPairStat_PN_", Sys.Date(), ".pdf", sep = ""))

a = ggplot(IsoPairStat_PN_df, aes(x = Isolation, y = AvgTraitCounts))+geom_point(pch = 19, fill = "grey60")
a = a + geom_errorbar(aes(ymin=AvgTraitCounts-VarTraitCounts, ymax=AvgTraitCounts+VarTraitCounts))
a = a + geom_hline(yintercept = 0, color = "red", linetype  = 2)
a = a + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"),axis.text.x = element_text(angle = 90))  
a = a + xlab("Isolation") + ylab("Average Trait Counts")
ggsave(paste("IsoPair_TraitCount_PN_", Sys.Date(), ".pdf", sep = ""))


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

# Limit Trait pairs in the isolation environment to trait pairs where both traits are significantly associated with the isolation environment - Negative (-/+) #
Both_IsoSigSub_PN_df = matrix(0, nrow = 1 ,ncol = ncol(AllIsoSum_PN_df))
colnames(Both_IsoSigSub_PN_df) = colnames(AllIsoSum_PN_df)

Both_IsoSigSub_PN_df = data.frame(Both_IsoSigSub_PN_df)

for(i in 1:length(Neg_SigTraitsIso_PN_list)){
  iso = Neg_SigTraitsIso_PN_list[[i]]$Isolation
  keepTraits = Neg_SigTraitsIso_PN_list[[i]]$Traits
  tempData = AllIsoSum_PN_df[which(AllIsoSum_PN_df$Isolation == iso),]
  dataKeep = tempData[which(tempData$Trait_A %in% keepTraits),] 
  dataKeep = dataKeep[which(dataKeep$Trait_B %in% keepTraits),]
  Both_IsoSigSub_PN_df = rbind(Both_IsoSigSub_PN_df, dataKeep)
}

Both_IsoSigSub_PN_df = Both_IsoSigSub_PN_df[-1,]

write.csv(Both_IsoSigSub_PN_df, file = paste("Both_IsoSigSub_PN_df_", Sys.Date(), ".csv", sep = ""), row.names = F)

pdf(paste("BothSig_Isolation_TraitPair_PN_", Sys.Date(), ".pdf", sep = ""))
for(iso in 1:length(Isolations)){
  tempDF = Both_IsoSigSub_PN_df[which(Both_IsoSigSub_PN_df$Isolation == Isolations[iso]),]
  lengthCheck = length(unique(tempDF$Color))
  Values = unique(tempDF$Color)
  if(lengthCheck == 3){
    colorValues = c("#5C0023", "#034A58", "grey80")
  }else if(lengthCheck == 2 & length(which(Values == 1) > 0)){
    colorValues = c("#5C0023", "grey80")
  }else if(lengthCheck == 2 & length(which(Values == 2) > 0)){
    colorValues = c("#034A58", "grey80")
  }
  Both_Mean_PN_Diff = mean(tempDF$Difference)
  
  a = ggplot(tempDF, aes(x = Trait_A, y = Difference))+geom_point()+scale_color_manual(values = colorValues)+geom_hline(yintercept = Both_Mean_PN_Diff, color = "red")
  a = a + theme_bw()+ ggtitle(Isolations[iso]) + xlab("Trait Pairs")+ylab("Count")+scale_fill_grey()+ theme(axis.title.x = element_text(size = 12), axis.text.x = element_text(size = 0, angle = 90),axis.text.y = element_text(size = 8))
  print(a)
}
dev.off()

Both_IsoStat_PN_df = data.frame(Isolation = character(), AvgDiff = numeric(), VarDiff = numeric(), AvgTraitCounts = numeric(), VarTraitCounts = numeric())
for(iso in 1:length(Isolations)){
  tempDF = Both_IsoSigSub_PN_df[which(Both_IsoSigSub_PN_df$Isolation == Isolations[iso]),]
  Both_IsoStat_PN_df[iso, 1] = Isolations[iso]
  Both_IsoStat_PN_df[iso, 2] = mean(tempDF$Difference, na.rm = TRUE)
  Both_IsoStat_PN_df[iso, 3] = var(tempDF$Difference, na.rm = TRUE)
  Both_IsoStat_PN_df[iso, 4] = mean(tempDF$Observed_PN, na.rm = TRUE)
  Both_IsoStat_PN_df[iso, 5] = var(tempDF$Observed_PN, na.rm = TRUE)
}

write.csv(Both_IsoStat_PN_df, file = paste("Both_IsoStat_PN_df_", Sys.Date(), ".csv", sep = ""), row.names = F)

a = ggplot(Both_IsoStat_PN_df, aes(x = Isolation, y = AvgDiff))+geom_point(pch = 19, fill = "grey60")
a = a + geom_errorbar(aes(ymin=AvgDiff-VarDiff, ymax=AvgDiff+VarDiff))
a = a + geom_hline(yintercept = 0, color = "red", linetype  = 2)
a = a + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"),axis.text.x = element_text(angle = 90))  
a = a + xlab("Isolation") + ylab("Average Difference")
ggsave(file = paste("Both_IsoStat_PN_", Sys.Date(), ".pdf", sep = ""))

a = ggplot(Both_IsoStat_PN_df, aes(x = Isolation, y = AvgTraitCounts))+geom_point(pch = 19, fill = "grey60")
a = a + geom_errorbar(aes(ymin=AvgTraitCounts-VarTraitCounts, ymax=AvgTraitCounts+VarTraitCounts))
a = a + geom_hline(yintercept = 0, color = "red", linetype  = 2)
a = a + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"),axis.text.x = element_text(angle = 90))  
a = a + xlab("Isolation") + ylab("Average Trait Counts")
ggsave(paste("Both_TraitCount_PN_", Sys.Date(), ".pdf", sep = ""))

All_IsoSigSub_PN_df = matrix(0, nrow = 1 ,ncol = ncol(AllIsoSum_PN_df))
colnames(All_IsoSigSub_PN_df) = colnames(AllIsoSum_PN_df)

All_IsoSigSub_PN_df = data.frame(All_IsoSigSub_PN_df)

for(i in 1:length(Neg_SigTraitsIso_PN_list)){
  iso = Neg_SigTraitsIso_PN_list[[i]]$Isolation
  keepTraits = Neg_SigTraitsIso_PN_list[[i]]$Traits
  tempData = AllIsoSum_PN_df[which(AllIsoSum_PN_df$Isolation == iso),]
  dataKeep = tempData[which(tempData$Trait_A %in% keepTraits | tempData$Trait_B %in% keepTraits),]
  All_IsoSigSub_PN_df = rbind(All_IsoSigSub_PN_df, dataKeep)
}

All_IsoSigSub_PN_df = All_IsoSigSub_PN_df[-1,]

write.csv(All_IsoSigSub_PN_df, file = paste("All_IsoSigSub_PN_df", Sys.Date(), ".csv", sep = ""), row.names = F)

pdf(paste("AllSig_Isolation_TraitPair_PN_", Sys.Date(), ".pdf", sep = ""))
for(iso in 1:length(Isolations)){
  tempDF = All_IsoSigSub_PN_df[which(All_IsoSigSub_PN_df$Isolation == Isolations[iso]),]
  lengthCheck = length(unique(tempDF$Color))
  Values = unique(tempDF$Color)
  if(lengthCheck == 3){
    colorValues = c("#5C0023", "#034A58", "grey80")
  }else if(lengthCheck == 2 & length(which(Values == 1) > 0)){
    colorValues = c("#5C0023", "grey80")
  }else if(lengthCheck == 2 & length(which(Values == 2) > 0)){
    colorValues = c("#034A58", "grey80")
  }
  All_Mean_PN_Diff = mean(tempDF$Difference)
  a = ggplot(tempDF, aes(x = Trait_A, y = Difference))+geom_point()+scale_color_manual(values = colorValues)+geom_hline(yintercept = All_Mean_PN_Diff, color = "red")
  a = a + theme_bw()+ ggtitle(Isolations[iso]) + xlab("Trait Pairs")+ylab("Count")+scale_fill_grey()+ theme(axis.title.x = element_text(size = 12), axis.text.x = element_text(size = 0, angle = 90),axis.text.y = element_text(size = 8))
  print(a)
}
dev.off()

All_IsoStat_PN_df = data.frame(Isolation = character(), AvgDiff = numeric(), VarDiff = numeric(), AvgTraitCounts = numeric(), VarTraitCounts = numeric())
for(iso in 1:length(Isolations)){
  tempDF = All_IsoSigSub_PN_df[which(All_IsoSigSub_PN_df$Isolation == Isolations[iso]),]
  All_IsoStat_PN_df[iso, 1] = Isolations[iso]
  All_IsoStat_PN_df[iso, 2] = mean(tempDF$Difference, na.rm = TRUE)
  All_IsoStat_PN_df[iso, 3] = var(tempDF$Difference, na.rm = TRUE)
  All_IsoStat_PN_df[iso, 4] = mean(tempDF$Observed_PN, na.rm = TRUE)
  All_IsoStat_PN_df[iso, 5] = var(tempDF$Observed_PN, na.rm = TRUE)
}

write.csv(All_IsoStat_PN_df, file = paste("All_IsoStat_PN_df_", Sys.Date(), ".csv", sep = ""), row.names = F)

a = ggplot(All_IsoStat_PN_df, aes(x = Isolation, y = AvgDiff))+geom_point(pch = 19, fill = "grey60")
a = a + geom_errorbar(aes(ymin=AvgDiff-VarDiff, ymax=AvgDiff+VarDiff))
a = a + geom_hline(yintercept = 0, color = "red", linetype  = 2)
a = a + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"),axis.text.x = element_text(angle = 90))  
a = a + xlab("Isolation") + ylab("Average Difference")
ggsave(file = paste("All_IsoStat_PN_", Sys.Date(), ".pdf", sep = ""))

a = ggplot(All_IsoStat_PN_df, aes(x = Isolation, y = AvgTraitCounts))+geom_point(pch = 19, fill = "grey60")
a = a + geom_errorbar(aes(ymin=AvgTraitCounts-VarTraitCounts, ymax=AvgTraitCounts+VarTraitCounts))
a = a + geom_hline(yintercept = 0, color = "red", linetype  = 2)
a = a + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"),axis.text.x = element_text(angle = 90))  
a = a + xlab("Isolation") + ylab("Average Trait Counts")
ggsave(paste("All_TraitCount_PN_", Sys.Date(), ".pdf", sep = ""))

###### Isolation Trait pair analysis for Negatively Associated traits (-,+) ######
# Determine the number of observed (-/+)'s for Trait pairs in an isolation environment for the observed data and compare it to the permuted data # 
IsoAnalysis_NP_df = data.frame(Isolation = character(),
                               Trait_A = character(),
                               Trait_B = character(),
                               Observed_NP = numeric(), 
                               Percent_NP = numeric())

k = 1
for(iso in 1:length(Isolations)){
  tempObs = IsoAllAssoc_list[[iso]]$Diff_NP
  for(i in 1:(length(Traits)-1)){
    for(j in (i+1):length(Traits)){
      temp_list = list()
      for(perm in 1:1000){
        Val1 = PermutedIsoAllAssoc_list[[perm]][[iso]]$Diff_NP[i,j]
        Val2 = PermutedIsoAllAssoc_list[[perm]][[iso]]$Diff_NP[j,i]
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

# Summary of permuted data trait pairs within an isolation environment - Postive (+,+) #
PermutedIsoSum_NP_df = data.frame(Isolation = character(),
                                  Trait_A = character(),
                                  Trait_B = character(),
                                  Mean_NP = numeric(), 
                                  SD_NP = numeric(), 
                                  Var_NP = numeric())

k = 1
for(iso in 1:length(Isolations)){
  for(i in 1:(length(Traits)-1)){
    for(j in (i+1):length(Traits)){
      temp_list = list()
      for(perm in 1:1000){
        Val1 = PermutedIsoAllAssoc_list[[perm]][[iso]]$Diff_NP[i,j]
        Val2 = PermutedIsoAllAssoc_list[[perm]][[iso]]$Diff_NP[j,i]
        if(Val1 != Val2 & Val1 > Val2){
          temp_list[perm] = Val1
        }else{ 
          temp_list[perm] = Val2
        }
      }
      templist = unlist(temp_list)
      PermutedIsoSum_NP_df[k,1] = Isolations[iso]
      PermutedIsoSum_NP_df[k,2] = Traits[i]
      PermutedIsoSum_NP_df[k,3] = Traits[j]
      PermutedIsoSum_NP_df[k,4] = mean(templist, na.rm = T)
      PermutedIsoSum_NP_df[k,5] =sd(templist, na.rm = T)
      PermutedIsoSum_NP_df[k,6] =var(templist, na.rm = T)
      k = k+1
    }
  } 
}

write.csv(PermutedIsoSum_NP_df, file = paste("PermutedIsoSum_NP_df_", Sys.Date(), ".csv", sep = ""), row.names = F)

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

####### Isolation Trait Pair Visualization - Positive (+,+) ####### 
AllIsoSum_NP_df = merge(IsoAnalysis_NP_df, PermutedIsoSum_NP_df)
AllIsoSum_NP_df$Difference = AllIsoSum_NP_df$Observed_NP - AllIsoSum_NP_df$Mean_NP
AllIsoSum_NP_df$Color = rep(0, nrow(AllIsoSum_NP_df))

for(i in 1:nrow(AllIsoSum_NP_df)){
  SigValue = NegIsoCutOffs_NP_df[which(NegIsoCutOffs_NP_df$Isolation == AllIsoSum_NP_df[i,"Isolation"]), 2] 
  if(AllIsoSum_NP_df[i,"Difference"] > 0 & AllIsoSum_NP_df[i,"Percent_NP"] >= SigValue){
    AllIsoSum_NP_df[i, "Color"] = 1
  }else if(AllIsoSum_NP_df[i,"Difference"] < 0 & AllIsoSum_NP_df[i,"Percent_NP"] >= SigValue){
    AllIsoSum_NP_df[i, "Color"] = 2
  }else{
    AllIsoSum_NP_df[i, "Color"] = 3
  }
}
write.csv(AllIsoSum_NP_df, file = paste("AllIsoSum_NP_df_", Sys.Date(), ".csv", sep = ""), row.names = F)

pdf(paste("Isolation_TraitPair_NP_", Sys.Date(), ".pdf", sep = ""))
for(iso in 1:length(Isolations)){
  tempDF = AllIsoSum_NP_df[which(AllIsoSum_NP_df$Isolation == Isolations[iso]),]
  lengthCheck = length(unique(tempDF$Color))
  Values = unique(tempDF$Color)
  if(lengthCheck == 3){
    colorValues = c("#5C0023", "#034A58", "grey80")
  }else if(lengthCheck == 2 & length(which(Values == 1) > 0)){
    colorValues = c("#5C0023", "grey80")
  }else if(lengthCheck == 2 & length(which(Values == 2) > 0)){
    colorValues = c("#034A58", "grey80")
  }
  Mean_NP_Diff = mean(tempDF$Difference)
  a = ggplot(tempDF, aes(x = Trait_A, y = Difference))+geom_point()+scale_color_manual(values = colorValues)+geom_hline(yintercept = Mean_NP_Diff, color = "red")
  a = a + theme_bw()+ ggtitle(Isolations[iso]) + xlab("Trait Pairs")+ylab("Count")+scale_fill_grey()+ theme(axis.title.x = element_text(size = 12), axis.text.x = element_text(size = 0, angle = 90),axis.text.y = element_text(size = 8))
  print(a)
}
dev.off()

AllIsoTraitPairStat_NP_df = data.frame(Isolation = character(), AvgDiff = numeric(), VarDiff = numeric(), AvgTraitCounts = numeric(), VarTraitCounts = numeric())
for(iso in 1:length(Isolations)){
  tempDF = AllIsoSum_NP_df[which(AllIsoSum_NP_df$Isolation == Isolations[iso]),]
  AllIsoTraitPairStat_NP_df[iso, 1] = Isolations[iso]
  AllIsoTraitPairStat_NP_df[iso, 2] = mean(tempDF$Difference, na.rm = TRUE)
  AllIsoTraitPairStat_NP_df[iso, 3] = var(tempDF$Difference, na.rm = TRUE)
  AllIsoTraitPairStat_NP_df[iso, 4] = mean(tempDF$Observed_NP, na.rm = TRUE)
  AllIsoTraitPairStat_NP_df[iso, 5] = var(tempDF$Observed_NP, na.rm = TRUE)
}

write.csv(AllIsoTraitPairStat_NP_df, file = paste("AllIsoTraitPairStat_NP_df_", Sys.Date(), ".csv", sep = ""), row.names = F)

a = ggplot(AllIsoTraitPairStat_NP_df, aes(x = Isolation, y = AvgDiff))+geom_point(pch = 19, fill = "grey60")
a = a + geom_errorbar(aes(ymin=AvgDiff-VarDiff, ymax=AvgDiff+VarDiff))
a = a + geom_hline(yintercept = 0, color = "red", linetype  = 2)
a = a + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"),axis.text.x = element_text(angle = 90))  
a = a + xlab("Isolation") + ylab("Average Difference")
ggsave(paste("AllIsoTraitPairStat_NP_", Sys.Date(), ".pdf", sep = ""))

a = ggplot(AllIsoTraitPairStat_NP_df, aes(x = Isolation, y = AvgTraitCounts))+geom_point(pch = 19, fill = "grey60")
a = a + geom_errorbar(aes(ymin=AvgTraitCounts-VarTraitCounts, ymax=AvgTraitCounts+VarTraitCounts))
a = a + geom_hline(yintercept = 0, color = "red", linetype  = 2)
a = a + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"),axis.text.x = element_text(angle = 90))  
a = a + xlab("Isolation") + ylab("Average Trait Counts")
ggsave(paste("IsoTrait_TraitCount_NP_", Sys.Date(), ".pdf", sep = ""))

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

# Limit Trait pairs in the isolation environment to trait pairs where both traits are significantly associated with the isolation environment - Negative (-/+) #
Both_IsoSigSub_NP_df = matrix(0, nrow = 1 ,ncol = ncol(AllIsoSum_NP_df))
colnames(Both_IsoSigSub_NP_df) = colnames(AllIsoSum_NP_df)

Both_IsoSigSub_NP_df = data.frame(Both_IsoSigSub_NP_df)

for(i in 1:length(Neg_SigTraitsIso_NP_list)){
  iso = Neg_SigTraitsIso_NP_list[[i]]$Isolation
  keepTraits = Neg_SigTraitsIso_NP_list[[i]]$Traits
  tempData = AllIsoSum_NP_df[which(AllIsoSum_NP_df$Isolation == iso),]
  dataKeep = tempData[which(tempData$Trait_A %in% keepTraits),] 
  dataKeep = dataKeep[which(dataKeep$Trait_B %in% keepTraits),]
  Both_IsoSigSub_NP_df = rbind(Both_IsoSigSub_NP_df, dataKeep)
}

Both_IsoSigSub_NP_df = Both_IsoSigSub_NP_df[-1,]

write.csv(Both_IsoSigSub_NP_df, file = paste("Both_IsoSigSub_NP_df_", Sys.Date(), ".csv", sep = ""), row.names = F)


pdf(paste("BothSig_Isolation_TraitPair_NP_", Sys.Date(), ".pdf", sep = ""))
for(iso in 1:length(Isolations)){
  tempDF = Both_IsoSigSub_NP_df[which(Both_IsoSigSub_NP_df$Isolation == Isolations[iso]),]
  lengthCheck = length(unique(tempDF$Color))
  Values = unique(tempDF$Color)
  if(lengthCheck == 3){
    colorValues = c("#5C0023", "#034A58", "grey80")
  }else if(lengthCheck == 2 & length(which(Values == 1) > 0)){
    colorValues = c("#5C0023", "grey80")
  }else if(lengthCheck == 2 & length(which(Values == 2) > 0)){
    colorValues = c("#034A58", "grey80")
  }
  Both_Mean_NP_Diff = mean(tempDF$Difference)
  a = ggplot(tempDF, aes(x = Trait_A, y = Difference))+geom_point()+scale_color_manual(values = colorValues)+geom_hline(yintercept = Both_Mean_NP_Diff, color = "red")
  a = a + theme_bw()+ ggtitle(Isolations[iso]) + xlab("Trait Pairs")+ylab("Count")+scale_fill_grey()+ theme(axis.title.x = element_text(size = 12), axis.text.x = element_text(size = 0, angle = 90),axis.text.y = element_text(size = 8))
  print(a)
}
dev.off()

Both_IsoStat_NP_df = data.frame(Isolation = character(), AvgDiff = numeric(), VarDiff = numeric(), AvgTraitCounts = numeric(), VarTraitCounts = numeric())
for(iso in 1:length(Isolations)){
  tempDF = Both_IsoSigSub_NP_df[which(Both_IsoSigSub_NP_df$Isolation == Isolations[iso]),]
  Both_IsoStat_NP_df[iso, 1] = Isolations[iso]
  Both_IsoStat_NP_df[iso, 2] = mean(tempDF$Difference, na.rm = TRUE)
  Both_IsoStat_NP_df[iso, 3] = var(tempDF$Difference, na.rm = TRUE)
  Both_IsoStat_NP_df[iso, 4] = mean(tempDF$Observed_NP, na.rm = TRUE)
  Both_IsoStat_NP_df[iso, 5] = var(tempDF$Observed_NP, na.rm = TRUE)
}

write.csv(Both_IsoStat_NP_df, file = paste("Both_IsoStat_NP_df_", Sys.Date(), ".csv", sep = ""), row.names = F)

a = ggplot(Both_IsoStat_NP_df, aes(x = Isolation, y = AvgDiff))+geom_point(pch = 19, fill = "grey60")
a = a + geom_errorbar(aes(ymin=AvgDiff-VarDiff, ymax=AvgDiff+VarDiff))
a = a + geom_hline(yintercept = 0, color = "red", linetype  = 2)
a = a + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"),axis.text.x = element_text(angle = 90))  
a = a + xlab("Isolation") + ylab("Average Difference")
ggsave(file = paste("Both_IsoStat_NP_", Sys.Date(), ".pdf", sep = ""))

a = ggplot(Both_IsoStat_NP_df, aes(x = Isolation, y = AvgTraitCounts))+geom_point(pch = 19, fill = "grey60")
a = a + geom_errorbar(aes(ymin=AvgTraitCounts-VarTraitCounts, ymax=AvgTraitCounts+VarTraitCounts))
a = a + geom_hline(yintercept = 0, color = "red", linetype  = 2)
a = a + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"),axis.text.x = element_text(angle = 90))  
a = a + xlab("Isolation") + ylab("Average Trait Counts")
ggsave(paste("Both_TraitCount_NP_", Sys.Date(), ".pdf", sep = ""))

All_IsoSigSub_NP_df = matrix(0, nrow = 1 ,ncol = ncol(AllIsoSum_NP_df))
colnames(All_IsoSigSub_NP_df) = colnames(AllIsoSum_NP_df)

All_IsoSigSub_NP_df = data.frame(All_IsoSigSub_NP_df)

for(i in 1:length(Neg_SigTraitsIso_NP_list)){
  iso = Neg_SigTraitsIso_NP_list[[i]]$Isolation
  keepTraits = Neg_SigTraitsIso_NP_list[[i]]$Traits
  tempData = AllIsoSum_NP_df[which(AllIsoSum_NP_df$Isolation == iso),]
  dataKeep = tempData[which(tempData$Trait_A %in% keepTraits | tempData$Trait_B %in% keepTraits),]
  All_IsoSigSub_NP_df = rbind(All_IsoSigSub_NP_df, dataKeep)
}

All_IsoSigSub_NP_df = All_IsoSigSub_NP_df[-1,]

write.csv(All_IsoSigSub_NP_df, file = paste("All_IsoSigSub_NP_df", Sys.Date(), ".csv", sep = ""), row.names = F)

pdf(paste("AllSig_Isolation_TraitPair_NP_", Sys.Date(), ".pdf", sep = ""))
for(iso in 1:length(Isolations)){
  tempDF = All_IsoSigSub_NP_df[which(All_IsoSigSub_NP_df$Isolation == Isolations[iso]),]
  lengthCheck = length(unique(tempDF$Color))
  Values = unique(tempDF$Color)
  if(lengthCheck == 3){
    colorValues = c("#5C0023", "#034A58", "grey80")
  }else if(lengthCheck == 2 & length(which(Values == 1) > 0)){
    colorValues = c("#5C0023", "grey80")
  }else if(lengthCheck == 2 & length(which(Values == 2) > 0)){
    colorValues = c("#034A58", "grey80")
  }
  All_Mean_NP_Diff = mean(tempDF$Difference)
  
  a = ggplot(tempDF, aes(x = Trait_A, y = Difference))+geom_point()+scale_color_manual(values = colorValues)+geom_hline(yintercept = All_Mean_NP_Diff, color = "Red")
  a = a + theme_bw()+ ggtitle(Isolations[iso]) + xlab("Trait Pairs")+ylab("Count")+scale_fill_grey()+ theme(axis.title.x = element_text(size = 12), axis.text.x = element_text(size = 0, angle = 90),axis.text.y = element_text(size = 8))
  print(a)
}
dev.off()

All_IsoStat_NP_df = data.frame(Isolation = character(), AvgDiff = numeric(), VarDiff = numeric(), AvgTraitCounts = numeric(), VarTraitCounts = numeric())
for(iso in 1:length(Isolations)){
  tempDF = All_IsoSigSub_NP_df[which(All_IsoSigSub_NP_df$Isolation == Isolations[iso]),]
  All_IsoStat_NP_df[iso, 1] = Isolations[iso]
  All_IsoStat_NP_df[iso, 2] = mean(tempDF$Difference, na.rm = TRUE)
  All_IsoStat_NP_df[iso, 3] = var(tempDF$Difference, na.rm = TRUE)
  All_IsoStat_NP_df[iso, 4] = mean(tempDF$Observed_NP, na.rm = TRUE)
  All_IsoStat_NP_df[iso, 5] = var(tempDF$Observed_NP, na.rm = TRUE)
}

write.csv(All_IsoStat_NP_df, file = paste("All_IsoStat_NP_df_", Sys.Date(), ".csv", sep = ""), row.names = F)

a = ggplot(All_IsoStat_NP_df, aes(x = Isolation, y = AvgDiff))+geom_point(pch = 19, fill = "grey60")
a = a + geom_errorbar(aes(ymin=AvgDiff-VarDiff, ymax=AvgDiff+VarDiff))
a = a + geom_hline(yintercept = 0, color = "red", linetype  = 2)
a = a + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"),axis.text.x = element_text(angle = 90))  
a = a + xlab("Isolation") + ylab("Average Difference")
ggsave(file = paste("All_IsoStat_NP_", Sys.Date(), ".pdf", sep = ""))

a = ggplot(All_IsoStat_NP_df, aes(x = Isolation, y = AvgTraitCounts))+geom_point(pch = 19, fill = "grey60")
a = a + geom_errorbar(aes(ymin=AvgTraitCounts-VarTraitCounts, ymax=AvgTraitCounts+VarTraitCounts))
a = a + geom_hline(yintercept = 0, color = "red", linetype  = 2)
a = a + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"),axis.text.x = element_text(angle = 90))  
a = a + xlab("Isolation") + ylab("Average Trait Counts")
ggsave(paste("All_TraitCount_NP_", Sys.Date(), ".pdf", sep = ""))


PermutedIsoSum_Neg_df = data.frame(Isolation = character(),
                                  Trait_A = character(),
                                  Trait_B = character(),
                                  Mean_Neg = numeric(), 
                                  SD_Neg = numeric(), 
                                  Var_Neg = numeric())

k = 1
for(iso in 1:length(Isolations)){
  for(i in 1:(length(Traits)-1)){
    for(j in (i+1):length(Traits)){
      temp_lista = list()
      temp_listb = list()
      for(perm in 1:1000){
        Val1a = PermutedIsoAllAssoc_list[[perm]][[iso]]$Diff_PN[i,j]
        Val2a = PermutedIsoAllAssoc_list[[perm]][[iso]]$Diff_PN[j,i]
        if(Val1a != Val2a & Val1a > Val2a){
          temp_lista[perm] = Val1a
        }else{ 
          temp_lista[perm] = Val2a
        }
      templista = unlist(temp_lista)
      
      Val1b = PermutedIsoAllAssoc_list[[perm]][[iso]]$Diff_NP[i,j]
      Val2b = PermutedIsoAllAssoc_list[[perm]][[iso]]$Diff_NP[j,i]
      if(Val1b != Val2b & Val1b > Val2b){
        temp_listb[perm] = Val1b
      }else{ 
        temp_listb[perm] = Val2b
      }
    }
    templistb = unlist(temp_listb)
    templist = c(templista, templistb)
      PermutedIsoSum_Neg_df[k,1] = Isolations[iso]
      PermutedIsoSum_Neg_df[k,2] = Traits[i]
      PermutedIsoSum_Neg_df[k,3] = Traits[j]
      PermutedIsoSum_Neg_df[k,4] = mean(templist, na.rm = T)
      PermutedIsoSum_Neg_df[k,5] =sd(templist, na.rm = T)
      PermutedIsoSum_Neg_df[k,6] =var(templist, na.rm = T)
      k = k+1
    }
  } 
  print(iso)
}
save(list = ls(), file = "AllData.csv")
write.csv(PermutedIsoSum_Neg_df, file = paste("PermutedIsoSum_Neg_df_", Sys.Date(), ".csv", sep = ""), row.names = F)


All_IsoSigSub_Neg_df = merge(All_IsoSigSub_NP_df, All_IsoSigSub_PN_df, by = c("Isolation", "Trait_A", "Trait_B"))
All_IsoSigSub_Neg_df$NegObs = All_IsoSigSub_Neg_df$Observed_NP + All_IsoSigSub_Neg_df$Observed_PN
All_IsoSigSub_Neg_df = All_IsoSigSub_Neg_df[,c(1, 2, 3, ncol(All_IsoSigSub_Neg_df))]

All_PooledIsoSigSub_Neg_df = merge(All_IsoSigSub_Neg_df, PermutedIsoSum_Neg_df, by = c("Isolation", "Trait_A", "Trait_B"))
All_PooledIsoSigSub_Neg_df$Difference = All_PooledIsoSigSub_Neg_df$NegObs - All_PooledIsoSigSub_Neg_df$Mean_Neg

All_IsoStat_Neg_df = data.frame(Isolation = character(), AvgDiff = numeric(), VarDiff = numeric(), AvgTraitCounts = numeric(), VarTraitCounts = numeric())
for(iso in 1:length(Isolations)){
  tempDF = All_PooledIsoSigSub_Neg_df[which(All_PooledIsoSigSub_Neg_df$Isolation == Isolations[iso]),]
  All_IsoStat_Neg_df[iso, 1] = Isolations[iso]
  All_IsoStat_Neg_df[iso, 2] = mean(tempDF$Difference, na.rm = TRUE)
  All_IsoStat_Neg_df[iso, 3] = var(tempDF$Difference, na.rm = TRUE)
  All_IsoStat_Neg_df[iso, 4] = mean(tempDF$Mean_Neg, na.rm = TRUE)
  All_IsoStat_Neg_df[iso, 5] = var(tempDF$Mean_Neg, na.rm = TRUE)
}

write.csv(All_PooledIsoSigSub_Neg_df, file = paste("All_IsoStat_NP_df_", Sys.Date(), ".csv", sep = ""), row.names = F)
write.csv(All_IsoStat_Neg_df, file = paste("All_IsoStat_Neg_df_", Sys.Date(), ".csv", sep = ""), row.names = F)

a = ggplot(All_IsoStat_Neg_df, aes(x = Isolation, y = AvgDiff))+geom_point(pch = 19, fill = "grey60")
a = a + geom_errorbar(aes(ymin=AvgDiff-VarDiff, ymax=AvgDiff+VarDiff))
a = a + geom_hline(yintercept = 0, color = "red", linetype  = 2)
a = a + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"),axis.text.x = element_text(angle = 90))  
a = a + xlab("Isolation") + ylab("Average Difference")
ggsave(file = paste("All_IsoStat_Neg_", Sys.Date(), ".pdf", sep = ""))
