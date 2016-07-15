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

Isolations <- rownames(envt)
PermutedTraitxTrait_IsoList = list()
for(iter in 42:1000){
  PermutedTraitxTrait_IsoList[[iter]] = list()
  tempMatrix <- randomizeMatrix(traits, null.model="independentswap", iterations=1000)
  CWM<-functcomp(tempMatrix, envt, CWM.type="all")
  tempMatrix = data.frame(tempMatrix)
  tempMatrix$Species = rownames(tempMatrix)
  CompiledData = merge(envs2, tempMatrix, by = "Species")
  PermutedTraitxTrait_IsoList[[iter]]$Permutation <- iter 
  PermutedTraitxTrait_IsoList[[iter]]$CWM <- CWM 
  TraitxTrait_IsoList = list()
  CompiledData= CompiledData[,-which(colnames(CompiledData) == "Species")]
  k = 1
  for(i in 1:length(Isolations)){
    TraitIsolation_m = matrix(0, ncol = length(CompiledData[56:ncol(CompiledData)]))
    colnames(TraitIsolation_m) = colnames(CompiledData)[56:ncol(CompiledData)]
    TraitIsolation_m = data.frame(TraitIsolation_m, check.names = F)
    TempCol = which(colnames(CompiledData) == Isolations[i])
    TempData = CompiledData[which(CompiledData[,TempCol] == 1),57:ncol(CompiledData)]
    PermutedTraitxTrait_IsoList[[iter]]$TraitxTrait_IsoList[[k]] = list()
    PermutedTraitxTrait_IsoList[[iter]]$TraitxTrait_IsoList[[k]]$Isolation = Isolations[i]
    PermutedTraitxTrait_IsoList[[iter]]$TraitxTrait_IsoList[[k]]$TraitMatrix = TempData      
    PermutedTraitxTrait_IsoList[[iter]]$TraitByTraitIso_co = cooccur(mat = t(TempData), type = "spp_site", thresh = T, spp_names = T, eff_matrix = T, only_effects = F)
    PermutedTraitxTrait_IsoList[[iter]]$TraitxTrait_IsoList[[k]]$CoOcc = TraitByTraitIso_co[[2]]
    PermutedTraitxTrait_IsoList[[iter]]$TraitxTrait_IsoList[[k]]$EffectSize = effect.sizes(TraitByTraitIso_co)
    k = k+1
  }
  print(iter)
}
