# R script for the analysis of the genomic distance to redundancy with MAGs of the Lake Dziani Dzaha
# Date of creation : 21 February 2025
# Last modified : 16 June 2026
# Author : Perrine Cruaud

# Loading packages --------------
library(stringr)

# Set working directory -----------
setwd("/path/to/my/directory")

# Importing data -------------
Results_AllBins_KO <- read.csv(file="Resultats_AllBins_KO_WithGeneNames.csv", sep="\t", header=TRUE,
                                 row.names=1)
CrisisImpact_PerMAG <- read.csv(file="Table_Effet_Crise_matriceAAI.csv", sep=";", row.names=1) 
AAIMatrix <- read.csv(file="matrice_AAI_All_good_MAGs_cor.csv", sep=";", header=TRUE, row.names=1)
KOList_MMDE <- read.table(file="Liste_KO_MMDE.txt")


# Data formatting ---------------------
# Data simplification : only the necessary information are retained for the dataframe Results_AllBins_KO 
Results_AllBins_KO[1:10,1:10]
Results_AllBins_KO_Simple <- Results_AllBins_KO[,-c(1,2)] # Elimination of the first 2 columns =>  "definition" "knum.1" 
Results_AllBins_KO_Simple_tr <- as.data.frame(t(Results_AllBins_KO_Simple))
Results_AllBins_KO_Simple_tr[1:10,1:10]


# Data analysis ---------------------------------

# KO elimination if duplicates are found in the list
KOList_Temp <- KOList_MMDE$V1
KOList <- unique(KOList_Temp) 

length(KOList_Temp)
length(unique(KOList_Temp))




## Creating a filtered KO list by eliminating KOs without bins and without double effects -----------------------------

KOList_Filtered <- c(NULL)
KOList_WithoutBin <- c(NULL)
KOList_OneEffect_BeforeAfter <- c(NULL)
KOTable_OneEffect_BeforeAfter <- data.frame(KO_ID=NA,MaxAvg=NA,MaxAvg_Without100=NA,Nb_EffectBefore=NA,Nb_WithoutEffect=NA,Nb_EffectAfter=NA)
KOTable_WithoutBin <-  data.frame(KO_ID=NA,MaxAvg=NA,MaxAvg_Without100=NA,Nb_EffectBefore=NA,Nb_WithoutEffect=NA,Nb_EffectAfter=NA)


# Search for KOs of the KOList in Results_AllBins_KO dataframe to have the corresponding bin for each KO
for (NumKO in 1:length(KOList)) {
  BinList_1KO <- row.names(Results_AllBins_KO_Simple_tr)[which(Results_AllBins_KO_Simple_tr[,which(colnames(Results_AllBins_KO_Simple_tr)==KOList[NumKO])] > 0)]
   
   ### 1 - Coding of the crisis effect (For each bin listed in BinList_1KO) ----------------
     Code_CrisisEffect_PerBin <- data.frame(BinID=NA,CodeEffect=NA) # Empty dataframe
  
      # Code :
      # 1 : Before crisis (first column of CrisisImpact_PerMAG)
      # 2 : Without any effect (second column of CrisisImpact_PerMAG)
      # 3 : After crisis (third column of CrisisImpact_PerMAG)
  
     for (Num in 1:length(BinList_1KO)) {
       NumBin_Code_CrisisEffect <- c(BinList_1KO[Num],which(CrisisImpact_PerMAG[which(row.names(CrisisImpact_PerMAG)==BinList_1KO[Num]),]==1))
       Code_CrisisEffect_PerBin[Num,] <- NumBin_Code_CrisisEffect
      }
      # Code_CrisisEffect_PerBin # For 1 KO => Bin ID (column 1) and crisis effect (column 2)

  
  
     
  ### 2 -  Elimination of bins with no match in CrisisImpact_PerMAG ---------------------------
  # NB : bins with no match in CrisisImpact_PerMAG have their ID in the second column of Code_CrisisEffect_PerBin => "BinXXX"
      if (sum(str_detect(Code_CrisisEffect_PerBin[,2],"Bin"))>=1) {
         Code_CrisisEffect_PerBin_WithInfoOfCrisisEffect <- Code_CrisisEffect_PerBin[-which(str_detect(Code_CrisisEffect_PerBin[,2],"Bin")),] # If no match => elimination
        } else {
         Code_CrisisEffect_PerBin_WithInfoOfCrisisEffect <- Code_CrisisEffect_PerBin 
       }
  
      ### Elimination of KOs with no bin (after the "Info for crisis effect" filtration step)
       if (nrow(Code_CrisisEffect_PerBin_WithInfoOfCrisisEffect)==0) {
          KOList_WithoutBin <- c(KOList_WithoutBin,KOList[NumKO])
          KOTable_WithoutBin <- rbind(KOTable_WithoutBin,c(KOList[NumKO],NA,NA,"NoBin","NoBin","NoBin")) 
        # If no bin in a KO => add in the KOList_WithoutBin and KOTable_WithoutBin, and restart for the next KO
        } else {

  
    
  ### 3 - Elimination of bins with no match in the AAI matrix ----------------------
      BinList_AAIMatrix <- c(colnames(AAIMatrix),row.names(AAIMatrix))
      BinList_WithMatch_AAIMatrix <- c(NULL)
  
      for (i in 1:length(Code_CrisisEffect_PerBin_WithInfoOfCrisisEffect$BinID)) {
        BinList_WithMatch_AAIMatrix[i] <- sum(BinList_AAIMatrix==Code_CrisisEffect_PerBin_WithInfoOfCrisisEffect$BinID[i])>0 # If >0 => bin with match in AAI Matrix => TRUE
        }
      Code_CrisisEffect_PerBin_WithInfoOfCrisisEffect_MatchInAAI <- Code_CrisisEffect_PerBin_WithInfoOfCrisisEffect[which(BinList_WithMatch_AAIMatrix),]
  
      ### Elimination of KO with no bin (after the "AAI matrix match" filtration step)
       if (nrow(Code_CrisisEffect_PerBin_WithInfoOfCrisisEffect_MatchInAAI)==0) {
           KOList_WithoutBin <- c(KOList_WithoutBin,KOList[NumKO])
           KOTable_WithoutBin <- rbind(KOTable_WithoutBin,c(KOList[NumKO],NA,NA,"NoBin","NoBin","NoBin"))
         # If no bin in a KO => add in the KOList_WithoutBin and KOTable_WithoutBin, and restart for the next KO
         } else {

  
           
           
  ### 4 - Elimination of KOs with only 1 effect among the bins with this KO -------------------------
      if (length(unique(Code_CrisisEffect_PerBin_WithInfoOfCrisisEffect_MatchInAAI$CodeEffect))>1){ # If TRUE => more than just 1 effect => KO kept for analysis
          KOList_Filtered <- c(KOList_Filtered,KOList[NumKO])
            } else {
              if (unique(Code_CrisisEffect_PerBin_WithInfoOfCrisisEffect_MatchInAAI$CodeEffect)=="2") { # If TRUE => only no effect => KO kept for analysis
                KOList_Filtered <- c(KOList_Filtered,KOList[NumKO])
                  } else {
                    KOList_OneEffect_BeforeAfter <- c(KOList_OneEffect_BeforeAfter,KOList[NumKO])
                    KOTable_OneEffect_BeforeAfter <- rbind(KOTable_OneEffect_BeforeAfter,c(KOList[NumKO],NA,NA,sum(Code_CrisisEffect_PerBin_WithInfoOfCrisisEffect_MatchInAAI$CodeEffect=="1"),sum(Code_CrisisEffect_PerBin_WithInfoOfCrisisEffect_MatchInAAI$CodeEffect=="2"),sum(Code_CrisisEffect_PerBin_WithInfoOfCrisisEffect_MatchInAAI$CodeEffect=="3")))
                    }
                 }
  } 
  } 
} 
   
length(KOList)
length(KOList_Filtered) 

KOTable_OneEffect_BeforeAfter
KOTable_WithoutBin












## Analysis using the filtered KO list "KOList_Filtered"  -----------------------------
AAIAverage <- data.frame(MaxAvg=NA,MaxAvg_Without100=NA)
Table_ApprovedKO <-  data.frame(KO_ID=NA,Nb_EffectBefore=NA,Nb_WithoutEffect=NA,Nb_EffectAfter=NA)

# Search for KOs from KOList_Filtered in Results_AllBins_KO dataframe to have the corresponding bins for each KO
  for (NumKO in 1:length(KOList_Filtered)) {
  BinList_1KO <- row.names(Results_AllBins_KO_Simple_tr)[which(Results_AllBins_KO_Simple_tr[,which(colnames(Results_AllBins_KO_Simple_tr)==KOList_Filtered[NumKO])] > 0)]
  
  
  
  
    ### 1 - Coding of the crisis effect (For each bin listed in BinList_1KO) ----------------
      Code_CrisisEffect_PerBin <- data.frame(BinID=NA,CodeEffect=NA) # Empty dataframe
  
        # Code :
        # 1 : Before crisis (first column of CrisisImpact_PerMAG)
        # 2 : Without any effect (second column of CrisisImpact_PerMAG)
        # 3 : After crisis (third column of CrisisImpact_PerMAG)
  
       for (Num in 1:length(BinList_1KO)) {
          NumBin_Code_CrisisEffect <- c(BinList_1KO[Num],which(CrisisImpact_PerMAG[which(row.names(CrisisImpact_PerMAG)==BinList_1KO[Num]),]==1))
          Code_CrisisEffect_PerBin[Num,] <- NumBin_Code_CrisisEffect
         }
       # Code_CrisisEffect_PerBin # For 1 KO => Bin ID (column 1) and crisis effect (column 2)
  
 
  

     ### 2 -  Elimination of bins with no match in CrisisImpact_PerMAG ---------------------------
      # NB : bins with no match in CrisisImpact_PerMAG have their ID in the second column of Code_CrisisEffect_PerBin => "BinXXX"
        if (sum(str_detect(Code_CrisisEffect_PerBin[,2],"Bin"))>=1) {
            Code_CrisisEffect_PerBin_WithInfoOfCrisisEffect <- Code_CrisisEffect_PerBin[-which(str_detect(Code_CrisisEffect_PerBin[,2],"Bin")),] # If no match => elimination
          } else {
            Code_CrisisEffect_PerBin_WithInfoOfCrisisEffect <- Code_CrisisEffect_PerBin 
         }
      
      
      
     ### 3 - Elimination of bins with no match in the AAI matrix ----------------------
        BinList_AAIMatrix <- c(colnames(AAIMatrix),row.names(AAIMatrix))
        BinList_WithMatch_AAIMatrix <- c(NULL)
      
        for (i in 1:length(Code_CrisisEffect_PerBin_WithInfoOfCrisisEffect$BinID)) {
          BinList_WithMatch_AAIMatrix[i] <- sum(BinList_AAIMatrix==Code_CrisisEffect_PerBin_WithInfoOfCrisisEffect$BinID[i])>0 # If >0 => bin with match in AAI Matrix => TRUE
          }
        Code_CrisisEffect_PerBin_WithInfoOfCrisisEffect_MatchInAAI <- Code_CrisisEffect_PerBin_WithInfoOfCrisisEffect[which(BinList_WithMatch_AAIMatrix),]
      
  
        Table_ApprovedKO <- rbind(Table_ApprovedKO,c(KOList_Filtered[NumKO],sum(Code_CrisisEffect_PerBin_WithInfoOfCrisisEffect_MatchInAAI$CodeEffect=="1"),sum(Code_CrisisEffect_PerBin_WithInfoOfCrisisEffect_MatchInAAI$CodeEffect=="2"),sum(Code_CrisisEffect_PerBin_WithInfoOfCrisisEffect_MatchInAAI$CodeEffect=="3")))
  
  
     ### 4 - For each KOs, search for AAI between bins detected before the crisis and bins detected after the crisis ---------------------
        
       ### Separation of the crisis effects on bins in 3 vectors => before, after and without effect
          Bins_BeforeCrisis <- Code_CrisisEffect_PerBin_WithInfoOfCrisisEffect_MatchInAAI[which(Code_CrisisEffect_PerBin_WithInfoOfCrisisEffect_MatchInAAI$CodeEffect==1),1]
          Bins_NoEffectCrisis <- Code_CrisisEffect_PerBin_WithInfoOfCrisisEffect_MatchInAAI[which(Code_CrisisEffect_PerBin_WithInfoOfCrisisEffect_MatchInAAI$CodeEffect==2),1]
          Bins_AfterCrisis <- Code_CrisisEffect_PerBin_WithInfoOfCrisisEffect_MatchInAAI[which(Code_CrisisEffect_PerBin_WithInfoOfCrisisEffect_MatchInAAI$CodeEffect==3),1]
          
       ### Clustering of : BinID before and without effect of the crisis (columns of DataframeAAI) and BinID after and without effect of the crisis (rows of DataframeAAI)
          Bins_BeforeCrisisANDWithoutEffect <- c(Bins_BeforeCrisis,Bins_NoEffectCrisis)
          Bins_AfterCrisisANDWithoutEffect <- c(Bins_AfterCrisis,Bins_NoEffectCrisis)
       
       ### Creation of a dataframe indicating AAI scores between bins detected before the crisis (and bins with no effect) and bins detected after the crisis (and bins with no effect)
          DataFrameAAI <- data.frame(matrix(ncol=length(Bins_BeforeCrisisANDWithoutEffect),nrow=length(Bins_AfterCrisisANDWithoutEffect)))
          ColumnNames <- c(NULL)
          RowNames <- c(NULL)
            for (NumRow in 1:length(Bins_AfterCrisisANDWithoutEffect)) {
                PerBinAfterCrisis <- c(NULL)
                for (NumColumn in 1:length(Bins_BeforeCrisisANDWithoutEffect)) {
                  PerBinAfterCrisis[NumColumn] <- AAIMatrix[which(row.names(AAIMatrix)==Bins_AfterCrisisANDWithoutEffect[NumRow]),which(colnames(AAIMatrix)==Bins_BeforeCrisisANDWithoutEffect[NumColumn])]
                  ColumnNames[NumColumn] <- Bins_BeforeCrisisANDWithoutEffect[NumColumn]
                  }
             DataFrameAAI[NumRow,] <- PerBinAfterCrisis
             RowNames[NumRow] <- Bins_AfterCrisisANDWithoutEffect[NumRow]
             }

      ### Adding Bin ID as column and row names
          colnames(DataFrameAAI) <- ColumnNames
          row.names(DataFrameAAI) <- RowNames
  
      ### Search of the AAI max value in each row of the DataframeAAI (for the KO considered in the loop)
          DataFrameAAI_MaxPerRow <- apply(DataFrameAAI,1,max)
          
      ### Calculation of the average max AAI score for each KO (ignoring 100% identity for bins with no effect of the crisis)
          PerKO_AvgAAI <- c(mean(DataFrameAAI_MaxPerRow),mean(DataFrameAAI_MaxPerRow[-which(DataFrameAAI_MaxPerRow==100)]))
          AAIAverage[NumKO,] <- PerKO_AvgAAI
          row.names(AAIAverage)[NumKO] <- KOList_Filtered[NumKO]
  }


    






    ### 5 - Compilation of the results to produce the final result table ---------------------
      Table_ApprovedKO <- Table_ApprovedKO[-1,]
      row.names(Table_ApprovedKO) <- Table_ApprovedKO[,1]
      Table_ApprovedKO <- Table_ApprovedKO[,-1]
      sum(row.names(Table_ApprovedKO) == row.names(AAIAverage)) == nrow(AAIAverage)

      FinalResults_AAIAverage_Temp <- cbind(AAIAverage,Table_ApprovedKO)

      KOTable_OneEffect_BeforeAfter <- KOTable_OneEffect_BeforeAfter[-1,]
      row.names(KOTable_OneEffect_BeforeAfter) <- KOTable_OneEffect_BeforeAfter[,1]
      KOTable_OneEffect_BeforeAfter <- KOTable_OneEffect_BeforeAfter[,-1]

      KOTable_WithoutBin <- KOTable_WithoutBin[-1,]
      row.names(KOTable_WithoutBin) <- KOTable_WithoutBin[,1]
      KOTable_WithoutBin <- KOTable_WithoutBin[,-1]

      FinalResults_AAIAverage <- rbind(FinalResults_AAIAverage_Temp,KOTable_OneEffect_BeforeAfter,KOTable_WithoutBin)

      write.csv(FinalResults_AAIAverage,file="ResultatsFinaux_MoyenneAAI_11Mars2025_Liste_KO_Adrien_859KO_MMDE_16Juin2026.csv",quote=FALSE)








