# Script R pour analyser l'effet de la crise sur les MAGs du Dziani
# Date de création : 21 février 2025
# Dernière modification : 10 mars 2025
# Autrice : Perrine Cruaud

# Chargement des librairies --------------
library(stringr)

# Choix du dossier de travail -----------
setwd("~/Documents/Adrien/Dziani/MAGS_AAI_ImpactCrise")

# Importation des fichiers -------------
Resultats_AllBins_KO <- read.csv(file="Resultats_AllBins_KO_WithGeneNames.csv", sep="\t", header=TRUE,
                                 row.names=1)
EffetCrise_ParMAG <- read.csv(file="Table_Effet_Crise_matriceAAI.csv", sep=";", row.names=1) 
MatriceAAI <- read.csv(file="matrice_AAI_All_good_MAGs_cor.csv", sep=";", header=TRUE, row.names=1)
Liste_KO_Adrien_859KO_MMDE_11Mars2025 <- read.table(file="Liste_KO_MMDE.txt")
dim(Liste_KO_Adrien_859KO_MMDE_11Mars2025) # 859 1


# Mise en forme des fichiers ---------------------
# Tranformation de la table résultats KO/Bins pour avoir juste les informations nécessaires
Resultats_AllBins_KO[1:10,1:10]
dim(Resultats_AllBins_KO) #9162 1187
colnames(Resultats_AllBins_KO)[1:10]
Resultats_AllBins_KO_Simple <- Resultats_AllBins_KO[,-c(1,2)] # Elimination des 2 premières colonnes =>  "definition" "knum.1" 
Resultats_AllBins_KO_Simple[1:10,1:10]
Resultats_AllBins_KO_Simple_tr <- as.data.frame(t(Resultats_AllBins_KO_Simple))
Resultats_AllBins_KO_Simple_tr[1:10,1:10]
str(Resultats_AllBins_KO_Simple_tr)

# Exploration des effets crise par MAGs
EffetCrise_ParMAG[1:10,]

# Exploration de la matrice AAI
MatriceAAI[1:10,1:10]



# Mise en place du Script d'analyse - Analyse pour plusieurs KO ---------------------------------

# Récupération du KO d'intérêt
Liste_KO_Temp <- Liste_KO_Adrien_859KO_MMDE_11Mars2025$V1
Liste_KO <- unique(Liste_KO_Temp) # Correction si présence de doublons de KO dans la liste

length(Liste_KO_Temp)
length(unique(Liste_KO_Temp))

## Création d'une liste filtrée (élimination des KO sans bins et sans double effet) -----------------------------

# Différentes problématiques :
# 1/ On extrait les bins qui correspondent aux KO => Certains ne sont pas présents dans la table effet crise => on supprime ces bins
# 2/ Certains ne sont pas présents dans la table matrice AAI => on supprime ces bins
# 3/ Certains KO n'ont plus aucun bin après étape 2 => on élimine ces KOs de l'analyse => nouvelle liste de recherche (mettre en mémoire ceux qu'on élimine)
# 4/ Certains KO n'ont que un ou des bins qui n'ont que des effets après ou que des effets avant => on élimine ces KOs de l'analyse => nouvelle liste de recherche (mettre en mémoire ceux qu'on élimine)

Liste_KO_New <- c(NULL)
Liste_KO_SansBinsCorrespondant <- c(NULL)
Liste_KO_UnSeulEffet_AvantApres <- c(NULL)
Table_KO_UnSeulEffet_AvantApres <- data.frame(KO_ID=NA,MoyenneMax=NA,MoyenneMax_Sans100=NA,Nb_EffetAvant=NA,Nb_SansEffet=NA,Nb_EffetApres=NA)
Table_KO_SansBinsCorrespondant <-  data.frame(KO_ID=NA,MoyenneMax=NA,MoyenneMax_Sans100=NA,Nb_EffetAvant=NA,Nb_SansEffet=NA,Nb_EffetApres=NA)


# Recherche de chacun des KO de la liste dans le fichiers résultats Bins pour avoir les Bins correspondant à chaque KO
for (NumKO in 1:length(Liste_KO)) {
  Liste_Bins_1KO <- row.names(Resultats_AllBins_KO_Simple_tr)[which(Resultats_AllBins_KO_Simple_tr[,which(colnames(Resultats_AllBins_KO_Simple_tr)==Liste_KO[NumKO])] > 0)]
  # Liste_Bins_1KO => Liste de tous les KO pour un bin
  
  # Création d'une table avec 1ère colonne => Nom du Bins pour ce KO et 2ème colonne => Effet Crise
  # Code :
  # 1 : Avant crise
  # 2 : Sans Effet
  # 3 : Apres crise
  
  ### Boucle pour chacun des bins listés dans Liste_Bins_1KO ----------------
  CodageEffetCriseParBin <- data.frame(BinID=NA,CodeEffet=NA) # Création dataframe vide
  
  for (Num in 1:length(Liste_Bins_1KO)) {
    NumBin_CodeEffetCrise <- c(Liste_Bins_1KO[Num],which(EffetCrise_ParMAG[which(row.names(EffetCrise_ParMAG)==Liste_Bins_1KO[Num]),]==1))
    CodageEffetCriseParBin[Num,] <- NumBin_CodeEffetCrise
  }
  
  # CodageEffetCriseParBin # Table correspondance pour un KO => Bins (colonne 1) et Effet Crise (colonne 2)

  
  
  ### Elimination des bins non présents dans le fichier "Table_MAG_EffetCrise.csv" ---------------------
  # Les bins qui n'ont pas de correspondance dans la table d'effet ont leur nom de bin dans la 2ème colonne => "BinXXX"
   if (sum(str_detect(CodageEffetCriseParBin[,2],"Bin"))>=1) {
    CodageEffetCriseParBin_SansBinAbsent <- CodageEffetCriseParBin[-which(str_detect(CodageEffetCriseParBin[,2],"Bin")),] # Si au moins 1 bin absent => enlève de la table
  } else {
    CodageEffetCriseParBin_SansBinAbsent <- CodageEffetCriseParBin # Sinon recopie tel quel
  }
  # CodageEffetCriseParBin_SansBinAbsent # Table correspondance pour un KO => Bins (colonne 1) et Effet Crise (colonne 2) sans les bins absents de la table effet crise
  
  nrow(CodageEffetCriseParBin_SansBinAbsent)==0 # Si TRUE => il faut éliminer ce KO (pas de correspondances de bins)
  
  if (nrow(CodageEffetCriseParBin_SansBinAbsent)==0) {
    Liste_KO_SansBinsCorrespondant <- c(Liste_KO_SansBinsCorrespondant,Liste_KO[NumKO])
    Table_KO_SansBinsCorrespondant <- rbind(Table_KO_SansBinsCorrespondant,c(Liste_KO[NumKO],NA,NA,"AucunBin","AucunBin","AucunBin")) # Si le KO est vide, on l'inscrit dans la liste et on recommence au début, sinon on continue
  } else {

  
    
  ### Elimination des bins non présents dans la matrice AAI ----------------------
    # Si le KO sur lequel on travaille a encore des bins dans CodageEffetCriseParBin_SansBinAbsent, on vérifie si les bins sont bien présents dans la matrice AAI
    # Création d'une liste de tous les bins présents dans la matrice AAi  
    Bins_MatriceAAI <- c(colnames(MatriceAAI),row.names(MatriceAAI))
  
    BinsPresents_MatriceAAI <- c(NULL)
  
  for (i in 1:length(CodageEffetCriseParBin_SansBinAbsent$BinID)) {
    BinsPresents_MatriceAAI[i] <- sum(Bins_MatriceAAI==CodageEffetCriseParBin_SansBinAbsent$BinID[i])>0 # Si c'est sup 0 ça veut dire que le bin est là => TRUE
  }
  
  # BinsAbsents_MatriceAAI => vecteur TRUE/FALSE indiquant pour chaque bin si présents dans AAI ou non => TRUE quand présent
  CodageEffetCriseParBin_SansBinAbsentDuTout <- CodageEffetCriseParBin_SansBinAbsent[which(BinsPresents_MatriceAAI),]
  
  nrow(CodageEffetCriseParBin_SansBinAbsentDuTout)==0 # Si TRUE => il faut éliminer ce KO (pas de correspondances de bins)
  
  if (nrow(CodageEffetCriseParBin_SansBinAbsentDuTout)==0) {
    Liste_KO_SansBinsCorrespondant <- c(Liste_KO_SansBinsCorrespondant,Liste_KO[NumKO])
    Table_KO_SansBinsCorrespondant <- rbind(Table_KO_SansBinsCorrespondant,c(Liste_KO[NumKO],NA,NA,"AucunBin","AucunBin","AucunBin"))
    # Si le KO est vide, on l'inscrit dans la liste et on recommence au début, sinon on continue
  } else {

  
  ### Elimination des KO qui n'ont pas au moins 2 effets différents => ensuite il faudra refaire un liste avec ces KOs en moins -------------------------
    length(unique(CodageEffetCriseParBin_SansBinAbsentDuTout$CodeEffet))>1 # Si c'est vrai => on a plus d'un effet, si c'est faux, on a qu'un seul effet
    length(unique(CodageEffetCriseParBin_SansBinAbsentDuTout$CodeEffet))==1 # si vrai ça veut dire qu'on a qu'une seule condition
  # Que des 1 => on enlève
  # Que des 3 => on enlève
  # Que des 2 => on garde (sans effet)
  unique(CodageEffetCriseParBin_SansBinAbsentDuTout$CodeEffet)=="2"
  # si c'est 2 (sans effet) => on garde, sinon on jette
  # if length ==1 TRUE => si c'est 2 on garde sinon on jette
  
  #=> Si TRUE, met le KO dans une nouvelle liste, si FALSE, ne met pas ce KO dans la nouvelle liste
  
  if (length(unique(CodageEffetCriseParBin_SansBinAbsentDuTout$CodeEffet))>1){
        Liste_KO_New <- c(Liste_KO_New,Liste_KO[NumKO])
          } else {
            if (unique(CodageEffetCriseParBin_SansBinAbsent$CodeEffet)=="2") {
              Liste_KO_New <- c(Liste_KO_New,Liste_KO[NumKO])
                } else {
                  Liste_KO_UnSeulEffet_AvantApres <- c(Liste_KO_UnSeulEffet_AvantApres,Liste_KO[NumKO])
                  Table_KO_UnSeulEffet_AvantApres <- rbind(Table_KO_UnSeulEffet_AvantApres,c(Liste_KO[NumKO],NA,NA,sum(CodageEffetCriseParBin_SansBinAbsentDuTout$CodeEffet=="1"),sum(CodageEffetCriseParBin_SansBinAbsentDuTout$CodeEffet=="2"),sum(CodageEffetCriseParBin_SansBinAbsentDuTout$CodeEffet=="3")))
                  }
          }
  } # Ferme boucle vérification présence dans matrice AAI
  } # Ferme boucle vérification présence dans table effet
} # Ferme boucle vérification pour tous les KOs
   
length(Liste_KO)
length(Liste_KO_New) # Avec la table des 859 KO => 716 KO restants (=> 143 KOs qui vont pas)

Table_KO_UnSeulEffet_AvantApres
Table_KO_SansBinsCorrespondant


## On peut maintenant refaire le même script mais en utilisant cette nouvelle liste --------------------------
MoyenneAAI <- data.frame(MoyenneMax=NA,MoyenneMax_Sans100=NA)
Table_BonsKO <-  data.frame(KO_ID=NA,Nb_EffetAvant=NA,Nb_SansEffet=NA,Nb_EffetApres=NA)


for (NumKO in 1:length(Liste_KO_New)) {
  Liste_Bins_1KO <- row.names(Resultats_AllBins_KO_Simple_tr)[which(Resultats_AllBins_KO_Simple_tr[,which(colnames(Resultats_AllBins_KO_Simple_tr)==Liste_KO_New[NumKO])] > 0)]
  
  # Création d'une table avec 1ère colonne => Nom du Bins pour ce KO et 2ème colonne => Effet Crise
  # Code :
  # 1 : Avant crise
  # 2 : Sans Effet
  # 3 : Apres crise
  
  ### Codage Effet Crise : Boucle pour chacun des bins listés dans Liste_Bins_1KO ---------------------
  CodageEffetCriseParBin <- data.frame(BinID=NA,CodeEffet=NA)
  
  for (Num in 1:length(Liste_Bins_1KO)) {
    NumBin_CodeEffetCrise <- c(Liste_Bins_1KO[Num],which(EffetCrise_ParMAG[which(row.names(EffetCrise_ParMAG)==Liste_Bins_1KO[Num]),]==1))
    CodageEffetCriseParBin[Num,] <- NumBin_CodeEffetCrise
  }
  
    # CodageEffetCriseParBin # Table correspondance pour un KO => Bins (colonne 1) et Effet Crise (colonne 2)
  
  
  
  ### Elimination des bins non présents dans le fichier "Table_MAG_EffetCrise.csv" ------------------------
   # Les bins qui n'ont pas de correspondance dans la table d'effet ont leur nom de bin dans la 2ème colonne => "BinXXX"
  if (sum(str_detect(CodageEffetCriseParBin[,2],"Bin"))>=1) {
    CodageEffetCriseParBin_SansBinAbsent <- CodageEffetCriseParBin[-which(str_detect(CodageEffetCriseParBin[,2],"Bin")),] # Si au moins 1 bin absent => enlève de la table
  } else {
    CodageEffetCriseParBin_SansBinAbsent <- CodageEffetCriseParBin # Sinon recopie tel quel
  }
  # CodageEffetCriseParBin_SansBinAbsent # Table correspondance pour un KO => Bins (colonne 1) et Effet Crise (colonne 2) sans les bins absents de la table effet crise
  
  
  
  ### Elimination des bins non présents dans la matrice AAI -----------------------------
  # On vérifie que les bins sont bien présents pour le KO sur lequel on travaill dans la matrice AAI
  # Création d'une liste de tous les bins présents dans la matrice AAi  
  Bins_MatriceAAI <- c(colnames(MatriceAAI),row.names(MatriceAAI))
  
  BinsPresents_MatriceAAI <- c(NULL)
  
  for (i in 1:length(CodageEffetCriseParBin_SansBinAbsent$BinID)) {
    BinsPresents_MatriceAAI[i] <- sum(Bins_MatriceAAI==CodageEffetCriseParBin_SansBinAbsent$BinID[i])>0 # Si c'est sup 0 ça veut dire que le bin est là => TRUE
  }
  
  # BinsAbsents_MatriceAAI => vecteur TRUE/FALSE indiquant pour chaque bin si présents dans AAI ou non => TRUE quand présent
  CodageEffetCriseParBin_SansBinAbsentDuTout <- CodageEffetCriseParBin_SansBinAbsent[which(BinsPresents_MatriceAAI),]
  Table_BonsKO <- rbind(Table_BonsKO,c(Liste_KO_New[NumKO],sum(CodageEffetCriseParBin_SansBinAbsentDuTout$CodeEffet=="1"),sum(CodageEffetCriseParBin_SansBinAbsentDuTout$CodeEffet=="2"),sum(CodageEffetCriseParBin_SansBinAbsentDuTout$CodeEffet=="3")))
  
  
  ### Separation dans 3 vecteurs des 3 catégories (avant, après et sans effet) ----------------------------
  Bins_AvantCrise <- CodageEffetCriseParBin_SansBinAbsentDuTout[which(CodageEffetCriseParBin_SansBinAbsentDuTout$CodeEffet==1),1]
  Bins_SansEffetCrise <- CodageEffetCriseParBin_SansBinAbsentDuTout[which(CodageEffetCriseParBin_SansBinAbsentDuTout$CodeEffet==2),1]
  Bins_ApresCrise <- CodageEffetCriseParBin_SansBinAbsentDuTout[which(CodageEffetCriseParBin_SansBinAbsentDuTout$CodeEffet==3),1]
  # Poolage des sans effet avec effets avant et effets après crise)
  Bins_AvantCriseEtSansEffet <- c(Bins_AvantCrise,Bins_SansEffetCrise)
  Bins_ApresCriseEtSansEffet <- c(Bins_ApresCrise,Bins_SansEffetCrise)
  
  
  ### Recherche des AAI dans matrice AAI -------------------------
  DataFrameAAI <- data.frame(matrix(ncol=length(Bins_AvantCriseEtSansEffet),nrow=length(Bins_ApresCriseEtSansEffet)))
  NomColonnes <- c(NULL)
  NomLignes <- c(NULL)
  
  for (NumLigne in 1:length(Bins_ApresCriseEtSansEffet)) {
    ParBinApresCrise <- c(NULL)
    for (NumColonne in 1:length(Bins_AvantCriseEtSansEffet)) {
      ParBinApresCrise[NumColonne] <- MatriceAAI[which(row.names(MatriceAAI)==Bins_ApresCriseEtSansEffet[NumLigne]),which(colnames(MatriceAAI)==Bins_AvantCriseEtSansEffet[NumColonne])]
      NomColonnes[NumColonne] <- Bins_AvantCriseEtSansEffet[NumColonne]
    }
    DataFrameAAI[NumLigne,] <- ParBinApresCrise
    NomLignes[NumLigne] <- Bins_ApresCriseEtSansEffet[NumLigne]
  }
  
  
  #Ajout des noms des bins en noms de colonnes et en noms de lignes
  colnames(DataFrameAAI) <- NomColonnes
  row.names(DataFrameAAI) <- NomLignes
  
  #Calcul du max par ligne pour chaque bin avec effet après crise
  # pour le KO considéré
  DataFrameAAI_MaxParLigne <- apply(DataFrameAAI,1,max)
  Pour1KO_MoyenneAAI <- c(mean(DataFrameAAI_MaxParLigne),mean(DataFrameAAI_MaxParLigne[-which(DataFrameAAI_MaxParLigne==100)]))
  
  MoyenneAAI[NumKO,] <- Pour1KO_MoyenneAAI
  row.names(MoyenneAAI)[NumKO] <- Liste_KO_New[NumKO]
  
}


dim(MoyenneAAI) # 716 KO avec la nouvelle liste de 859 KO du 11 mars 2025
length(Liste_KO_New) # 716 KO avec la nouvelle liste de 859 KO du 11 mars 2025
MoyenneAAI
Table_BonsKO <- Table_BonsKO[-1,]
row.names(Table_BonsKO) <- Table_BonsKO[,1]
Table_BonsKO <- Table_BonsKO[,-1]
sum(row.names(Table_BonsKO) == row.names(MoyenneAAI)) == nrow(MoyenneAAI)

ResultatsFinaux_MoyenneAAI_Temp <- cbind(MoyenneAAI,Table_BonsKO)


Table_KO_UnSeulEffet_AvantApres
Table_KO_UnSeulEffet_AvantApres <- Table_KO_UnSeulEffet_AvantApres[-1,]
row.names(Table_KO_UnSeulEffet_AvantApres) <- Table_KO_UnSeulEffet_AvantApres[,1]
Table_KO_UnSeulEffet_AvantApres <- Table_KO_UnSeulEffet_AvantApres[,-1]

Table_KO_SansBinsCorrespondant
Table_KO_SansBinsCorrespondant <- Table_KO_SansBinsCorrespondant[-1,]
row.names(Table_KO_SansBinsCorrespondant) <- Table_KO_SansBinsCorrespondant[,1]
Table_KO_SansBinsCorrespondant <- Table_KO_SansBinsCorrespondant[,-1]

ResultatsFinaux_MoyenneAAI <- rbind(ResultatsFinaux_MoyenneAAI_Temp,Table_KO_UnSeulEffet_AvantApres,Table_KO_SansBinsCorrespondant)

nrow(ResultatsFinaux_MoyenneAAI) # 859
nrow(MoyenneAAI) # 716
nrow(Table_KO_UnSeulEffet_AvantApres) # 79
nrow(Table_KO_SansBinsCorrespondant) #64
nrow(MoyenneAAI)+nrow(Table_KO_UnSeulEffet_AvantApres)+nrow(Table_KO_SansBinsCorrespondant) #859

ResultatsFinaux_MoyenneAAI[1:10,]

write.csv(ResultatsFinaux_MoyenneAAI,file="ResultatsFinaux_MoyenneAAI_11Mars2025_Liste_KO_Adrien_859KO_MMDE_11Mars2025.csv",quote=FALSE)








