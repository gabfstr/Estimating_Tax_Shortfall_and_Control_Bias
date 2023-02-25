rm(list=ls())
graphics.off()

library(fastDummies)
library(MatchIt)
library(glmnet)
set.seed(1)

#Import et visualisation du jeu de donnees de depart
base_mag<-read.csv("~/MAG/Data_MSA_MAG/base_mag.txt", sep='\t', header=TRUE, stringsAsFactors = TRUE)


#On veut se debarasser des variables avec :
# - des NA (sauf les variables cibles)
# - Pas d'explication
# - 1 seule valeur prise

withNA<-c()
only1va<-c()
for(ele in names(base_mag)){
  if(sum(is.na(base_mag[ele]))>0){
    withNA<-append(withNA,ele)
  }
  if(nrow(unique(base_mag[ele]))==1){
    only1va<-append(only1va,ele)
  }
}

#on garde la variable cible et son log
withNA<-withNA[1:21]


#A dichotomiser / deja fait (donc a enlever)
aDichotomiser<-c("catJur_cl", "num_ma_reg")
dejaDicho<-c("ccdrpdta7cl","segment", "nbsal","nbsal4cl","saison3cl", "nbsaison", "type_entrep3cl")

#A enlever
toDrop<-c("ccdrpdta","regionhc4",
          "dept2", "nmsa", "EETA","cdpteta", "region", "ma_region", "regionhc2",
          "saison4cl", "saison2cl",
          "ccdrpdta71",  "ccdrpdta72",  "ccdrpdta73",  "ccdrpdta74",  
          "ccdrpdta75",  "ccdrpdta76",  "ccdrpdta77",
          "cr1","cr2","cr3","cr4","cr5","cr6","cr7","cr8","cr9","cr10",
          "cr11","cr12","cr13","cr14","cr15",
          "ct1","ct2","ct3","ct4","ct5","ct6","ct7","ct8","ct9","ct10",
          "ct11","ct12","ct13","ct14","ct15", "ct16",
          "cr_reg1", "cr_reg3", "cr_reg4", "cr_reg5", "cr_reg8", 
          "cr_reg9", "cr_reg10", "cr_reg11", "cr_reg12", "cr_reg13", 
          "cr_reg14", "cr_reg15", "cr_reg16",
          "ccdrpdta15cl", "ccdrpdta7cl2", "ccdrpdta13cl", "saison4cl", "saison2cl")


#pas d'explications
noExplanation<-c("zone1", "zone2", "zone21", "zone22", "zone23", "zone23", "zone24", 
                 "zone25", "zone3", "zone4", "z1","z21", "z22")

#certaines var sont juste inutiles :
noUse<-c("trim1","trim2","trim3","trim4", "redres_legal_2014")



drop<-c(withNA, only1va, toDrop, noExplanation, noUse, dejaDicho,"id")
base_mag2<-base_mag[!(names(base_mag) %in% drop)]

#Discretisation cfcatjur : 
base_mag2$catJur_cl<-trunc(base_mag2$cfcatjur * (10**(-3)))
base_mag2<-base_mag2[!(names(base_mag2) == "cfcatjur")]

#Dichotomisations
base_mag_clean<-dummy_cols(base_mag2, select_columns = aDichotomiser)
base_mag_clean<-base_mag_clean[!(names(base_mag_clean) %in% aDichotomiser)]

#On supprime un outlier non contrôlé (10x plus que 2eme sur presque toutes les valeurs numeriques
base_mag_clean<-base_mag_clean[-160388,]
#Correction : 2 obs où lmontant = NA
base_mag_clean<-base_mag_clean[!(is.na(base_mag_clean$lmontant) & base_mag_clean$redres_2014==1),]
rownames(base_mag_clean) <- NULL

#BASE DES INDIVIDUS CONTROLES (5,466 obs)
base_controle<-subset(base_mag_clean,base_mag_clean$ctrlex_2014==1)
base_controle<-base_controle[!(names(base_controle)=="ctrlex_2014")]

#BASE DES INDIVIDUS CONTROLES ET FRAUDEURS (1,612 obs)
base_fraude<-subset(base_controle,!(is.na(base_controle$Mt_redress_tot)))
base_fraude<-base_fraude[!(names(base_fraude) %in% c("ctrlex_2014", "redres_2014"))]
rownames(base_fraude)<-NULL

############################### Linear Regression Base ########################

#On supprime MTRedressTot (la regression se fait sur lmontant, le log)
base_fraude<- base_fraude[!(names(base_fraude)=="Mt_redress_tot")]

#### Gestion des singularités 

#Singularities : 
aSuppr<-c("catJur_cl_4", "cotitotsq", "lcoti_tot","lcotitotsq", 
          "TXCOUV_REG1", "TXCOUV_REG2", "TXCOUV_REG3", "TXCOUV_REG4", 
          "TXCOUV_REG5", "TXCOUV_REG6", "TXCOUV_REG7", "TXCOUV_REG8", 
          "TXCOUV_REG9", "TXCOUV_REG10", 
          "LCOTITOT_REG1", "LCOTITOT_REG2", "LCOTITOT_REG3", "LCOTITOT_REG4", 
          "LCOTITOT_REG5", "LCOTITOT_REG6", "LCOTITOT_REG7", "LCOTITOT_REG8", 
          "LCOTITOT_REG9", "LCOTITOT_REG10",
          "LCOTITOT_REG_sq1", "LCOTITOT_REG_sq2", "LCOTITOT_REG_sq3", "LCOTITOT_REG_sq4", 
          "LCOTITOT_REG_sq5", "LCOTITOT_REG_sq6", "LCOTITOT_REG_sq7", "LCOTITOT_REG_sq8", 
          "LCOTITOT_REG_sq9", "LCOTITOT_REG_sq10")
defColinearities<-c("heu_tec","heucdd","mrmcdd","coti_deb","salhom_nbsal", 
                    "heutec_nbheures","joucdd", "cotideb_cotitot", "nbcdd", 
                    "embcdd", "embcdi_nbtot", "entrep_bur", "nbtot_nbsal", 
                    "COTITOT_REG5", "COTITOT_REG7")


discretCompl<-c( "ETP3112", "segment2", "nbsal42", "saison3cl_2", "cd3",
                 "catJur_cl_5", "num_ma_reg_4", "num_ma_reg_10")




#### Delete singularities
singularities <- c(aSuppr, defColinearities, discretCompl)
base_fraude<-base_fraude[!(names(base_fraude) %in% singularities)]

########################### END OF LM DATA BASE PART


#BASE POUR LE LOGIT
base_logit<-base_mag_clean[!(names(base_mag_clean) %in% c("Mt_redress_tot", "lmontant"))]


####### LOGIT
#### Gestion des singularités 

#Singularities : 
aSuppr<-c("cotitotsq", "lcoti_tot","coti_tot", 
          "TXCOUV_REG1", "TXCOUV_REG2", "TXCOUV_REG3", "TXCOUV_REG4", 
          "TXCOUV_REG5", "TXCOUV_REG6", "TXCOUV_REG7", "TXCOUV_REG8", 
          "TXCOUV_REG9", "TXCOUV_REG10", 
          "LCOTITOT_REG1", "LCOTITOT_REG2", "LCOTITOT_REG3", "LCOTITOT_REG4", 
          "LCOTITOT_REG5", "LCOTITOT_REG6", "LCOTITOT_REG7", "LCOTITOT_REG8", 
          "LCOTITOT_REG9", "LCOTITOT_REG10",
          "COTITOT_REG1", "COTITOT_REG2", "COTITOT_REG3", "COTITOT_REG4", 
          "COTITOT_REG5", "COTITOT_REG6", "COTITOT_REG7", "COTITOT_REG8", 
          "COTITOT_REG9", "COTITOT_REG10")
defColinearities<-c("heu_bur","heucdi","mrmcdi","coti_deb","salhom_nbsal", 
                    "heutec_nbheures","joucdi", "cotideb_cotitot", "nbcdi", 
                    "embcdi", "embcdi_nbtot", "entrep_bur", "nbtot_nbsal", 
                    "LCOTITOT_REG_sq1", "LCOTITOT_REG_sq9")


#DISCRET EN TROP :
discretCompl<-c("segment1", "nbsal43", "saison3cl_3", "cd7",
                "catJur_cl_3", "num_ma_reg_6", "num_ma_reg_9", "ETP0101")


####### Delete singularities
singularities <- c(aSuppr, defColinearities, discretCompl)
base_logit<-base_logit[!(names(base_logit) %in% singularities)]



###### Lasso feature selection results
lassoKept<-c("ancentr","pcadre","pbureau","age_ec", "nombre_mensuel_moyen_CDI","heu_tec",
             "nbviecdi", "cotiemp_cotitot", "entrep_tec", "zeroturn","nbsaison_nbtot",
             "EFF0101", "entrep_controleur", "segment2", "nbsal41", "saison3cl_2","cd2", 
             "cd3", "cd6", "meancoti2", "sumcotitot", "LCOTITOT_REG_sq2","LCOTITOT_REG_sq3", 
             "LCOTITOT_REG_sq6", "LCOTITOT_REG_sq7", "lcotitotsq", "catJur_cl_1", 
             "catJur_cl_5", "catJur_cl_8" ,"num_ma_reg_1","num_ma_reg_7" ,"num_ma_reg_8")
#32 features



########## POST-LASSO

###### Matching

#inverting control (treated has to be uncontrolled)
base_test<-base_logit[names(base_logit) %in% c(lassoKept,"ctrlex_2014", "redres_2014")]
fraud<-which(base_test$redres_2014 ==1)

base_test["invert_ctrlex_2014"]= 1 - base_test["ctrlex_2014"]
base_test<-base_test[!(names(base_test) %in% c("ctrlex_2014", "redres_2014"))]

# Initialisation des bases finales
base_finale<-base_mag_clean
base_finale$id<-(1:nrow(base_finale))
Resume<-data.frame(ctrlex=base_finale$ctrlex_2014, redress = base_finale$Mt_redress_tot)

## Preparation de l'ajout de la Regression lineaire conditionnellement a la fraude potentielle
# reprise du modele de regression sur lmontant
disturbing<-c(63, 876, 1336, 1562)
base_fraudeb <- base_fraude[-disturbing,]
modele1 <- lm(lmontant ~ mensu + sal_hom + cotiemp_cotitot + segment3 + entrep_controleur + num_ma_reg_5 + num_ma_reg_8 + embcdd_nbtot + num_ma_reg_2 , data = base_fraudeb)


# on testera un modele k-nn pour chacun des k suivants : 
ratios <-c(3,5,7)

for(nbRatio in ratios){
  cat(paste("Matching en cours pour ratio ",nbRatio,"nn",sep=""), end="\n")
  cat("Cela peut prendre jusqu'à plusieurs minutes", end="\n")
  start_time<-Sys.time()
  response<- matchit( invert_ctrlex_2014 ~., data = base_test, method = "nearest", distance = "logit", 
                      replace = TRUE, ratio = nbRatio)
  summary(response)
  cat("Matching terminé",end="\n")
  cat(paste("Temps d'exécution :",Sys.time()-start_time),end="\n")
  
  
  ######### Matching exploitation
  data2<-get_matches(response)
  data2$id<-as.integer(data2$id)
  data2$subclass<-as.integer(data2$subclass)
  nbMaj=ifelse(nbRatio%%2==1, as.integer(nbRatio/2),(nbRatio/2) -1)
  fraud<-which(base_test$redres_2014 ==1)
  
  #### For Knn on redres_2014
  
  #line to get info about known fraud
  data2<-merge(data2,data.frame(id =base_finale$id, redres_2014 = base_finale$redres_2014, Mt_redress_tot=base_finale$Mt_redress_tot))
  #Counting occurences of matching with fraudeurs
  data_fraud<-aggregate(data2$redres_2014, by=list(data2$subclass), FUN=sum)
  
  ### mean
  data_for_mean<-subset(data2,data2$invert_ctrlex_2014==0)
  data_mean_redress<-aggregate(ifelse( is.na(data_for_mean$Mt_redress_tot),0,data_for_mean$Mt_redress_tot ), by = list(data_for_mean$subclass), FUN=mean)
  ###############
  data2<-merge(data2, data.frame(subclass=data_fraud$Group.1, maj_fraud=data_fraud$x))
  data2<-merge(data2, data.frame(subclass=data_mean_redress$Group.1, max_redress = data_mean_redress$x))
  
  
  data2$fraudPrediction<-ifelse(data2$maj_fraud > nbMaj, data2$max_redress, 0)
  data3<-subset(data2,data2$invert_ctrlex_2014==1)
  idNotThere<-c(subset(base_finale, !(base_finale$id  %in% data3$id))$id)
  
  str_fraud<-paste("fraud_",nbRatio,"nn", sep="")
  str_fraudPrediction<-paste("fraudPrediction_",nbRatio,"nn", sep="")
  str_linear_fraud<-paste("linear_prediction_",nbRatio,"nn", sep="")
  
  base_finale[str_fraud]<-ifelse(base_finale$id %in% idNotThere, 0, 1*(data3$fraudPrediction>0))
  base_finale[str_fraudPrediction]<-ifelse(base_finale$id %in% idNotThere, 0, data3$fraudPrediction)
  base_finale[str_linear_fraud]<-ifelse(base_finale[,str_fraud]==1,exp(predict.lm(modele1,base_finale)),0)
  
  
  Resume[str_fraud]<-base_finale[,str_fraud]
  Resume[str_fraudPrediction]<-base_finale[,str_fraudPrediction]
  Resume[str_linear_fraud]<-base_finale[,str_linear_fraud]
}


# Résumé des predictions
View(Resume)



################ Estimation du nombre de fraudeurs :

######## 3-nn :
potentialFraud_2014_3nn <- sum(Resume$fraud_3nn)
potentialFraud_2014_3nn
# 40,261 


######## 5-nn :
potentialFraud_2014_5nn <- sum(Resume$fraud_5nn)
potentialFraud_2014_5nn
# 30,980 


######## 7-nn :
potentialFraud_2014_7nn <- sum(Resume$fraud_7nn)
potentialFraud_2014_7nn
# 23,246 




##################### Estimation du manque a gagner :

#Montant de redressement effectif : 
redress_2014 <- sum(subset(Resume,!is.na(Resume$redress))$redress)
#11,525,959 euros


############# k-nn puis Regression :

######## 3-nn :
mag_lin_2014_3nn <- sum(Resume$linear_prediction_3nn)
mag_lin_2014_3nn
# 78,135,182


######## 5-nn :
mag_lin_2014_5nn <- sum(Resume$linear_prediction_5nn)
mag_lin_2014_5nn
# 36,245,498


######## 7-nn :
mag_lin_2014_7nn <- sum(Resume$linear_prediction_7nn)
mag_lin_2014_7nn
# 54,933,906



############# k-nn puis k-nn :

######## 3-nn :
mag_2014_3nn <- sum(Resume$fraudPrediction_3nn)
mag_2014_3nn
# 156,876,663 


######## 5-nn :
mag_2014_5nn <- sum(Resume$fraudPrediction_5nn)
mag_2014_5nn
# 108,155,659 


######## 7-nn :
mag_2014_7nn <- sum(Resume$fraudPrediction_7nn)
mag_2014_7nn
# 76,853,934

