rm(list=ls())
graphics.off()

library(glmnet)
library(randomForest)
library(caret)
library(smotefamily)
library(fastDummies)

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


#BASE DES INDIVIDUS CONTROLES
base_controle<-subset(base_mag_clean,base_mag_clean$ctrlex_2014==1)
base_controle<-base_controle[!(names(base_controle) =="ctrlex_2014")]

#BASE DES INDIVIDUS NON CONTROLES
base_non_controle<-subset(base_mag_clean,base_mag_clean$ctrlex_2014==0)
base_non_controle<-base_non_controle[!(names(base_non_controle) =="ctrlex_2014")]


#BASE DES INDIVIDUS CONTROLES ET FRAUDEURS
base_fraude<-subset(base_controle,!(is.na(base_controle$Mt_redress_tot)))
base_fraude<-base_fraude[!(names(base_fraude) %in% c("ctrlex_2014", "redres_2014"))]
#remettre lmontant

#petite correction : 2 NA dans lmontant -> delete obs
base_fraude<-base_fraude[!is.na(base_fraude$lmontant),]
base_controle<-base_controle[!(is.na(base_controle$lmontant) & base_controle$redres_2014==1),]


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

######### POST LASSO

lassoKept<-c("lucea", "mensu", "pcadre", "age_moy", "sal_moy", "sal_ec", "nbcdd", 
             "nbviecdd", "coti_etat", "coti_deb", "emb_mrmcdi", "emb_mrmcdd", "nbcdi_nbsal", 
             "nbsaison_nbsal", "salhom_nbsal","salfem_nbsal", "embcdi_nbsal",
             "cotietat_cotitot", "cotideb_cotitot", "entrep_bur", "zeroturn", 
             "nbcdd_nbtot", "embauche", "embauchecdi","EFF0101", "ETP0101", "ETP3112",
             "entrep_controleur","segment1","nbsal42", "nbsal44", "anc_cl", "saison3cl_2",
             "cd2","cd3", "cd4", "cd5", "meancoti", "sumcotitot", "cotitotsq", 
             "TXCOUV_REG1", "TXCOUV_REG2", "TXCOUV_REG7", "TXCOUV_REG10", "COTITOT_REG1",
             "COTITOT_REG2", "COTITOT_REG7", "COTITOT_REG10", "LCOTITOT_REG4", 
             "LCOTITOT_REG_sq2", "LCOTITOT_REG_sq3", "LCOTITOT_REG_sq6", 
             "LCOTITOT_REG_sq8", "LCOTITOT_REG_sq9", "lcotitotsq", "catJur_cl_0",
             "catJur_cl_3", "catJur_cl_4", "catJur_cl_5", "catJur_cl_6", "catJur_cl_7",
             "catJur_cl_8", "num_ma_reg_1", "num_ma_reg_4", "num_ma_reg_5", 
             "num_ma_reg_8", "num_ma_reg_9", "num_ma_reg_10")

#Base controle post lasso 
base_ctrl_RF<-base_controle[names(base_controle) %in% c(lassoKept, "redres_2014")]
ctrl_equilibre <- SMOTE(base_ctrl_RF, base_ctrl_RF$redres_2014,K=5)
base_ctrl_equilibre <- ctrl_equilibre$data
base_ctrl_equilibre<-base_ctrl_equilibre[!(names(base_ctrl_equilibre)=="class")]
#base_non_controle_RF<-base_non_controle[names(base_non_controle) %in% c(lassoKept)]


###### FEATURE SELECTION FOR RF (variable importance)
MostImportant<-c("pcadre","age_moy","sal_moy","sal_ec","nbcdd","coti_etat",
                 "coti_deb","emb_mrmcdi","nbcdi_nbsal","nbsaison_nbsal","salfem_nbsal",
                 "embcdi_nbsal","cotietat_cotitot","cotideb_cotitot","zeroturn",
                 "embauche","embauchecdi","ETP3112","entrep_controleur","nbsal42",
                 "anc_cl","cd2","cd5","meancoti","sumcotitot","LCOTITOT_REG_sq6",
                 "LCOTITOT_REG_sq9","catJur_cl_4","catJur_cl_5","redres_2014")
base_rf<-base_ctrl_equilibre[names(base_ctrl_equilibre)%in% c("redres_2014", MostImportant)]
#### 29 features
base_rf$redres_2014 <-as.factor(base_rf$redres_2014)
levels(base_rf$redres_2014) <- c("non_fraudeur", "fraudeur") 

################### POST NTREE, MTRY, CUTOFF TUNNING : we select :

n <- 600
m <- 5
caliper <- c(0.58,0.42)

# Rf classification Model 5-fold cv
train.control <- trainControl(method = "cv", number = 5,verboseIter=TRUE,classProbs = TRUE, summaryFunction = twoClassSummary)
mod<-train(redres_2014 ~ ., data = base_rf, method = "rf", trControl = train.control, metric ="ROC", maximize = TRUE,
           tuneGrid = data.frame(mtry = c(m)), cutoff=caliper,
           ntree=n)
print(mod)
print(mod$finalModel$confusion)
# ROC : 0.849
# Sensitivity : 0.774
# Specificity : 0.755
# Error I : 0.22
# Error II : 0.23


# Predictions :
base_non_controle$isFraud<-predict(mod, base_non_controle)
nb_fraud<-sum(base_non_controle$isFraud=="fraudeur")
nb_fraud
# 50,130 potential fraudeurs detected



############### predictions regression RF
#Feature selection for linear RF
lasso_lin_kept <- c("mensu","pcadre","sal_moy","sal_hom","coti_etat","embcdi",
                    "cotiemp_cotitot","embcdd_nbtot","embauchecdi","entrep_controleur",
                    "segment1","segment3","nbsal44","cd2","COTITOT_REG6","num_ma_reg_2",
                    "num_ma_reg_5","num_ma_reg_8")
base_lin_RF<-base_fraude[names(base_fraude) %in% c(lasso_lin_kept, "lmontant")]
#18 predictors

#train 5-fold cv 
train.control <- trainControl(method = "cv", number = 5,verboseIter=TRUE)
mod_lin<-train(lmontant ~ ., data = base_lin_RF, method = "rf", trControl = train.control, metric ="Rsquared", maximize = TRUE,
               tuneGrid = data.frame(mtry = c(4)),
               ntree=n)
print(mod_lin)
# MTRY = 4
# RMSE : 1.617
# Rsquared : 0.255
# MAE : 1.278

base_non_controle$potentialFraud<-ifelse(base_non_controle$isFraud=="fraudeur",exp(predict(mod_lin, base_non_controle)),0)




############### predictions regression lineaire
disturbing<-c(63, 876, 1336, 1562)
base_fraudeb <- base_fraude[-disturbing,]
modele1 <- lm(lmontant ~ mensu + sal_hom + cotiemp_cotitot + segment3 + entrep_controleur + num_ma_reg_5 + num_ma_reg_8 + embcdd_nbtot + num_ma_reg_2 , data = base_fraudeb)
summary(modele1)
RMSE<-sqrt(mean(modele1$residuals^2))
RMSE
# RMSE : 1.608

# Linear Prediction
base_non_controle$linearFraud<-ifelse(base_non_controle$isFraud=="fraudeur",exp(predict.lm(modele1,base_non_controle)),0)






##################### Estimation du manque a gagner :

#Montant de redressement effectif : 
redress_2014 <- sum(subset(Resume,!is.na(Resume$redress))$redress)
#11,525,959 euros

mag_2014_RF<-sum(base_non_controle$potentialFraud)
mag_2014_RF
# 80,196,795 of fraud : cohérent, Améliorer ?
# 79,103,441



mag_2014_lin<-sum(base_non_controle$linearFraud)
mag_2014_lin
# 142,090,248
# 140,666,075



#resume les predictions
Resume<-data.frame(isFraud = base_non_controle$isFraud, predictionFraud = base_non_controle$potentialFraud,
                   linear_prediction = base_non_controle$linearFraud)
View(Resume)
#RF has max fraud : 151,000
# Whereas linear prediction has ~28,000,000
# => difference because of largest observations, RF can't extrapolate
