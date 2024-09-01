################################################################################
############# Structural Equations: Finql Project ##############################
################################################################################
###################### By Manuel Esteban Arias #################################
################################################################################
################################################################################

# Importing libraries
library(dplyr)
library(shape)
library(amap)
library(diagram)
library(turner)
library(tester)
library(plspm)
library(mnormt)
library(pbivnorm)
library(lavaan)
library(corrplot)

# Importing the data
df0 <- read.table("C:/Users/manue/Downloads/mobil_init.txt",header=TRUE)
#Non Standarized data
df<-df0[,-1]
# Standarized data
df_stand <- scale(df[,-1])

################################################################################
#################  Question 1 :  ###############################################
################################################################################
#  Explorer les données avec les statistiques descriptives et fournir un état 
#  des lieux de celles-ci.


cor_df <- cor(df[,-1])
corrplot(cor_df)

col1 <- colorRampPalette(c("#7F0000", "red", "#FF7F00", "yellow", "white",
                                    "cyan", "#007FFF", "blue", "#00007F"))
col2 <- colorRampPalette(c("#67001F", "#B2182B", "#D6604D", "#F4A582",
                           "#FDDBC7", "#FFFFFF", "#D1E5F0", "#92C5DE",
                            "#4393C3", "#2166AC", "#053061"))
col3 <- colorRampPalette(c("red", "white", "blue")) 
col4 <- colorRampPalette(c("#7F0000", "red", "#FF7F00", "yellow", "#7FFF7F",
                           "cyan", "#007FFF", "blue", "#00007F"))
whiteblack <- c("white", "black")
                                                                                                            
corrplot(cor_df, method = "color", addrect = 2, col = col1(100))
                                                                                                            
ggpairs(df[,-1])

################################################################################
#################  Question 3 :  ###############################################
################################################################################
#  Vérifier l’unidimensionnalité de chacun des blocs. 
# Fonction pour le test d'unidimensionnalité des bloc
unidim_new=function(data,name_blocks,alpha)
                    {n=nrow(data)
                    res_unidim <- unidim(data,blocks = name_blocks)
                    seuil <- matrix(0,nrow=nrow(res_unidim))
                    pval <- matrix(0,nrow=nrow(res_unidim))
                    for (k in 1:nrow(res_unidim))
                    {seuil[k]=1+qnorm(1-alpha/2)*sqrt((res_unidim[k,2]-1)/(n-1))
                    pval[k]=1-pnorm((res_unidim[k,6]-1)/sqrt((res_unidim[k,2]-1)/(n-1)),0,1)}
                    out_unidim <- data.frame(res_unidim,seuil,pval)
                    print(out_unidim)
                    }
PerQual = c(0,0,0,0,0)
Perval = c(1,0,0,0,0)
CusSat = c(1,1,0,0,0)
CusExp=c(1,1,1,0,0)
Loyal = c(1,1,1,1,0)
sat_path = rbind(PerQual, Perval, CusSat,CusExp,Loyal)

# plot diagram of path matrix
innerplot(sat_path,box.cex=2,txt.col="blue",lcol="blue",colpos="green",colneg="orange")

# blocks of outer model
sat_blocks = list(1:3, 4:10,11:12,13:15,16:18)

unidim_new(df,sat_blocks,alpha=0.05)


################################################################################
#################  Question 4 :  ###############################################
################################################################################
#  Appliquer la méthode LISREL sur ce modèle à l’aide du package lavaan. 
# Dérouler l’ensemble de la démarche de validation statistique vue en cours. 
# Puis commenter de façon détaillée cette démarche et les résultats obtenus.
########## Model 1
model <- '
# Modèle de mesure
PerQual =~ PERQ1 + PERQ2 + PERQ3+ PERQ4 + PERQ5 + PERQ6 + PERQ7 
PerVal =~ PERV1 + PERV2
CusSat =~ CUSA1 + CUSA2 + CUSA3 
Loyal =~  CUSL1 + CUSL2 + CUSL3
CusExp =~ CUEX1 + CUEX2 + CUEX3 

# Modèle de structure
PerQual ~ CusExp
PerVal ~ CusExp + PerQual 
CusSat ~ PerVal + CusExp + PerQual 
Loyal ~ CusSat
'

FIT <- sem(model, data = df) 
result <- summary(FIT, standardized = TRUE, fit.measures=TRUE)
result

mod_ind <- modindices(fit, minimum.value = 0, sort = TRUE)
p_val <- 1-pchisq(mod_ind$mi,1)
mod_ind <- data.frame(mod_ind,p_val)
mod_ind_signif <- mod_ind %>% filter(p_val < 0.025)
mod_ind_signif
 
########## Model 2

model <- '
# Modèle de mesure
PerQual =~ PERQ1 + PERQ2 + PERQ3+ PERQ4 + PERQ5 + PERQ6 + PERQ7 
PerVal =~ PERV1 + PERV2
CusSat =~ CUSA1 + CUSA2 + CUSA3 
Loyal =~  CUSL1  + CUSL3
CusExp =~ CUEX1 + CUEX2 + CUEX3 

# Modèle de structure
PerQual ~ CusExp
PerVal ~ CusExp + PerQual 
CusSat ~ PerVal + CusExp + PerQual 
Loyal ~ CusSat

'

fit <- sem(model, data = df) 
summary(fit, standardized = TRUE, fit.measures=TRUE)

mod_ind <- modindices(fit, minimum.value = 0, sort = TRUE)
p_val <- 1-pchisq(mod_ind$mi,1)
mod_ind <- data.frame(mod_ind,p_val)
mod_ind_signif <- mod_ind %>% filter(p_val < 0.025)
mod_ind_signif
########## Model 3

model <- '
# Modèle de mesure
PerQual =~ PERQ1 + PERQ2 + PERQ3+ PERQ4 + PERQ5 + PERQ6 + PERQ7 
PerVal =~ PERV1 + PERV2
CusSat =~ CUSA1 + CUSA2 + CUSA3 
Loyal =~  CUSL1  + CUSL3
CusExp =~ CUEX1 + CUEX2 + CUEX3 

# Modèle de structure
PerQual ~ CusExp
PerVal ~ CusExp + PerQual 
CusSat ~ PerVal + CusExp + PerQual 
Loyal ~ CusSat
# Covariance entre les erreurs
PERQ1 ~~ CUEX3

'

fit <- sem(model, data = df) 
summary(fit, standardized = TRUE, fit.measures=TRUE)

mod_ind <- modindices(fit, minimum.value = 0, sort = TRUE)
p_val <- 1-pchisq(mod_ind$mi,1)
mod_ind <- data.frame(mod_ind,p_val)
mod_ind_signif <- mod_ind %>% filter(p_val < 0.025)
mod_ind_signif
########## Model 4
model <- '
# Modèle de mesure
PerQual =~ PERQ1 + PERQ2 + PERQ3+ PERQ4 + PERQ5 + PERQ6 + PERQ7 
PerVal =~ PERV1 + PERV2
CusSat =~ CUSA1 + CUSA2 + CUSA3 
Loyal =~  CUSL1  + CUSL3
CusExp =~ CUEX1 + CUEX2 + CUEX3 

# Modèle de structure
PerQual ~ CusExp
PerVal ~ CusExp + PerQual 
CusSat ~ PerVal + CusExp + PerQual 
Loyal ~ CusSat
# Covariance entre les erreurs
PERQ1 ~~ CUEX3
PERQ3 ~~   PERQ7
# Impose la variance de la fidélité à 0
Loyal ~~ 0*Loyal
'

fit <- sem(model, data = df) 
summary(fit, standardized = TRUE, fit.measures=TRUE)

mod_ind <- modindices(fit, minimum.value = 0, sort = TRUE)
p_val <- 1-pchisq(mod_ind$mi,1)
mod_ind <- data.frame(mod_ind,p_val)
mod_ind_signif <- mod_ind %>% filter(p_val < 0.025)
mod_ind_signif
########## Model 5
model <- '
# Modèle de mesure
PerQual =~ PERQ1 + PERQ2 + PERQ3+ PERQ4 + PERQ5 + PERQ6 + PERQ7 
PerVal =~ PERV1 + PERV2
CusSat =~ CUSA1 + CUSA2 + CUSA3 
Loyal =~  CUSL1  + CUSL3
CusExp =~ CUEX1 + CUEX2 + CUEX3 

# Modèle de structure
PerQual ~ CusExp
PerVal ~ CusExp + PerQual 
CusSat ~ PerVal + CusExp + PerQual 
Loyal ~ CusSat
# Covariance entre les erreurs
PERQ1 ~~ CUEX3
# Impose la variance de la fidélité à 0
Loyal ~~ 0*Loyal
'

fit <- sem(model, data = df) 
summary(fit, standardized = TRUE, fit.measures=TRUE)

mod_ind <- modindices(fit, minimum.value = 0, sort = TRUE)
p_val <- 1-pchisq(mod_ind$mi,1)
mod_ind <- data.frame(mod_ind,p_val)
mod_ind_signif <- mod_ind %>% filter(p_val < 0.025)
mod_ind_signif
########## Model 5
model <- '# Modèle de mesure
PerQual =~ PERQ1 + PERQ2 + PERQ3+ PERQ4 + PERQ5 + PERQ6 + PERQ7 
PerVal =~ PERV1 + PERV2
CusSat =~ CUSA1 + CUSA2 + CUSA3 
Loyal =~  CUSL1  + CUSL3
CusExp =~ CUEX1 + CUEX2 + CUEX3 
# Modèle de structure
PerQual ~ CusExp
PerVal ~ CusExp + PerQual 
CusSat ~ PerVal + CusExp + PerQual 
Loyal ~ CusSat
# Covariance entre les erreurs
PERQ1 ~~ CUEX3
# Impose la variance de la fidélité à 0
Loyal ~~ 0*Loyal
'

fit <- sem(model, data = df) 
summary(fit, standardized = TRUE, fit.measures=TRUE)

###Not convergent, we hold model 5

lavResiduals(fit)


################################################################################
#################  Question 6 :  ###############################################
################################################################################
#  Utiliser le package plspm pour ajuster ce modèle en choisissant. Dérouler 
# l’ensemble de la démarche de validation statistique vue en cours et commenter 
# de façon détaillée les résultats obtenus.
# Fonction pour le test d'unidimensionnalité des bloc
unidim_new=function(data,name_blocks,alpha)
                    {n=nrow(data)
                    res_unidim <- unidim(data,blocks = name_blocks)
                    seuil <- matrix(0,nrow=nrow(res_unidim))
                    pval <- matrix(0,nrow=nrow(res_unidim))
                    for (k in 1:nrow(res_unidim))
                    {seuil[k]=1+qnorm(1-alpha/2)*sqrt((res_unidim[k,2]-1)/(n-1))
                    pval[k]=1-pnorm((res_unidim[k,6]-1)/sqrt((res_unidim[k,2]-1)/(n-1)),0,1)}
                    out_unidim <- data.frame(res_unidim,seuil,pval)
                    print(out_unidim)
                    }

# path matrix
CusExpl =c(0,0,0,0,0)
PerQua = c(1,0,0,0,0)
Perval = c(1,1,0,0,0)
CusSat=  c(1,1,1,0,0)
Loyal =  c(1,1,1,1,0)
sat_path = rbind(CusExpl,PerQua, Perval, CusSat,Loyal)

# plot diagram of path matrix
innerplot(sat_path,box.cex=2,txt.col="blue",lcol="blue",colpos="green",colneg="orange")
##Taking agawy CUSL2 becasues is not correlated to anything
df2=df[,-17]
# blocks of outer model
sat_blocks = list(1:3, 4:10,11:12,13:15,16:17)

unidim_new(df2,sat_blocks,alpha=0.05)
# vector of modes (reflective indicators)
sat_mod = rep("A", 5)

# apply plspm

satpls3 = plspm(df2, sat_path, sat_blocks, modes = sat_mod, scaled = TRUE, boot.val=TRUE, br=100)
summary(satpls3)
# plot diagram of the inner model
innerplot(satpls3,box.cex=2,txt.col="blue",lcol="blue",colpos="green",colneg="orange",cex.txt=2)

# plot loadings
outerplot(satpls3,box.size=0.1,box.cex=1,txt.col="red",lcol="red",colpos="green",colneg="orange",cex.txt=1,what = "loadings")

# plot outer weights
outerplot(satpls3,box.size=0.1,box.cex=1,txt.col="red",lcol="red",colpos="green",colneg="orange",cex.txt=1,what = "weights")
################################################################################
#################  Question 7 :  ###############################################
################################################################################
#Appliquer sur ce même modèle, l’approche RFPC. 


RFPC_pm <- function(data,path,blocks){
  
  # Calcul des nombres d'observations, de variables et de blocs
  n=nrow(data)
  nb_blocks <- nrow(as.matrix(blocks))
  decal <- min(data.frame(blocks[1]))-1
  p=ncol(data)-decal
  nb_blocks_endo=0
  
  
  # Modèle externe : Calcul des premières composantes principales pour chaque bloc et, des corrélations simples entre chaque variable manifeste et de leur variable latente
  redundancy <- matrix(0,nrow=p,ncol=1)
  R2 <- matrix(0,nrow=nb_blocks,ncol=1)
  for (k in 1:nb_blocks)
    {blocks_ind <- data.frame(blocks[k])
    deb <- blocks_ind[1,1]
    fin <- blocks_ind[nrow(blocks_ind),1]
    r <- cor(data[,deb:fin])
    res.eigen <- eigen(r)
    v_kj <- res.eigen$vectors
    z_1 <- scale(data[,deb:fin])%*%v_kj[,1]
    data_z_1 <- data.frame(data[,deb:fin],z_1)
    col <- ncol(data_z_1)
    nb_man <- col-1
    loading <- cor(data_z_1)[1:nb_man,col]
    name <- names(loading)
    weight <- v_kj[,1]
    communality <- loading^2
    if (k == 1) {var_lat <- z_1}
      else {var_lat <- data.frame(var_lat,z_1)}
    block <- matrix(rownames(path)[k],nrow=nb_man,ncol=1)
    if (k == 1) {outer_model <- data.frame(name,block,weight,loading,communality)
                  colnames(z_1) <- rownames(path)[k]
                  z_all <- z_1
                  x_all <- data[,deb:fin]
                  }
      else {outer_model <- rbind(outer_model,data.frame(name,block,weight,loading,communality))
            colnames(z_1) <- rownames(path)[k]
            z_all <- data.frame(z_all,z_1)
            x_all <- data.frame(x_all,data[,deb:fin])
            }
    }
    
  
  # Modèle interne : régression entre les variables endogènes/(endogènes,exogènes)
  
  colnames(var_lat) <- rownames(path)
  for (k in 1:nb_blocks)
      {eq_reg <- paste(rownames(path)[k],"~")
      s=0
      for (m in 1:nb_blocks)
          {if (m != k & path[k,m] == 1) {s=s+1
                if (s > 1) {eq_reg <- paste(eq_reg,"+")
                eq_reg <- paste(eq_reg,rownames(path)[m])
                }                                
      else {eq_reg <- paste(eq_reg,rownames(path)[m])}
            }
  }    
  if (s >= 1) {nb_blocks_endo=nb_blocks_endo+1
  reg <- lm(eq_reg,data=var_lat)
  r2 <- summary(reg)$r.squared
  blocks_ind <- data.frame(blocks[k])
  deb <- blocks_ind[1,1]-decal
  fin <- blocks_ind[nrow(blocks_ind),1]-decal
  redundancy[deb:fin] <- outer_model$communality[deb:fin]*r2
  res <- summary(reg)$coef
  res[1,1:4]=NA
  rownames(res)[1] <- rownames(path)[k]
  if (nb_blocks_endo == 1) {inner_model <- as.matrix(res)}
  else {inner_model <- rbind(inner_model,as.matrix(res))}
  R2[k]=r2
  }
  }
  
  outer_model <- data.frame(outer_model,redundancy)
  titre <- c("OUTER MODEL")
  print(titre)
  print(outer_model)
  
  titre <- c("INNER MODEL")
  print(titre)
  print(inner_model)
  
  # Calcul des corrélations entre variables latentes et variables manifestes
  
  all_man_lat <- data.frame(x_all,z_all)
  p1=p+1
  pq=p+nb_blocks
  crossloadings <- cor(all_man_lat)[1:p,p1:pq]
  titre <- c("CROSSLOADINGS")
  print(titre)
  print(crossloadings)
  
  
  # Calcul des corrélations entre variables latentes
  
  cor_var_lat <- cor(z_all)
  titre <- c("CORRELATIONS BETWEEN LVs")
  print(titre)
  print(cor_var_lat)
  
  
  # Calcul des R2, communautés et redondances
  
  mean_communality <- by(outer_model$communality,outer_model$block,mean)
  mean_redundancy <- by(outer_model$redundancy,outer_model$block,mean)
  summary_inner_model <- cbind(R2[1],mean_communality[1],mean_redundancy[1])
  for (k in 2:nb_blocks)
  {summary_inner_model <- rbind(summary_inner_model,cbind(R2[k],mean_communality[k],mean_redundancy[k]))}
  colnames(summary_inner_model) <- c("R2","Mean communality","Mean Redundancy")
  titre <- c("SUMMARY INNER MODEL")
  print(titre)
  print(summary_inner_model)
  
  
  # Calcul du GoF
  
  GoF <- sqrt(mean(summary_inner_model[,2])*sum(summary_inner_model[,1]/nb_blocks_endo))
  titre <- c("GOODNESS-OF-FIT")
  print(titre)
  print(GoF)
  
  return(list(all_man_lat=all_man_lat,outer_model=outer_model,inner_model=inner_model,summary_inner_model=summary_inner_model,crossloadings=crossloadings,cor_var_lat=cor_var_lat,GoF=GoF))
  
}



res.RFPC <- RFPC_pm(df2,sat_path,sat_blocks)


################################################################################
#################  Question 7 :  ###############################################
################################################################################
#Comparer tout d’abord les pouvoirs explicatifs des trois modèles de mesure et des trois modèles de structure (LISREL, PLS et RFPC). Que remarquez-vous ? 

##################
# Méthode LISREL #
##################

# Ecriture du modèle

model <- '
# Modèle de mesure
PerQual =~ PERQ1 + PERQ2 + PERQ3+ PERQ4 + PERQ5 + PERQ6 + PERQ7 
PerVal =~ PERV1 + PERV2
CusSat =~ CUSA1 + CUSA2 + CUSA3 
Loyal =~  CUSL1  + CUSL3
CusExp =~ CUEX1 + CUEX2 + CUEX3 

# Modèle de structure
PerQual ~ CusExp
PerVal ~ CusExp + PerQual 
CusSat ~ PerVal + CusExp + PerQual 
Loyal ~ CusSat
# Covariance entre les erreurs
PERQ1 ~~ CUEX3
# Impose la variance de la fidélité à 0
Loyal ~~ 0*Loyal
'


# Application de la méthode LISREL
fit <- sem(model, data = df2) 
summary(fit, standardized = TRUE, fit.measures=TRUE)

# Calcul des loadings normalisés pour LISREL
result_teta <- data.frame(parameterEstimates(fit))
man_var <- result_teta$rhs[1:17]
lat_var <- result_teta$lhs[1:17]
lambda_LISREL <- result_teta$est[1:17]
outer_model_LISREL <- data.frame(lat_var,man_var,lambda_LISREL)
cum_lambda <- by(outer_model_LISREL$lambda_LISREL^2,outer_model_LISREL$lat_var,sum)
deno_lambda <- as.matrix(sqrt(by(outer_model_LISREL$lambda_LISREL^2,outer_model_LISREL$lat_var,sum)))
lat_var <- rownames(deno_lambda)
deno_lambda <- data.frame(lat_var,deno_lambda)
outer_model_LISREL <- merge(outer_model_LISREL,deno_lambda,by="lat_var")
lambda_stand_LISREL <- round(outer_model_LISREL$lambda_LISREL/outer_model_LISREL$deno_lambda,4)
outer_model_LISREL <- data.frame(outer_model_LISREL,lambda_stand_LISREL)
outer_model_LISREL$lambda <- round(outer_model_LISREL$lambda_LISREL,4)
#outer_model_LISREL <- outer_model_LISREL[,-4]
#outer_model_LISREL <- outer_model_LISREL[,-5]

# Variables latentes pour LISREL
var_lat_LISREL <- data.frame(predict(fit))
colnames(var_lat_LISREL) <- c("CusExpl","PerQua", "Perval", "CusSat","Loyal")


################
# Approche PLS #
################

# path matrix
CusExpl =c(0,0,0,0,0)
PerQua = c(1,0,0,0,0)
Perval = c(1,1,0,0,0)
CusSat=  c(1,1,1,0,0)
Loyal =  c(1,1,1,1,0)
Russet_path = rbind(CusExpl,PerQua, Perval, CusSat,Loyal)

# plot diagram of path matrix
innerplot(Russet_path,box.cex=2,txt.col="blue",lcol="blue",colpos="green",colneg="orange")

# blocks of outer model
Russet_blocks <- list(1:3, 4:10,11:12,13:15,16:17)

# vector of modes (reflective indicators)
Russet_mod <- rep("A", 5)

# apply plspm
Russet_pls <- plspm(df2, Russet_path, Russet_blocks, modes = Russet_mod, scaled = TRUE)

# Calcul des loadings normalisés pour PLS
outer_model_PLS <- data.frame(Russet_pls$outer_model)
cum_lambda <- by(outer_model_PLS$loading^2,outer_model_PLS$block,sum)
deno_lambda <- as.matrix(sqrt(cum_lambda))
block <- rownames(deno_lambda)
block <- data.frame(block)
deno_lambda <- data.frame(block,deno_lambda)
outer_model_PLS <- merge(outer_model_PLS,deno_lambda,by="block")
lambda_PLS <- round(outer_model_PLS$loading,4)
lambda_stand_PLS <- outer_model_PLS$loading/outer_model_PLS$deno_lambda
outer_model_PLS <- data.frame(outer_model_PLS[,1:2],lambda_PLS,lambda_stand_PLS)
outer_model_PLS$loading_PLS <- round(outer_model_PLS$lambda_PLS,4)
outer_model_PLS$lambda_stand_PLS <- round(outer_model_PLS$lambda_stand_PLS,4)
outer_model_PLS <- outer_model_PLS[,-5]

# Variables latentes pour PLS
var_lat_PLS <- Russet_pls$scores
colnames(var_lat_PLS) <- c("CusExpl","PerQua", "Perval", "CusSat","Loyal")


#################
# Approche RFPC #
#################
data=df2
RFPC_pm <- function(ddata,path,blocks){
  
  # Calcul des nombres d'observations, de variables et de blocs
  
  n=nrow(data)
  nb_blocks <- nrow(as.matrix(blocks))
  decal <- min(data.frame(blocks[1]))-1
  p=ncol(data)-decal
  nb_blocks_endo=0
  
  
  # Modèle externe : Calcul des premières composantes principales pour chaque bloc et, des corrélations simples entre chaque variable manifeste et de leur variable latente
  
  redundancy <- matrix(0,nrow=p,ncol=1)
  R2 <- matrix(0,nrow=nb_blocks,ncol=1)
  for (k in 1:nb_blocks)
  {blocks_ind <- data.frame(blocks[k])
  deb <- blocks_ind[1,1]
  fin <- blocks_ind[nrow(blocks_ind),1]
  r <- cor(data[,deb:fin])
  res.eigen <- eigen(r)
  v_kj <- res.eigen$vectors
  z_1 <- scale(data[,deb:fin])%*%v_kj[,1]
  data_z_1 <- data.frame(data[,deb:fin],z_1)
  col <- ncol(data_z_1)
  nb_man <- col-1
  loading <- cor(data_z_1)[1:nb_man,col]
  name <- names(loading)
  weight <- v_kj[,1]
  communality <- loading^2
  if (k == 1) {var_lat <- z_1}
  else {var_lat <- data.frame(var_lat,z_1)}
  block <- matrix(rownames(path)[k],nrow=nb_man,ncol=1)
  if (k == 1) {outer_model <- data.frame(name,block,weight,loading,communality)
  colnames(z_1) <- rownames(path)[k]
  z_all <- z_1
  x_all <- data[,deb:fin]
  }
  else {outer_model <- rbind(outer_model,data.frame(name,block,weight,loading,communality))
  colnames(z_1) <- rownames(path)[k]
  z_all <- data.frame(z_all,z_1)
  x_all <- data.frame(x_all,data[,deb:fin])
  }
  }
  
  
  # Modèle interne : régression entre les variables endogènes/(endogènes,exogènes)
  
  colnames(var_lat) <- rownames(path)
  for (k in 1:nb_blocks)
  {eq_reg <- paste(rownames(path)[k],"~")
  s=0
  for (m in 1:nb_blocks)
  {if (m != k & path[k,m] == 1) {s=s+1
  if (s > 1) {eq_reg <- paste(eq_reg,"+")
  eq_reg <- paste(eq_reg,rownames(path)[m])
  }                                
  else {eq_reg <- paste(eq_reg,rownames(path)[m])}
  }
  }    
  if (s >= 1) {nb_blocks_endo=nb_blocks_endo+1
  reg <- lm(eq_reg,data=var_lat)
  r2 <- summary(reg)$r.squared
  blocks_ind <- data.frame(blocks[k])
  deb <- blocks_ind[1,1]-decal
  fin <- blocks_ind[nrow(blocks_ind),1]-decal
  redundancy[deb:fin] <- outer_model$communality[deb:fin]*r2
  res <- summary(reg)$coef
  res[1,1:4]=NA
  rownames(res)[1] <- rownames(path)[k]
  if (nb_blocks_endo == 1) {inner_model <- as.matrix(res)}
  else {inner_model <- rbind(inner_model,as.matrix(res))}
  R2[k]=r2
  }
  }
  
  outer_model <- data.frame(outer_model,redundancy)
  outer_model$weight <- round(outer_model$weight,4)
  outer_model$loading <- round(outer_model$loading,4)
  outer_model$communality <- round(outer_model$communality,4)
  outer_model$redundancy <- round(outer_model$redundancy,4)
  titre <- c("OUTER MODEL")
  print(titre)
  print(outer_model)
  
  inner_model[,1] <- round(inner_model[,1],4)
  inner_model[,2] <- round(inner_model[,2],4)
  inner_model[,3] <- round(inner_model[,3],4)
  inner_model[,4] <- round(inner_model[,4],4)
  titre <- c("INNER MODEL")
  print(titre)
  print(inner_model)
  
  # Calcul des corrélations entre variables latentes et variables manifestes
  
  all_man_lat <- data.frame(x_all,z_all)
  p1=p+1
  pq=p+nb_blocks
  crossloadings <- round(cor(all_man_lat)[1:p,p1:pq],4)
  titre <- c("CROSSLOADINGS")
  print(titre)
  print(crossloadings)
  
  
  # Calcul des corrélations entre variables latentes
  
  cor_var_lat <- round(cor(z_all),4)
  titre <- c("CORRELATIONS BETWEEN LVs")
  print(titre)
  print(cor_var_lat)
  
  
  # Calcul des R2, communautés et redondances
  
  mean_communality <- by(outer_model$communality,outer_model$block,mean)
  mean_redundancy <- by(outer_model$redundancy,outer_model$block,mean)
  summary_inner_model <- cbind(R2[1],mean_communality[1],mean_redundancy[1])
  for (k in 2:nb_blocks)
  {summary_inner_model <- rbind(summary_inner_model,cbind(R2[k],mean_communality[k],mean_redundancy[k]))}
  colnames(summary_inner_model) <- c("R2","Mean communality","Mean Redundancy")
  summary_inner_model <- round(summary_inner_model,4)
  titre <- c("SUMMARY INNER MODEL")
  print(titre)
  print(summary_inner_model)
  
  
  # Calcul du GoF
  
  GoF <- round(sqrt(mean(summary_inner_model[,2])*sum(summary_inner_model[,1]/nb_blocks_endo)),4)
  titre <- c("GOODNESS-OF-FIT")
  print(titre)
  print(GoF)
  
  return(list(all_man_lat=all_man_lat,outer_model=outer_model,inner_model=inner_model,summary_inner_model=summary_inner_model,crossloadings=crossloadings,cor_var_lat=cor_var_lat,GoF=GoF))
  
}

Russet_bis <- df2

# path matrix
CusExpl =c(0,0,0,0,0)
PerQua = c(1,0,0,0,0)
Perval = c(1,1,0,0,0)
CusSat=  c(1,1,1,0,0)
Loyal =  c(1,1,1,1,0)
Russet_path = rbind(CusExpl,PerQua, Perval, CusSat,Loyal)

# blocks of outer model
Russet_blocks <- list(1:3, 4:10,11:12,13:15,16:17)


res.RFPC <- RFPC_pm(Russet_bis,Russet_path,Russet_blocks)

# Calcul des loadings normalisés pour RFPC
outer_model_RFPC <- res.RFPC$outer_model
cum_lambda <- by(outer_model_RFPC$loading^2,outer_model_RFPC$block,sum)
deno_lambda <- as.matrix(sqrt(cum_lambda))
block <- rownames(deno_lambda)
block <- data.frame(block)
deno_lambda <- data.frame(block,deno_lambda)
outer_model_RFPC <- merge(outer_model_RFPC,deno_lambda,by="block")
lambda_RFPC <- round(outer_model_RFPC$loading,4)
lambda_stand_RFPC <- outer_model_RFPC$loading/outer_model_RFPC$deno_lambda
outer_model_RFPC <- data.frame(outer_model_RFPC[,1:2],lambda_RFPC,lambda_stand_RFPC)
outer_model_RFPC$lambda_stand_RFPC <- round(outer_model_RFPC$lambda_stand_RFPC,4)

# Variables latentes pour RFPC
var_lat_RFPC <- res.RFPC$all_man_lat[,11:15]
colnames(var_lat_RFPC) <- c("CusExpl","PerQua", "Perval", "CusSat","Loyal")


################################################
# Fusion des loadings pour les trois méthodes  #
################################################

outer_model_3_methodes <- data.frame(outer_model_LISREL,outer_model_PLS[,3:4],outer_model_RFPC[,3:4])
outer_model_3_methodes

##############################################
# Comparaison des corrélations des variables 
#latentes pour les trois méthodes ############
##############################################

#var_lat_PLS <- read.table("C:/Users/C57318/Documents/Cours SEM/Docs étudiants/Comparaison des méthodes en R/var_lat_PLS.txt",header=TRUE)

all_var_lat <- data.frame(var_lat_LISREL,var_lat_PLS,var_lat_RFPC)
round(cor(all_var_lat),4)

#pairs(~ineg_agri_LISREL+ineg_indus_LISREL+inst_pol_LISREL+ineg_agri_PLS+ineg_indus_PLS+inst_pol_PLS+ineg_agri_RFPC+ineg_indus_RFPC+inst_pol_RFPC,data=all_var_lat,main="Simple Scatterplot Matrix")
ggpairs(all_var_lat)

cor_all_var_lat <- cor(all_var_lat)
corrplot(cor_all_var_lat)

col1 <- colorRampPalette(c("#7F0000", "red", "#FF7F00", "yellow", "white",
                                    "cyan", "#007FFF", "blue", "#00007F"))
col2 <- colorRampPalette(c("#67001F", "#B2182B", "#D6604D", "#F4A582",
                           "#FDDBC7", "#FFFFFF", "#D1E5F0", "#92C5DE",
                           "#4393C3", "#2166AC", "#053061"))
col3 <- colorRampPalette(c("red", "white", "blue")) 
col4 <- colorRampPalette(c("#7F0000", "red", "#FF7F00", "yellow", "#7FFF7F",
                            "cyan", "#007FFF", "blue", "#00007F"))
whiteblack <- c("white", "black")
corrplot(cor_all_var_lat, method = "color", addrect = 2, col = col1(100))                                                                                                            

