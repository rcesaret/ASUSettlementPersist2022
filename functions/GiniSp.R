#### GiniSp.R
#### Rudolf Cesaretti, 7/4/2022

#### "LQ" = Location Quotient
#### 
#### 
#### 
#### 

require(tidyverse)
require(ineq)

###############################################################
########################    GiniSp.knn    #####################
###############################################################

GiniSp.knn <- function(df, var, m, k){
  
  Global <- ineq::Gini(as.numeric(unlist(df[,var])))
  v=NA
  for (i in 1:nrow(m)){
    v[i] <- df %>% filter(AggSite %in% rownames(m)[order(m[,i])[1:(k+1)]]) %>% 
      summarize(spGini = ineq::Gini(!!sym(var))/Global) %>% pull(spGini)
  }
  #v <- ifelse(Global == 0, rep(-1,nrow(m)), v)
  return(v)
}

###############################################################
######################    GiniSp.dist    ######################
###############################################################


GiniSp.dist <- function(df, var, m, r){

  Global <- ineq::Gini(as.numeric(unlist(df[,var])))
  v=NA
  for (i in 1:nrow(m)){
    v[i] <- df %>% filter(AggSite %in% subset(df,m[,i]<r)$AggSite) %>% 
      summarize(spGini = ineq::Gini(!!sym(var))/Global) %>% pull(spGini)
  }
  #v <- ifelse(Global == 0, rep(-1,nrow(m)), v)
  return(v)
}



