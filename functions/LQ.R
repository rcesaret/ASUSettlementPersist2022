#### LQ.R
#### Rudolf Cesaretti, 7/4/2022

#### "LQ" = Location Quotient
#### 
#### 
#### 
#### 

require(tidyverse)

###############################################################
########################       LQ       #######################
###############################################################

LQ <- function(df, e, E){
  Denom <- df %>% summarize(Denom = sum(!!sym(e), na.rm=T)/sum(!!sym(E), na.rm=T)) %>% pull(Denom)
  LQ <- df %>% mutate(LQ = (!!sym(e)/!!sym(E))/Denom) %>% pull(LQ)
                        
  return(LQ)
}

###############################################################
########################    LQ.knn    #########################
###############################################################

FLQ.knn <- function(df, e, E, m, k){
  #k = 3
  #m <- TrsprtNet_CDmatList[[1]]
  #diag(m) <- 0
  Denom <- df %>% summarize(Denom = sum(!!sym(e), na.rm=T)/sum(!!sym(E), na.rm=T)) %>% pull(Denom)
  v=NA
  for (i in 1:nrow(m)){
    v[i] <- df %>% filter(AggSite %in% rownames(m)[order(m[,i])[1:(k+1)]]) %>% 
      mutate(LQ = (sum(!!sym(e), na.rm=T)/sum(!!sym(E), na.rm=T)/Denom)) %>% pull(LQ)
  }
  #v <- ifelse(Denom == 0, rep(-1,nrow(m)), v)
  return(v)
}

###############################################################
########################    LQ.dist    ########################
###############################################################

FLQ.dist <- function(df, e, E, m, r){
  Denom <- df %>% summarize(Denom = sum(!!sym(e), na.rm=T)/sum(!!sym(E), na.rm=T)) %>% pull(Denom)
  v=NA
  for (i in 1:nrow(m)){
    v[i] <- df %>% filter(AggSite %in% subset(df,m[,i]<r)$AggSite) %>% 
      mutate(LQ = (sum(!!sym(e), na.rm=T)/sum(!!sym(E), na.rm=T)/Denom)) %>% pull(LQ)
  }
  #v <- ifelse(Denom == 0, rep(-1,nrow(m)), v)
  return(v)
}
