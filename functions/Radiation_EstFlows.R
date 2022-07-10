#### Radiation_EstFlows.R
#### Rudolf Cesaretti, 6/31/2022

#### "Radiation_EstFlows" 
#### 
#### 
#### 
#### https://book.archnetworks.net/spatialinteraction

pak <- c("compiler", "msm", "snowfall", "parallel", "doParallel", "rgdal", "sp", 
         "sf", "GISTools", "lwgeom", "tidyverse", "tidyr", "data.table", "zoo")
# Install packages not yet installed
ip <- pak %in% rownames(installed.packages())
if (any(ip == FALSE)) {
  install.packages(pak[!ip])
}
# load packages
invisible(lapply(pak, library, character.only = TRUE))
rm(pak,ip)


library(raster)
library(terra)
library(scales)
library(pracma)
library(modelsummary)
library(cowplot)

# Multicore setup
cl <- makeCluster(detectCores())
registerDoParallel(cl)

###############################################################
##########################  Radiation_EstFlows  #########################
###############################################################

d_mat <- CD.mats[[1]]
pop <- Poly_List[[1]]@data$Population.s2

x = radiation(Poly_List[[7]]@data$Population.s2, CD.mats[[7]], scale = "invariant", commuters="pop")



radiation <- function(pop, 
                      d_mat, 
                      #prob = c("original", "extended"),
                      #alpha = NULL,
                      scale = c("invariant", "variant", "none"),
                      commuters = c("pop", "input", "none"),
                      Ti_input = NULL){
  if (commuters == "pop"){T_i <- pop}
  if (commuters == "input"){T_i <- Ti_input}
  if (commuters == "none"){T_i <- rep(1, length(pop))}
  if (scale == "invariant" | scale == "none"){scalevar <- rep(1, length(pop))}
  if (scale == "variant"){scalevar <- (1 - pop/sum(pop))}
  
  ## create square matrix with rows and columns for every site
  out <-matrix(0, length(pop), length(pop))
  for (i in seq_len(nrow(d_mat))) {
    # start loop on rows = i origins
    for (j in seq_len(length(pop))) {
      # start loop on columns = j destinations
      if (i == j) next()# skip diagonal of matrix
      
      m <- pop[i] # set population value for site i
      n <- pop[j] # set population value for site j
      
      # find radius as distance between sites i and j
      r_ij <- d_mat[i, j]
        
      # find all sites within the distance from i to j centered on i
      sel_circle <- which(d_mat[i, ] <= r_ij)
      
      # remove the site i and j from list
      sel_circle <-sel_circle[-which(sel_circle %in% c(i, j))]
      
      s <- sum(pop[sel_circle]) # sum population within radius
      
      # calculate T_i and output to matrix
      temp <- T_i[i]/scalevar[i] * ((m * n) / ((m + s) * (m + n + s)))
      if (is.na(temp)) temp <- 0
      out[i, j] <- temp
    }
  }
  colnames(out)<- colnames(d_mat)
  rownames(out)<- rownames(d_mat)
  return(out)
}




































#per.short <- foreach(i=1:(nrow(period)-1)) %dopar% (period[i+1,1]-period[i,1])
#cer.per <- foreach(i=1:length(start)) %dopar% (seq((start[i]:end[i]))+(start[i]-1))
#per.len <- foreach(i=1:nrow(period)) %dopar% (seq(period[i,1]:period[i,2])+(period[i,1])-1)
#prob <- foreach(i=1:length(cer.per),.combine='rbind') %:% foreach(j=1:length(per.len),.combine='c') %dopar% (bpcalc(cp=cer.per[[i]], pl=per.len[[j]], a=alpha[i], b=beta[i]))
#prob <- foreach(i=1:length(cer.per),.combine='rbind') %:% foreach(j=1:length(per.len),.combine='c') %dopar% (tnpcalc(cp=cer.per[[i]], pl=per.len[[j]], zq=z[i]))
#prob <- foreach(a=cer.per,.combine='rbind') %:% foreach(b=per.len,.combine='c') %dopar% (ucalc(a,b))

d_mat <- CD.mats[[1]]
pop <- Poly_List[[1]]@data$Population.s2
Ti <- NULL


radiation <- function(pop, 
                      d_mat, 
                      prob = c("original", "extended"),
                      alpha = NULL,
                      scale = c("invariant", "variant", "none"),
                      commuters = c("pop", "input"),
                      Ti = NULL){
  
  ## create square matrix with rows and columns for every site
  out <- matrix(0, length(pop), length(pop))
  
  #prob <- foreach(i=1:length(cer.per),.combine='rbind') %:% foreach(j=1:length(per.len),.combine='c') %dopar% (bpcalc(cp=cer.per[[i]], pl=per.len[[j]], a=alpha[i], b=beta[i]))
  
  pop1.df <- pop2.df <- d_mat
  pop1.df <- pop2.df <- data.frame(AggSite = rownames(d_mat), Pop1 = pop)
  colnames(pop2.df) <- c("AggSite","Pop2")
  m2[lower.tri(m2)] <- NA
  
  
  d_mat_df <- as.data.frame(d_mat)
  s_df <- d_mat_df <- d_mat_df %>% rownames_to_column(var = "Site") %>% 
    tidyr::pivot_longer(!Site,names_to="Site2",values_to ="dist_12") %>% 
    mutate(r_ij = dist_12) %>%
    left_join(pop1.df, by = c("Site" = "AggSite")) %>% 
    left_join(pop2.df, by = c("Site2" = "AggSite"))# %>% 
  s_df %>% group_by(Site) %>%  mutate(test = sum(filter(which(dist_12 < dist_12[row_number()]))),
    s = sum(Pop2[dist_12 > 0]) - sum(Pop2[dist_12 > row_number(dist_12)]))
    
    #Pop2[r_ij < dist_12[Site2]]))
    #row_number(dist_12)
    #df1 %>% [.$dist_12 < r_ij]
    mutate(Matt = rowSums(select(cur_data(), 
                                 where(~ is.numeric(.) &&
                                         (.[Name == 'Matt'] > 25| .[Name == 'Matt'] < -25) ))))  
    
  #upper tri = i_orig
  #lower tri = j_dest
  
  # 
  
  out <- x %>% rownames_to_column(var = "Site1") %>% 
    tidyr::pivot_longer(!Site1,names_to="Site2",values_to ="CDist") %>% 
    filter(CDist <= r) %>% group_by(Site1) %>% 
    left_join(select(df, AggSite, !!sym(e), !!sym(E)), by = c("Site2" = "AggSite")) %>% 
    summarize(LQ = (sum(!!sym(e), na.rm=T)/sum(!!sym(E), na.rm=T)/Denom)) %>% 
    mutate(LQ=ifelse(is.nan(LQ),NA,LQ)) %>% 
    arrange(match(Site1, df$AggSite)) %>% pull(LQ)
  
  for (i in seq_len(length(pop))) {
    # start loop on rows
    for (j in seq_len(length(pop))) {
      # start loop on columns
      if (i == j) next() # skip diagonal of matrix
      
      
      m <- pop[i] # set population value for site i
      n <- pop[j] # set population value for site j
      # find radius as distance between sites i and j
      r_ij <- d_mat[i, j]
      # find all sites within the distance from i to j centered on i
      sel_circle <- which(d_mat[i, ] <= r_ij)
      
      # remove the site i and j from list
      sel_circle <- sel_circle[-which(sel_circle %in% c(i, j))]
      
      s <- sum(pop[sel_circle]) # sum population within radius
      # calculate T_i and output to matrix
      temp <-Ti[i] * ((m * n) / ((m + s) * (m + n + s)))
      
      if (is.na(temp)) temp <- 0
      out[i, j] <- temp
    }
  }
  return(out)
}








radiation <- function(pop, 
                      d_mat, 
                      prob = c("original", "extended"),
                      alpha = NULL,
                      scale = c("invariant", "variant", "none"),
                      commuters = c("pop", "input"),
                      Ti = NULL){
  
  ## create square matrix with rows and columns for every site
  out <- matrix(0, length(pop), length(pop))
    
  for (i in seq_len(length(pop))) {
    # start loop on rows
    for (j in seq_len(length(pop))) {
      # start loop on columns
      if (i == j) next() # skip diagonal of matrix
        
      
      m <- pop[i] # set population value for site i
      n <- pop[j] # set population value for site j
      # find radius as distance between sites i and j
      r_ij <- d_mat[i, j]
      # find all sites within the distance from i to j centered on i
      sel_circle <- which(d_mat[i, ] <= r_ij)
        
      # remove the site i and j from list
      sel_circle <- sel_circle[-which(sel_circle %in% c(i, j))]
        
      s <- sum(pop[sel_circle]) # sum population within radius
      # calculate T_i and output to matrix
      temp <-Ti[i] * ((m * n) / ((m + s) * (m + n + s)))
        
      if (is.na(temp)) temp <- 0
      out[i, j] <- temp
    }
  }
  return(out)
}