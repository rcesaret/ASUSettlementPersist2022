#### Radiation_EstFlows.R
#### Rudolf Cesaretti, 6/31/2022

#### "Radiation_EstFlows" 
#### 
#### 
#### 
#### https://book.archnetworks.net/spatialinteraction

pak <- c("rgdal", "sp", "sf", "GISTools", "lwgeom", "tidyverse", "tidyr", "data.table", "zoo")
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


###############################################################
##########################  Radiation_EstFlows  #########################
###############################################################






radiation <- function(pop, d_mat) {
  ## create square matrix with rows and columns for every site
  out <-
    matrix(0, length(pop), length(pop))
  for (i in seq_len(length(pop))) {
    # start loop on rows
    for (j in seq_len(length(pop))) {
      # start loop on columns
      if (i == j)
        next()
      # skip diagonal of matrix
      m <- pop[i] # set population value for site i
      n <- pop[j] # set population value for site j
      # find radius as distance between sites i and j
      r_ij <-
        d_mat[i, j]
      # find all sites within the distance from i to j centered on i
      sel_circle <-
        which(d_mat[i, ] <= r_ij)
      # remove the site i and j from list
      sel_circle <-
        sel_circle[-which(sel_circle %in% c(i, j))]
      s <- sum(pop[sel_circle]) # sum population within radius
      # calculate T_i and output to matrix
      temp <-
        pop[i] * ((m * n) / ((m + s) * (m + n + s)))
      if (is.na(temp)) temp <- 0
      out[i, j] <- temp
    }
  }
  return(out)
}