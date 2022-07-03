#### RS_Acoef.R
#### Rudolf Cesaretti, 6/10/2022

#### "RS_Acoef" 
#### 
#### 
#### 
#### 

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
##########################  RS_Acoef  #########################
###############################################################

RS_Acoef <- function(z, 
                     ids,
                     return.data = F,
                     plot.results = T,
                     plots = "Both", #RS, Acoef
                     plot_title = "Rank-Size A-Coefficient",
                     xaxis_title = "Log Rank",
                     yaxis_title = "Log Size"){
  
  #z = c(5000, 3000, 2000, 1500, 1000, 800, 500, 500, 500, 200)
  #ids = c(LETTERS[1:length(X)])
  #setup data
  df <- data.frame(ID = ids, Z = z)
  dfc <- df[complete.cases(df), ]
  dfc$LogZ <- log(dfc$Z)
  dfc$Rank <- rank(-dfc$LogZ, ties.method = "random")
  dfc$LogRank <- log(dfc$Rank)
  
  dfc$C <- max(dfc$LogZ)
  fit <- lm(LogZ ~ I(LogRank) - 1, data = dfc, offset = C)
  
  fitsum <- summary(fit)
  
  dfc$LogZ.pred <- fit$fitted.values
  
  #standardize the log-log plot so that log(S1)-log(Sn)=sqrt(2)
  #and log(n)-log(1)=sqrt(2)
  dfc$LogZ.pred.resc <- rescale(dfc$LogZ.pred, to = c(0, sqrt(2)))
  dfc$LogRank.resc <- rescale(dfc$LogRank, to = c(0, sqrt(2)))
  
  #rescale LogZ
  f=approxfun(range(dfc$LogZ.pred), c(0, sqrt(2)),rule=1)
  x=f(dfc$LogZ)
  q=x[-which(is.na(x))]
  q2=dfc$LogZ[-which(is.na(x))]
  x[which(is.na(x))] <- spline(q2,q,xout=dfc$LogZ[which(is.na(x))])$y
  dfc$LogZ.resc <- x
  
  #Compute the area above the Zipfs law diagonal and below the observed
  #rank-size curve (A1), and then the area below the diagonal and above the
  #empirical data (A2)
  un = match(unique(dfc$LogZ.resc), dfc$LogZ.resc)
  vv = 1:length(dfc$LogZ.resc)
  vv = setdiff(vv,un)
  if (length(vv)>0){
    dfc$LogZ.resc[vv] <- seq(-0.0001,length(vv)*-0.0001, by=-0.0001)+dfc$LogZ.resc[vv]
  }
  
  zipf <- approxfun(dfc$LogRank.resc, dfc$LogZ.pred.resc)
  obs <- approxfun(dfc$LogRank.resc, dfc$LogZ.resc)
  
  # defining x range and dividing it to sections (for example n=500)
  i <- seq(0, sqrt(2), length.out=500)
  
  # calculating the distance between the density curves
  h1 <- obs(i)-zipf(i)
  h2 <- zipf(i)-obs(i)
  h1 <- ifelse(h1 < 0, 0, h1)
  h2 <- ifelse(h2 < 0, 0, h2)
  
  #and using the formula for the area of a trapezoid we add up the areas
  A1 <- trapz(i, h1)  # for the regions where zipf>obs
  A2 <- trapz(i, h2)  # for the regions where obs>zipf
  
  #Finally, compute A as the difference A1 - A2
  A <- A1 - A2 
  
  
  m <- c("Slope", "Slope_se", "Slope_t", "Slope_p", "Intercept", "R2", "n", "Ymin", "Ymax", "A1", "A2", "A")
  v <- c(fitsum$coefficients[1], fitsum$coefficients[2], fitsum$coefficients[3], fitsum$coefficients[4], fit$offset[1], 
         fitsum$adj.r.squared, length(dfc$LogZ.resc), range(dfc$LogZ)[1], range(dfc$LogZ)[2], A1, A2, A)
  model.df <- data.frame(Metric = m, Value = v)
  data.df <- dfc
  out_list <- list()
  j=1
  
  if (plot.results == T){
    
    if (plots == "RS" | plots == "Both"){
      
      x <- c(dfc$LogRank,dfc$LogRank,dfc$LogRank)
      y1 <- -1*dfc$LogRank + fit$offset[1]
      y2 <- fitsum$coefficients[1]*dfc$LogRank + fit$offset[1]
      y3 <- dfc$LogZ
      l1 <- rep("Zipf's Law\n(slope = -1)", length(dfc$LogRank))
      lll <- paste0("OLS Estimate\n(slope = ",round(fitsum$coefficients[1],2),")")
      l2 <- rep(lll, length(dfc$LogRank))
      l3 <- rep("Empirical\nRank-Size", length(dfc$LogRank))
      lines <- data.frame(x = x, y = c(y1, y2, y3), lab = c(l1, l2, l3))
      
      RS_plt <- ggplot() + geom_line(data=lines, aes(x,y,color=factor(lab)), size=1.2) + 
        scale_colour_manual(values = c("black", "firebrick1", "dodgerblue"))+
        guides(color = guide_legend(nrow = 3, byrow = TRUE))+
        geom_point(data=dfc, aes(x=LogRank, y=LogZ), shape=19, size=1.5, color="black") + 
        labs(x = xaxis_title, y = yaxis_title) +
        theme_bw()+
        theme(
          axis.text.x = element_text(color="black", size=12, face="bold"), 
          axis.text.y = element_text(color="black", size=12, face="bold"), 
          axis.title.x = element_text(color="black", size=14, face="bold"),
          axis.title.y = element_text(color="black", size=14, face="bold"),
          legend.position="right", legend.title=element_blank(),
          legend.spacing.y = unit(0.5,"cm"),
          legend.key = element_rect(colour = "transparent", fill = "white"),
          legend.text = element_text(color="black", size=10, face="bold"))
      
    }
    
    if (plots == "Acoef" | plots == "Both"){
      
      dat = data.frame(x=c(i,i), y=c(obs(i),zipf(i)), y2=c(pmax(obs(i),zipf(i)),pmin(obs(i),zipf(i))),
                       line = c(rep("obs",length(i)),rep("zipf",length(i))), h=c(h1,h1))
      
      dat2 = data.frame(x=i, obs=obs(i), zipf=zipf(i), h1=h1, h2=h2)
      
      dat2 <- dat2 %>% rowwise() %>% 
        mutate(ymin = min(obs,zipf), ymax= max(obs,zipf),
               group = ifelse(h1 >= h2, "A1", "A2")) %>% ungroup()
      
      dat2$id = ave(dat2$h2,dat2$group, FUN=seq_along)
      dat2$id2 = cumsum(c(1, abs(diff(dat2$id)) > 1))
      
      xs = c(0, sqrt(2))
      beta = c(sqrt(2), -1)
      ys = cbind(1, xs) %*% beta
      
      a1t= paste0("A1=",round(A1,3))
      a2t= paste0("A2=",round(A2,3))
      at = paste0("A = ",round(A,4))
      
      xaxt <- paste0("Rescaled ", xaxis_title, " [0,sqrt(2)]")
      yaxt <- paste0("Rescaled ", yaxis_title, " [0,sqrt(2)]")
      
      Acoef_plt <- ggplot() + 
        geom_ribbon(data = dat2 %>% filter(obs >= zipf),
                    aes(x=x, ymin = zipf, ymax = obs, fill="dodgerblue"), 
                    alpha=0.7) +
        geom_ribbon(data = dat2 %>% filter(obs <= zipf),
                    aes(x=x, ymin = obs, ymax = zipf, fill = "coral1"), 
                    alpha=0.7) +
        scale_fill_identity(name = at,
                            breaks = c("dodgerblue", "coral1"),
                            labels = c(a1t, a2t),
                            guide = "legend")+
        geom_line(data=dfc, aes(x=LogRank.resc,y=LogZ.resc), size=1.2, color="black") + 
        geom_point(data=dfc, aes(x=LogRank.resc,y=LogZ.resc), shape=19, size=1.7, color="black") +
        geom_segment(aes(x = xs[1], xend = xs[2], y = ys[1], yend = ys[2]), size=1, color="black", lty = "twodash")+
        guides(fill = guide_legend(nrow = 2, byrow = TRUE))+
        labs(x = xaxt, y = yaxt) +
        theme_bw()+
        theme(
          axis.text.x = element_text(color="black", size=12, face="bold"), 
          axis.text.y = element_text(color="black", size=12, face="bold"), 
          axis.title.x = element_text(color="black", size=14, face="bold"),
          axis.title.y = element_text(color="black", size=14, face="bold"),
          legend.position="right", legend.spacing.y = unit(0.5,"cm"),
          legend.key = element_rect(colour = "transparent", fill = "white"),
          legend.title = element_text(color="black", size=14, face="bold"),
          legend.text = element_text(color="black", size=12, face="bold"))
      
    }
    if (plots == "RS"){
      out_plot <- RS_plt + labs(title = plot_title, subtitle = "Rank-Size Plot") + 
        theme(plot.title = element_text(hjust = 0.5, face="bold"), 
              plot.subtitle = element_text(hjust = 0.5, face="bold"))
      
      
    }
    if (plots == "Acoef"){
      out_plot <- Acoef_plt + labs(title = plot_title, subtitle = "A-Coefficient Plot") + 
        theme(plot.title = element_text(hjust = 0.5, face="bold"), 
              plot.subtitle = element_text(hjust = 0.5, face="bold"))
      
      
    }
    if (plots == "Both"){
      RS_plt <- RS_plt + labs(subtitle = "Rank-Size Plot") + 
        guides(color = guide_legend(nrow = 1))+
        theme(plot.subtitle = element_text(hjust = 0.5, face="bold", size=14),
              legend.position = "bottom")
      Acoef_plt <- Acoef_plt + labs(subtitle = "A-Coefficient Plot") + 
        guides(fill = guide_legend(nrow = 1))+
        theme(plot.subtitle = element_text(hjust = 0.5, face="bold", size=14),
              legend.position = "bottom")
      plot_row <- plot_grid(RS_plt, Acoef_plt, align = "h", nrow = 1)
      title <- ggdraw() + 
        draw_label(
          plot_title,
          fontface = 'bold',
          x = 0.5,
          hjust = 0.5,
          size=24
        ) +
        theme(
          # add margin on the left of the drawing canvas,
          # so title is aligned with left edge of first plot
          plot.margin = margin(0, 0, 0, 7)
        )
      out_plot <- plot_grid(
        title,
        plot_row,
        ncol = 1,
        # rel_heights values control vertical title margins
        rel_heights = c(0.1, 1)
      )
      
    }
    
    out_list[[j]] <- out_plot
    j=j+1
    
  }
  
  out_list[[j]] <- model.df
  j=j+1
  
  if (return.data == T){
    out_list[[j]] <- data.df
  }
  
  return(out_list)
  
}


