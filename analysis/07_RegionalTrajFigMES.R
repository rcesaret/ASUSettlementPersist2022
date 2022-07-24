library(ggplot2)
library(ggridges)
library(tidyverse)
library(ineq)
source(paste0(wd$funcs,"RS_Acoef.R"))


ModellingPeriod_AxisBreaks <- readRDS(file=paste0(wd$data_p,"ModellingPeriod_AxisBreaks.RData"))
ModellingPeriod_AxisLabels <- readRDS(file=paste0(wd$data_p,"ModellingPeriod_AxisLabels.RData"))
ModellingPeriods_df <- readRDS(file=paste0(wd$data_p,"ModellingPeriods_df.RData"))

want <- c("AggSite","AggID","Period","PeriodType","PeriodNum","SubOccSeqLoc",
  "ComponentNum","Area_ha","Population.s2","Log_Population.s2","UrbanPop.s2",
  "PopDens.s2","r12_Pert","rPert","Urb_rPert","Abs_rPert","LogAbs_rPert","PopBwCont",
  "PopFwCont","AreaBwCont","AreaFwCont","Found","FoundInit","Abandon","Persist",
  
  "OccuTime","OccuInertia","UrbOccuTime","UrbOccuInertia",
  "density","connectiv","cntrlz_btw","cntrlz_eig","cntrlz_clo","cntrlz_avg",
  "Catch_Popdens","SetHierLevel") 

df <- All_AggPoly@data %>% 
        select(!!!syms(want)) %>% 
        left_join(ModellingPeriods_df[,c(1,3,5:7)],by="Period") %>% 
        mutate(Grow = ifelse(rPert > 0, 1, 0),
               Decline = ifelse(rPert < 0, 1, 0)) %>% 
        group_by(PeriodNum,Period,Begin,End,Mid,Length) %>% 
        summarize(Sites = n(),
                  Pop = sum(Population.s2),
                  MaxPop = max(Population.s2),
                  sdPop = sd(Population.s2),
                  MedPop = median(Population.s2),
                  MeanPop = mean(Population.s2),
                  LogPop = log(sum(Population.s2)),
                  UrbPop = sum(UrbanPop.s2),
                  LogUrbPop = log(sum(UrbanPop.s2)),
                  UrbRatio = UrbPop/Pop,
                  MeanPopdens = mean(PopDens.s2),
                  GiniPop = ineq(Population.s2, type = "Gini", na.rm = TRUE),
                  AtkinsonPop = ineq(Population.s2, type = "Atkinson", na.rm = TRUE),
                  TheilPop = ineq(Population.s2, type = "Theil", na.rm = TRUE),
                  LAsymPop = Lasym(Population.s2),
                  MeanPert = mean(rPert, na.rm=T),
                  sdPert = sd(rPert, na.rm=T),
                  MedPert = median(rPert, na.rm=T),
                  MaxPert = max(rPert, na.rm=T),
                  MinPert = min(rPert, na.rm=T),
                  MeanAbsPert = mean(rPert, na.rm=T),
                  sdAbsPert = sd(rPert, na.rm=T),
                  MedAbsPert = median(rPert, na.rm=T),
                  MaxAbsPert = max(rPert, na.rm=T),
                  MinAbsPert = min(rPert, na.rm=T),
                  MeanLogAbsPert = mean(rPert, na.rm=T),
                  sdLogAbsPert = sd(rPert, na.rm=T),
                  MedLogAbsPert = median(rPert, na.rm=T),
                  MaxLogAbsPert = max(rPert, na.rm=T),
                  MinLogAbsPert = min(rPert, na.rm=T),
                  Persist.n = sum(Persist),
                  Found.n = sum(Found),
                  FoundInit.n = sum(FoundInit),
                  Abandon.n = sum(Abandon),
                  Persist.Pct = Persist.n/Sites,
                  Found.Pct = Found.n/Sites,
                  FoundInit.Pct = FoundInit.n/Sites,
                  Abandon.Pct = Abandon.n/Sites,
                  Persist.Pop = sum((Persist*Population.s2)),
                  Found.Pop = sum((Found*Population.s2)),
                  FoundInit.Pop = sum((FoundInit*Population.s2)),
                  Abandon.Pop = sum((Abandon*Population.s2)),
                  Persist.PctPop = Persist.Pop/Pop,
                  Found.PctPop = Found.Pop/Pop,
                  FoundInit.PctPop = FoundInit.Pop/Pop,
                  Abandon.PctPop = Abandon.Pop/Pop,
                  Persist.PopMean = Persist.Pop/Persist.n,
                  Found.PopMean = Found.Pop/Found.n,
                  FoundInit.PopMean = FoundInit.Pop/FoundInit.n,
                  Abandon.PopMean = Abandon.Pop/Abandon.n,
                  Grow.n = sum(Grow),
                  Decline.n = sum(Decline), 
                  Grow.Pct = Grow.n/Sites,
                  Decline.Pct = Decline.n/Sites, 
                  Grow.Pop = sum((Grow*Population.s2)),
                  Decline.Pop = sum((Decline*Population.s2)), 
                  Grow.PctPop = Grow.Pop/Pop,
                  Decline.PctPop = Decline.Pop/Pop,
                  Grow.PopMean = Grow.Pop/Grow.n,
                  Decline.PopMean = Decline.Pop/Decline.n,
                  TotPopFwCont = sum((PopFwCont*Population.s2/Pop), na.rm=T),
                  TotPopBwCont = sum((PopBwCont*Population.s2/Pop), na.rm=T),
                  TotAreaFwCont = sum((AreaFwCont*Population.s2/sum(Area_ha)), na.rm=T),
                  TotAreaBwCont = sum((AreaBwCont*Population.s2/sum(Area_ha)), na.rm=T),
                  AvgOccuTime = mean(OccuTime,na.rm=T),
                  sdOccuTime = sd(OccuTime,na.rm=T),
                  AvgUrbOccuTime = mean(UrbOccuTime,na.rm=T),
                  TotOccuInertia = sum(OccuInertia,na.rm=T),
                  TotUrbOccuInertia = sum(UrbOccuInertia,na.rm=T),
                  TranspNet_Density = mean(density,na.rm=T),
                  TranspNet_connectiv = mean(connectiv,na.rm=T),
                  TranspNet_cntrlz_btw = mean(cntrlz_btw,na.rm=T),
                  TranspNet_cntrlz_eig = mean(cntrlz_eig,na.rm=T),
                  TranspNet_cntrlz_clo = mean(cntrlz_clo,na.rm=T),
                  TranspNet_cntrlz_avg = mean(cntrlz_avg,na.rm=T),
                  AvgCatchPopdens = mean(Catch_Popdens,na.rm=T),
                  GiniSetHierLevel = ineq(SetHierLevel, type = "Gini", na.rm = TRUE),
                  AtkinsonSetHierLevel = ineq(SetHierLevel, type = "Atkinson", na.rm = TRUE),
                  TheilSetHierLevel = ineq(SetHierLevel, type = "Theil", na.rm = TRUE),
                  LAsymSetHierLevel = Lasym(SetHierLevel),
                  MeanSetHierLevel = mean(SetHierLevel, na.rm = TRUE),
                  sdSetHierLevel = sd(SetHierLevel, na.rm = TRUE),
                  nSetHierLevel = max(SetHierLevel, na.rm = TRUE),
                  MeanrPert12 = mean(r12_Pert, na.rm=T),
                  MedrPert12 = median(r12_Pert, na.rm=T)) %>% ungroup()
  
  
nam <- df$Period

RSA_df_list <- list()
RSA_plot_list <- list()

for (i in 1:length(Poly_List)){
  
  tmp <- RS_Acoef(z=Poly_List[[i]]$Population.s2, ids=Poly_List[[i]]$AggSite, plot_title = paste0(i,". ",nam[i]),yaxis_title = "Log Population")
  
  RSA_plot_list[[i]] <- tmp[[1]]
  
  x <- tmp[[2]]
  xx=as.data.frame(t(as.matrix(x$Value)))
  colnames(xx) = x$Metric
  xx$Period <- nam[i]
  xx$PeriodNum <- i
  xx$meanlog <- mean(Poly_List[[i]]$Log_Population.s2, na.rm=T)
  xx$sdlog <- sd(Poly_List[[i]]$Log_Population.s2, na.rm=T)
  RSA_df_list[[i]] <- xx
}

names(RSA_plot_list) <- nam
RSA_df = do.call(rbind, RSA_df_list)

rm(x,xx,tmp)

df$RS_Slope <- as.numeric(RSA_df$Slope)
df$A1 <- as.numeric(RSA_df$A1)
df$A2 <- as.numeric(RSA_df$A2)
df$A <- as.numeric(RSA_df$A)
df$A1A2 <- as.numeric(RSA_df$A1)+as.numeric(RSA_df$A2)
df$RS_R2 <- as.numeric(RSA_df$R2)





#Urban growth rate
#Pop growth rate
df$TotPop_r12_Pert <- NA
df$len <- NA
#Q$UrbPop_r12_Pert <- NA
for (i in 1:(nrow(df)-1)) {
  df$len[i] <- abs((df$Mid[i]) - (df$Mid[i+1]))
  df$TotPop_r12_Pert[i] <- ((df$Pop[i+1]/df$Pop[i])^(1/(abs((df$Mid[i]) - (df$Mid[i+1])))))-1
}
#for (i in 1:(nrow(Q)-1)) {
#  Q$UrbPop_r12_Pert[i] <- ((Q$UrbPop[i+1]/Q$UrbPop[i])^(1/(abs((Q$Mid[i]) - (Q$Mid[i+1])))))-1
#}


df2 = df %>% select(PeriodNum:Mid,Pop,UrbPop,UrbRatio,AtkinsonPop,nSetHierLevel,Grow.PctPop,
                    AvgOccuTime,TotOccuInertia,Persist.PctPop,Found.PctPop) %>%
  pivot_longer(
  cols = Pop:Found.PctPop,
  names_to = "Metric",
  values_to = "Value",
  values_drop_na = F)%>%
  arrange(Metric,PeriodNum)#%>%
  #mutate(Set = Metric) %>%
  #rowwise()%>%
  #mutate(Set = ifelse(Metric == "UrbPop", 7, Set))%>%
  #ungroup()
ModellingPeriods_df2 = ModellingPeriods_df
ModellingPeriods_df2[1,1] <- -1375

labz = data.frame(
  Metric=c("Pop","UrbPop","UrbRatio","AtkinsonPop","nSetHierLevel","Grow.PctPop",
  "AvgOccuTime","TotOccuInertia","Persist.PctPop","Found.PctPop"),
  MyLabels=c("Total Population\n& Urban Population","Total Population\n& Urban Population","Urbanization Ratio\n(Urban Pop / Total Pop)",
             "Settlement Population Inequality\n(Atkinson Index)","Settlement Hierarchy\n(Ordinal Levels)",
             "Growing Settlements\n(Percent of Population)",
             "Average Occupation Length\n(Years Since Foundation)","Integral of Settlement Population\n(Person Years since Founding)",
             "Persisting Settlements\n(Percent of Population)","Newly Founded Settlements\n(Percent of Population)"),
  y=c(0,0,0,0,0,0,0,0,0,0))

df2=df2%>%left_join(labz, by="Metric")%>% 
  mutate(MyLabels=factor(MyLabels,levels=c("Total Population\n& Urban Population","Urbanization Ratio\n(Urban Pop / Total Pop)",
                                           "Settlement Population Inequality\n(Atkinson Index)","Settlement Hierarchy\n(Ordinal Levels)",
                                           "Growing Settlements\n(Percent of Population)",
                                           "Average Occupation Length\n(Years Since Foundation)","Integral of Settlement Population\n(Person Years since Founding)",
                                           "Persisting Settlements\n(Percent of Population)","Newly Founded Settlements\n(Percent of Population)")))

#Theil Generalised Entropy Inequality Index for Settlement Population
ModellingPeriod_AxisLabels2 = ModellingPeriod_AxisLabels
ModellingPeriod_AxisLabels2[c(5,7,9,11,13,15,17)] = ""
#
  ggp=ggplot()+
  geom_rect(data=ModellingPeriods_df2[c(1,3,5,7,9,11,13,15,17),], 
            aes(xmin=Begin, xmax=End, ymin=-Inf, ymax=Inf),alpha=0.2,color="gray80")+
  geom_ridgeline(data=df2[df2$Metric != "UrbPop",], aes(x=Mid, y=y,height=Value,fill=MyLabels, 
                                                        group=MyLabels),color="black",alpha=0.7)+
  #scale_fill_viridis(option = "G", discrete = TRUE, na.value = NA)+
  #geom_point(data=df2[df2$Metric != "UrbPop",], aes(x=Mid, y=Value,group = MyLabels),color="black", size=0.75)+
  guides(fill="none")+
  scale_x_continuous(limits= c(-1400,1520),breaks = ModellingPeriod_AxisBreaks,labels = ModellingPeriod_AxisLabels2)+
  facet_wrap(~MyLabels, ncol=3, scales='free_y')+
  labs(y='',x= "Time (Years BC/AD)", title="Time Series of Selected Regional-Scale Variables",
       subtitle="Southern Basin of Mexico (SBOM), c.1600 BC - AD 1520")+
  theme_bw()+
  theme(strip.text.x = element_text(color="black", size=11, face="bold"),
    legend.position = "none",strip.background = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, color="black", size=9), 
        axis.text.y = element_text(color="black"), 
        axis.title.x = element_text(color="black", size=14, face="bold"),
        axis.title.y = element_text(color="black", size=14, face="bold"),
        plot.subtitle = element_text(hjust = 0.5, face="bold", size=18),
        plot.title = element_text(hjust = 0.5, face="bold", size=18), 
        plot.background = element_rect(fill = "white", colour = NA))

  ggp=ggp+geom_ridgeline(data=df2[df2$Metric == "UrbPop",], aes(x=Mid, y=y,height=Value,group=MyLabels),color="black",fill= "indianred3",alpha=0.5)
  
  ggpmap<-ggsave("RegionalTimeSeries.png", plot = ggp, device = "png", path = wd$figs, scale = 1, width = 9, height = 7.5,   units = "in",  dpi = 1200)


rm(df2, df, RSA_df_list, RSA_plot_list, RSA_df, ggp, labz)

