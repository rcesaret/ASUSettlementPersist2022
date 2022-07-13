library(ggplot2)
library(ggridges)
library(tidyverse)

ModellingPeriod_AxisBreaks <- readRDS(file=paste0(wd$data_p,"ModellingPeriod_AxisBreaks.RData"))
ModellingPeriod_AxisLabels <- readRDS(file=paste0(wd$data_p,"ModellingPeriod_AxisLabels.RData"))
ModellingPeriods_df <- readRDS(file=paste0(wd$data_p,"ModellingPeriods_df.RData"))

ordering <- c("AggSite","AggID","Period","PeriodType","PeriodNum""SubOccSeqLoc",
  "ComponentNum","Area_ha","Population.s2","Log_Population.s2","UrbanPop.s2",
  "rPert","Urb_rPert","Abs_rPert","LogAbs_rPert","PopBwCont","PopFwCont",
  "Found","FoundInit","Abandon","Persist","DewarType",
  "OccuTime","OccuInertia","UrbOccuTime","UrbOccuInertia",
  "density","connectiv","cntrlz_btw","cntrlz_eig","cntrlz_clo","cntrlz_avg",
  "Catch_Popdens","SetHierLevel") 




Agg2b = left_join(Agg2,MP2[,c(1,3,5:8)],by="Period")

Q <- Agg2b %>% 
  mutate(
    Grow = ifelse(r12_Pert > 0, 1, 0),
    Decline = ifelse(r12_Pert < 0, 1, 0)) %>% 
  group_by(PeriodNum, Period, Begin.y,End,Mid,Length) %>% 
  summarize(
    Sites = n(),
    Pop = sum(Population.s2),
    LogPop = log(sum(Population.s2)),
    UrbPop = sum(UrbanPop.s2),
    LogUrbPop = log(sum(UrbanPop.s2)),
    UrbRatio = UrbPop/Pop,
    MeanPopdens = mean(PopDens.s2),
    GiniPop = ineq(Population.s2, type = "Gini", na.rm = TRUE),
    AtkinsonPop = ineq(Population.s2, type = "Atkinson", na.rm = TRUE),
    TheilPop = ineq(Population.s2, type = "Theil", na.rm = TRUE),
    LAsymPop = Lasym(Population.s2),
    MeanPop = Pop/Sites,
    MeanPert = mean(r12_Pert, na.rm=T),
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
    Decline.PopMean = Decline.Pop/Decline.n) %>% ungroup()

#Urban growth rate
#Pop growth rate
Q$TotPop_r12_Pert <- NA
Q$UrbPop_r12_Pert <- NA
for (i in 1:(nrow(Q)-1)) {
  Q$TotPop_r12_Pert[i] <- ((Q$Pop[i+1]/Q$Pop[i])^(1/(abs((Q$Mid[i]) - (Q$Mid[i+1])))))-1
}
for (i in 1:(nrow(Q)-1)) {
  Q$UrbPop_r12_Pert[i] <- ((Q$UrbPop[i+1]/Q$UrbPop[i])^(1/(abs((Q$Mid[i]) - (Q$Mid[i+1])))))-1
}

QL = Q %>% pivot_longer(
  cols = Sites:UrbPop_r12_Pert,
  names_to = "Metric",
  values_to = "Value",
  values_drop_na = F)

QL2a <- QL %>% filter(Metric == "Pop"|Metric =="UrbPop") %>% mutate(Line = "Pop &\nUrban Pop", Value=Value*3)

QL2b <- QL %>% filter(Metric == "LogPop"|Metric =="LogUrbPop") %>% mutate(Line = "Log Pop &\nUrban Pop", Value = ifelse(!is.finite(Value),0,Value*5000))

QL2c <- QL %>% filter(Metric == "GiniPop") %>% mutate(Line = "Gini Pop", Value=Value*200000)

QL2d <- QL %>% filter(Metric == "Persist.PctPop") %>% mutate(Line = "Persisting &\nFounded\n%Pop", Value=Value*200000)

QL2e <- QL %>% filter(Metric == "Found.PctPop") %>% mutate(Line = "Persisting &\nFounded\n%Pop", Value=Value*200000)

QL2f <- QL %>% filter(Metric == "Abandon.PctPop") %>% mutate(Line = "Abandoned\n%Pop", Value=Value*200000)

QL2g <- QL %>% filter(Metric == "UrbRatio") %>% mutate(Line = "Urbanization\nRatio", Value=Value*250000)

QL2 <- rbind(QL2a,  QL2c, QL2d, QL2e, QL2f, QL2g)#QL2b,

QL2 <- QL2 %>%
  mutate(Line = fct_relevel(Line, levels = "Gini Pop", "Urbanization\nRatio", "Abandoned\n%Pop", "Persisting &\nFounded\n%Pop", "Pop &\nUrban Pop"))



RegionalVars = ggplot(QL2, aes(Mid,Value,  fill=Metric)) + 
  geom_ridgeline(aes(y=Line, height=Value), color="black",alpha=0.5, scale=.000004, min_height=-Inf) +
  geom_vline(data=MP, aes(xintercept = years), linetype = "dashed", alpha=0.6)+
  labs(title="Time Series of Regional-Scale Variables", subtitle = "Southern Basin of Mexico", x ="Years BC/AD", y = "Metric")+
  scale_y_discrete(expand = c(0.01, 0)) +
  scale_x_continuous(limits= c(-1400,1520),breaks = MyBreaks,labels = MyLabs) +
  theme_bw()+
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, color="black", size=9), 
        axis.text.y = element_text(color="black", face="bold"), 
        axis.title.x = element_text(color="black", size=14, face="bold"),
        axis.title.y = element_text(color="black", size=14, face="bold"),
        plot.title = element_text(hjust = 0.5, face="bold", size=16), 
        plot.subtitle = element_text(hjust = 0.5, face="bold", size=14))+
  theme(plot.background = element_rect(fill = "white", colour = NA))

ggsave("RegionalVars.png", plot = RegionalVars, 
       device = "png", path = wd$figs, scale = 1, 
       width = 7, height = 5.5,   units = "in",  dpi = 1500)

rm(QL2a, QL2b, QL2c, QL2d, QL2e, QL2f, QL2g,RegionalVars, QL, Q, Agg2b, QL2, Agg21)

