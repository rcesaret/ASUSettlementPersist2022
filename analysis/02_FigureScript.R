#### 02_FigureScript.R
#### Rudolf Cesaretti, 7/5/2022


# Package names
packages <- c("rgdal", "sp", "sf", "GISTools", "lwgeom", "tidyverse", "tidyr", 
              "era", "ggnewscale", "gridExtra", "cowplot", "datplot", "ggridges", 
              "scales", "ggstatsplot", "hrbrthemes", "ineq")

# Install packages not yet installed
installed_packages <- packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
  install.packages(packages[!installed_packages])
}

# load packages
invisible(lapply(packages, library, character.only = TRUE))

rm(packages,installed_packages)


########################### Calc Bayesian Data for comparison #################################


sites <- unique(Agg_NoOvlp$SubOccSeqLoc)

out.list <- list()

for (i in 1:length(sites)){
  
  site.qv <- Agg[which(Agg$SubOccSeqLoc==sites[i]),]
  site.qv_NoOvlp <- Agg_NoOvlp[which(Agg_NoOvlp$SubOccSeqLoc==sites[i]),]
  
  app.out <- spda(site = site.qv$SubOccSeqLoc,  
                  ID = site.qv$AggID, 
                  ph.sites = site.qv$AggSite,
                  ph.periods = site.qv$Period,
                  cer.type = site.qv_NoOvlp$Period,
                  ct = site.qv_NoOvlp$Tot.Assemb,
                  start = site.qv_NoOvlp$Begin,
                  start.era = site.qv_NoOvlp$Begin.era,
                  end = site.qv_NoOvlp$End,
                  end.era = site.qv_NoOvlp$End.era,
                  pop.ests = site.qv$Population, 
                  pc.input = PC,
                  obs.input = site.qv$Tot.Assemb,
                  interval=1,
                  cutoff = 0.05,
                  min.period = 25, 
                  pc.method = "input", 
                  method = "bayesian",
                  alpha = site.qv_NoOvlp$param1,
                  beta = site.qv_NoOvlp$param2)
  out.list[[i]] <- app.out
  
}

out.tbl = bind_rows(out.list)

out.tbl <- out.tbl[order(out.tbl$AggID),]

AggB <- Agg[order(Agg$AggID),]
identical(AggB$AggSite, out.tbl$AggSite)
identical(AggB$AggID, out.tbl$AggID)

names(out.tbl)[names(out.tbl) == 'Interval'] <- 'PeriodInterval'
names(out.tbl)[names(out.tbl) == 'Log_Population'] <- 'Log_Population.s2'
names(out.tbl)[names(out.tbl) == 'Population'] <- 'Population.s2'
names(out.tbl)[names(out.tbl) == 'UrbanPop'] <- 'UrbanPop.s2'
names(out.tbl)[names(out.tbl) == 'Assemb'] <- 'ApportAssemb'

out.tbl <- out.tbl %>% select(-SubOccSeqLoc,-Period, -AggSite)

namez <- c("Population","PopDens","UrbanScale","UrbanPop","RuralPop", "PctUrban","PctRural")

for (i in 1:length(namez)){
  names(AggB)[names(AggB) == namez[i]] <- paste0(namez[i],".s1")
}

#AggB <- AggB %>% select(-OccuIntertia)
#AggB = AggB[,-c(43:56)]
AggB2 <- left_join(AggB,out.tbl,by="AggID")

#Agg2 <- cbind(Agg,out.tbl)

rm(out.list, out.tbl, sites, site.qv_NoOvlp, app.out, site.qv, out.list, sites)

UrbanThresh = 1500

AggB2 <- AggB2 %>% rowwise() %>% mutate(
  RuralPop.s2 = ifelse(Population.s2 < 1500, Population.s2, UrbanThresh)) %>%
  ungroup() %>% mutate(
    PopDens.s2 = Population.s2 / Area_ha,
    UrbanScale.s2 = Population.s2 / UrbanThresh,
    PctUrban.s2 = UrbanPop.s2 / Population.s2,
    PctRural.s2 = RuralPop.s2 / Population.s2)

################################Integrate MeanOccuProb and Bayesian Data  ########################################

AggB2$Method <- "Bayesian"
AggA2 <- Agg2
AggA2$Method <- "MeanOccuProb"
colnames(AggB2) <- colnames(AggA2)
MethodCompare.df = rbind(AggA2,AggB2)


############ OccuProb_Method_PeriodType.png ##################

OccuProb_Method_PeriodType <- grouped_ggbetweenstats(
  data = MethodCompare.df,
  x = Method,
  y = MeanOccuProb,
  grouping.var = PeriodType,
  xlab = "SPDA Method",
  ylab = "Probability (Posterior or MeanOccuProb)",
  annotation.args = list(title = "Occupational Probabilities by SPDA Method and Period Type", subtitle="Overlap/Transition Periods vs Main Ceramic Complex Phases")
)+theme(plot.background = element_rect(fill = "white", colour = NA))

#save figure
ggsave("OccuProb_Method_PeriodType.png", plot = OccuProb_Method_PeriodType, 
       device = "png", path = wd$figs, scale = 1, 
       width = 10.25, height = 6,   units = "in",  dpi = 1800)

rm(OccuProb_Method_PeriodType,MethodCompare.df,AggA2,AggB2)

############## OldNewPopsPlot.png ##################

Agg21 <- Agg2
Agg21$Log_Population.s1 <- log(Agg21$Population.s1)

OldNewPopsPlot <- ggscatterstats(
  data = Agg21, ## dataframe from which variables are taken
  x = Log_Population.s1, ## predictor/independent variable
  y = Log_Population.s2, ## dependent variable
  xlab = "Log Population, Step 2", ## label for the x-axis
  ylab = "Log Population, Step 1", ## label for the y-axis
  xfill = "#CC79A7", ## fill for marginals on the x-axis
  yfill = "#009E73", ## fill for marginals on the y-axis
  title = "Step #1 vs. New Step #2 Population Estimates"
)+theme(plot.background = element_rect(fill = "white", colour = NA))

#save figure
#ggsave("OldNewPopsPlot.png", plot = OldNewPopsPlot, 
#       device = "png", path = wd$figs, scale = 1, 
#       width = 6.5, height = 6,   units = "in",  dpi = 1500)

rm(OldNewPopsPlot)

##################### DurationHists.png #########################

Agg212 <- Agg21 %>% group_by(SubOccSeqLoc) %>% summarize(
  count = n(),
  Length = sum(PeriodLength),
  maxPop = max(Population.s2),
  sumPop = sum(Population.s2))

h1 = ggplot(Agg212, aes(count)) +
  geom_histogram(color = "#000000", fill = "#0099F8", bins = 9) +
  labs(x = "Modelling Periods (n)", y = "Frequency") +
  theme_classic()+
  theme(axis.title.x = element_text(color="black", size=14, face="bold"),
        axis.title.y = element_text(color="black", size=14, face="bold"),
        axis.text.y = element_text(face="bold", color="black", size=12),
        axis.text.x = element_text(face="bold", color="black", size=12))

h2 = ggplot(Agg212, aes(Length)) +
  geom_histogram(color = "#000000", fill = "orangered", bins = 9) +
  labs(x = "Length of Occupation (Years)", y = "Frequency") +
  theme_classic()+
  theme(axis.title.x = element_text(color="black", size=14, face="bold"),
        axis.title.y = element_text(color="black", size=14, face="bold"),
        axis.text.y = element_text(face="bold", color="black", size=12),
        axis.text.x = element_text(face="bold", color="black", size=12))

plot_row <- plot_grid(h1, h2)

title <- ggdraw() + draw_label("SubOccSeqLoc Site Occupational Duration in the SBOM",
    fontface = 'bold',x = 0.5,hjust = 0.5) +
  theme(plot.margin = margin(0, 0, 0, 7))

DurationHists <- plot_grid(title, plot_row,ncol = 1,  rel_heights = c(0.1, 1))+
  theme(plot.background = element_rect(fill = "white", colour = NA))

ggsave("DurationHists.png", plot = DurationHists, 
       device = "png", path = wd$figs, scale = 1, 
       width = 6, height = 3,   units = "in",  dpi = 1500)

rm(DurationHists, title, plot_row, h1, h2, Agg212)

######################### UrbanPopHists1500.png #############################

hpop = Agg21 %>% filter(Population.s2 > 1499) %>% 
  group_by(SubOccSeqLoc) %>% summarize(
    p1000 = n(),
    l1000 = sum(PeriodLength))


h3 = ggplot(hpop, aes(p1000)) +
  geom_histogram(color = "#000000", fill = "aquamarine", bins = 7) +
  labs(x = "Num Periods with Pop > 1500", y = "Frequency") +
  theme_classic()+
  theme(axis.title.x = element_text(color="black", size=12, face="bold"),
        axis.title.y = element_text(color="black", size=14, face="bold"),
        axis.text.y = element_text(face="bold", color="black", size=12),
        axis.text.x = element_text(face="bold", color="black", size=12))

h4 = ggplot(hpop, aes(l1000)) +
  geom_histogram(color = "#000000", fill = "bisque", bins = 11) +
  labs(x = "Length of Periods with Pop > 1500", y = "Frequency") +
  theme_classic()+
  theme(axis.title.x = element_text(color="black", size=12, face="bold"),
        axis.title.y = element_text(color="black", size=14, face="bold"),
        axis.text.y = element_text(face="bold", color="black", size=12),
        axis.text.x = element_text(face="bold", color="black", size=12))

plot_row <- plot_grid(h3, h4)

title <- ggdraw() + draw_label("SubOccSeqLoc Site Urban Pop Duration in the SBOM",
    fontface = 'bold',x = 0.5,hjust = 0.5) +   theme(plot.margin = margin(0, 0, 0, 7))
subtitle <- ggdraw() + draw_label("n = 40 (of 1417)",  fontface = 'bold', x = 0.5,hjust = 0.5,size=12) + 
  theme(plot.margin = margin(0, 0, 0, 7))


UrbanPopHists1500 <- plot_grid(title, subtitle, plot_row,ncol = 1,rel_heights = c(0.08, 0.08, 1))+
  theme(plot.background = element_rect(fill = "white", colour = NA))

ggsave("UrbanPopHists1500.png", plot = UrbanPopHists1500, 
       device = "png", path = wd$figs, scale = 1, 
       width = 6, height = 3,   units = "in",  dpi = 1500)

rm(UrbanPopHists1500,subtitle, title, plot_row, h4, h3, hpop)

########################## PertBoxes.png ##############################

MP2$Width = MP2$Length/2
Agg2b = left_join(Agg2,MP2[,c(3,5:8)],by="Period")

df_mean <- Agg2b %>% filter(PeriodNum<17) %>% 
  group_by(PeriodNum) %>% 
  summarize(average = mean(r12_Pert)) %>%
  ungroup()

labz = c("EF\nto\nEF_MF",       "EF_MF\nto\nMF",    "MF\nto\nMF_LF",       "MF_LF\nto\nLF",    "LF\n to\nLF_TF",       "LF_TF\nto\nTF",    "TF\nto\nTF_CL",       "TF_CL\nto\nCL",    "CL\nto\nCL_ET",      
         "CL_ET\nto\nET",    "ET\nto\nET_LTAzI",       "ET_LTAzI\nto\nLTAzI", "LTAzI\nto\nLTAzI_EA",    "LTAzI_EA\n to\nEA", "EA\nto\nEA_LA",       "EA_LA\nto\nLA")  

xxx = Agg2b %>% filter(PeriodNum<17)

PertBoxes = ggplot() + 
  geom_boxplot(data=xxx, mapping = aes(x=as.factor(PeriodNum), y=r12_Pert))+
  geom_point(data = df_mean, 
             mapping = aes(x = PeriodNum, y = average),
             color="blue") +
  geom_hline(yintercept=0, color="red",linetype = "solid", size=1)+
  geom_line(data = df_mean, 
            mapping = aes(x =PeriodNum, y = average, group=1),color="blue", size=1.2)+
  scale_x_discrete(labels = labz)+
  labs(title="Projected Exponential Population Growth Rate Between Periods",x= "Time Period", y="Popularion Growth Rates (Pe^rt)")+
  theme_ipsum() +
  theme(axis.text.x = element_text(color="black", size=9), 
        axis.text.y = element_text(color="black"), 
        axis.title.x = element_text(color="black", size=16, face="bold"),
        axis.title.y = element_text(color="black", size=12, face="bold"),
        plot.title = element_text(hjust = 0.5, face="bold",size=14))+
  theme(plot.background = element_rect(fill = "white", colour = NA))

ggsave("PertBoxes.png", plot = PertBoxes, 
       device = "png", path = wd$figs, scale = 1, 
       width = 8.5, height = 5,   units = "in",  dpi = 1600)

rm(PertBoxes,labz,df_mean,Agg2b,xxx)

####################### RegionalVars.png ######################


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




