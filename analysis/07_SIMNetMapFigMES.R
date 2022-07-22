# MES Revised SIM Network Map Figure

require(ggsn)

pnum <- 9

titl <- bquote(bold('Radiation'~Model^1~'Estimation of Socioeconomic Flows'))

subtitl <- "Classic Period S. Basin of Mexico (c. AD 100-550)"#"Scale Invariant Original Model with Commuters as Population"

footnote <- bquote('1. Original `radiation` spatial interaction model (Simini et al., 2012) calibrated \n     using settlement population and transport network cost-distance\n** Map displaying the natural log of estimated flows greater than ~7.4 persons '(log('T'[ij]) > 2))

scalecoords = c("x" = 496500, "y" = 2108000)

test <- Radiation_EstFlows(df = Poly_List[[pnum]]@data,
                           var = "Population.s2",
                           d_mat = TrsprtNet_CDmatList[[pnum]], 
                           xcoords = Poly_List[[pnum]]@data$East,
                           ycoords = Poly_List[[pnum]]@data$North,
                           crs_coords = 26914,
                           extended = F, #T
                           alpha = 1.2,
                           scale = "invariant", #c("invariant", "variant", "none")
                           commuters = "var", #c("var", "input", "none")
                           Ti_input = NULL,
                           outputs = c("sflines", "sfpts_netstats"))#,

test$sflines$Log_Tij <- log(test$sflines$Tij)

sflines2 <- test$sflines[test$sflines$Log_Tij > 2,]

base <- ggplot() +  geom_stars(data = Hillshade.s)+
              scale_fill_gradientn(colours = c("turquoise3", "black", "gray70"), 
                  values = scales::rescale(c(-9999, -161, 254)), guide="none")+
              coord_sf(datum = sf::st_crs(26914)) +
              xlim(484500,529378)+ylim(2106800,2150500)+
              theme_void()

outmap <- base +
  geom_sf(data = CatchLims.sf, color="black", size=1, alpha = 0) +
  ggspatial::annotation_north_arrow(location = "tr", which_north = "true",
    pad_x = unit(0.7, "cm"), pad_y = unit(1.0, "cm"), height = unit(2.0, "cm"),
    width = unit(2.0, "cm"),
    style = ggspatial::north_arrow_fancy_orienteering(fill = c("black", "white"),
      line_col = "black", text_size = 16, text_face = "bold", text_col = "black"))+
  #north(data = test$sfpts_netstats, symbol = 1, scale = 0.15)+
  scalebar(data = test$sfpts_netstats, dist = 5, dist_unit = "km", transform = F, 
           model = "WGS84", st.color = "black", height = 0.03, st.dist = 0.03, 
           st.size = 5, st.bottom = F, location = "bottomleft", anchor = scalecoords)+
  ggnewscale::new_scale_fill()+
  geom_sf(data = sflines2, mapping=aes(color=Log_Tij, alpha=Log_Tij), size=1.2) +
  scale_color_viridis(option = "turbo", begin = 0.62, na.value = NA, name = expression('log'('T'[ij]))) +
  scale_alpha_continuous(range = c(0.3, 1), guide="none") +
  #geom_sf(data = Catch.sf, color="black", size=0.3, alpha = 0) +
  geom_sf(data=test$sfpts_netstats, mapping=aes(size=Population.s2), shape=21, fill = "darkblue", color="black", stroke = 1, alpha=1)+
  scale_size(breaks = c(50,500,1500,3000), range=c(2,6), name="Population")+
  #guides(fill = guide_legend(order = 2),size = guide_legend(order = 1))+
  labs(title = titl, subtitle = subtitl, caption = footnote)+
  coord_sf(datum = sf::st_crs(26914)) +
  #xlim(484500,529378)+ylim(2106800,2150500)+
  theme_void() + 
  theme(plot.caption = element_text(hjust=0),
        #legend.position="right", 
        legend.position=c(0.165,0.32),
        legend.justification=c(0.5,0.5),  
        legend.margin = margin(1, 4, 1, 4),
        legend.spacing.y = unit(0.1,"line"),legend.spacing.x = unit(0.55,"line"), 
        legend.text=element_text(size=11),legend.title=element_text(size=14),
        legend.key = element_rect(colour = "transparent", fill = "white"), 
        legend.background = element_rect(colour = 'black', fill = 'white', linetype='solid'),
        plot.title = element_text(hjust = 0.5, face="bold", size=14),
        plot.subtitle = element_text(hjust = 0.5, face="bold", size=14),
        plot.background = element_rect(fill = "white", colour = NA)
  )

ggsave("RadSIMNetMapMES.png", plot = outmap, device = "png", path = wd$figs, 
           scale = 1, width = 7, height = 6,   units = "in",  dpi = 1000)

