# This script reproduces the figures and analysis from "Repeated Fires Reduce Plant Diversity in Wyoming Big Sagebrush Ecosystems (1984-2014)" by Adam L. Mahood and Jennifer K. Balch. Figures are saved to the "figures" directory, and data from which everything is derived is in the "Data" directory. 
# Author: Adam Mahood
# Last update: October 2018


# setup ------------------------------------------------------------------------
libs <- c("ggplot2", "raster",  "rgdal", "plyr", "nlme", "broom", "rgeos", "maptools",
          "doBy", "tidyr", "ggpubr", "agricolae",  "lme4", "remef", "vegan", 
          "knitr", "cowplot", "sf", "dplyr", "MuMIn", "r2glmm")
if (!"devtools" %in% rownames(installed.packages())) install.packages("devtools")
if (!"remef" %in% rownames(installed.packages())) devtools::install_github("hohenstein/remef")

iini <-function(x){
  #stands for install if not installed
  if (!x %in% rownames(installed.packages())) install.packages(x)
}
lapply(libs, iini)
lapply(libs, library, character.only = TRUE, verbose = FALSE)

set.seed(43534)

# data prep for map figure------------------------------------------------------
plots = readOGR("Data/Plots/", "FF_plots", verbose = FALSE)
fire_frequency = raster("Data/FF_Sev_2-4_8414.tif")
roads = readOGR("Data/roads/", "tl_2016_32_prisecroads", verbose = FALSE)


# reprojecting
plots_nad83 = sp::spTransform(plots, crs(fire_frequency))
roads_nad83 = sp::spTransform(roads, crs(fire_frequency))

plots_nad83@data$Block = paste(substr(plots_nad83@data$Name, 2,3))
plots_nad83@data$Block = as.factor(ifelse(plots_nad83@data$Block == "03", "05", plots_nad83@data$Block)) 
# block 3 and 5 are merged

# clipping FF to plots

# First expand the exent of the points so that we can clip the raster to 
# an extent that doesn't cut out the point symbols on the edges
extent_fig = extent(plots_nad83) 
extent_fig@xmin = extent_fig@xmin - 5000
extent_fig@xmax = extent_fig@xmax + 5000
extent_fig@ymin = extent_fig@ymin - 5000
extent_fig@ymax = extent_fig@ymax + 5000

fire_frequency_cl = raster::crop(fire_frequency, extent_fig)
roads_cl = raster::crop(roads_nad83, extent_fig)

ff.p = rasterToPoints(fire_frequency_cl)
ff.df =data.frame(ff.p)
colnames(ff.df) = c("easting", "northing", "fire_frequency")
ff.df$fire_frequency = as.factor(ff.df$fire_frequency)
my_colors = c("white", "grey80", "grey60", "grey40", "grey20", "black")

roads_cl@data$id = rownames(roads_cl@data)
roads.points = tidy(roads_cl, region="id")
roads.df = join(roads.points, roads_cl@data, by="id")

plotpoints = as.data.frame(plots_nad83@coords)
plotpoints$`Plot Locations` = "Plot Locations"

### for the box around the main map #

extent_df = as.data.frame(matrix(nrow = 5, ncol =2), col.names = c(x,y))
extent_df$x = c(extent_fig@xmin,extent_fig@xmin,extent_fig@xmax,extent_fig@xmax,extent_fig@xmin)
extent_df$y = c(extent_fig@ymax,extent_fig@ymin,extent_fig@ymin,extent_fig@ymax,extent_fig@ymax)
extent_df = extent_df[,3:4]

# for the zoom-in #

zoom_fig = extent(plots_nad83[plots_nad83$Block == "06",]) 
zoom_fig@xmin = zoom_fig@xmin - 300
zoom_fig@xmax = zoom_fig@xmax + 300
zoom_fig@ymin = zoom_fig@ymin - 300
zoom_fig@ymax = zoom_fig@ymax + 300

ffzoom_cl = raster::crop(fire_frequency, zoom_fig)
#rzoom_cl = raster::crop(roads_nad83, zoom_fig)

ffz.p = rasterToPoints(ffzoom_cl)
ffz.df = data.frame(ffz.p)
colnames(ffz.df) = c("easting", "northing", "fire_frequency")
ffz.df$fire_frequency = as.factor(ffz.df$fire_frequency)
my_colors = c("white", "grey80", "grey60", "grey40", "grey20", "black")

# rzoom_cl@data$id = rownames(rzoom_cl@data)
# rz.points = tidy(rzoom_cl, region="id")
# rz.df = join(rz.points, roads_cl@data, by="id")

pzpoints = as.data.frame(plots_nad83[plots_nad83$Block == "06",]@coords)
pzpoints$`Plot Locations` = "Plot Locations"

zoom_df = as.data.frame(matrix(nrow = 5, ncol =2), col.names = c(x,y))
zoom_df$x = c(zoom_fig@xmin,zoom_fig@xmin,zoom_fig@xmax,zoom_fig@xmax,zoom_fig@xmin)
zoom_df$y = c(zoom_fig@ymax,zoom_fig@ymin,zoom_fig@ymin,zoom_fig@ymax,zoom_fig@ymax)
zoom_df = zoom_df[,3:4]


# for the extent indicator

zoom1_df = list()
for(i in levels(plots_nad83$Block)){
  zoom1_fig = extent(plots_nad83[plots_nad83$Block == i,]) 
  zoom1_fig@xmin = zoom1_fig@xmin - 3000
  zoom1_fig@xmax = zoom1_fig@xmax + 3000
  zoom1_fig@ymin = zoom1_fig@ymin - 3000
  zoom1_fig@ymax = zoom1_fig@ymax + 3000
  
  zoom1_df[[i]] = as.data.frame(matrix(nrow = 5, ncol =2), col.names = c(x,y))
  zoom1_df[[i]]$x = c(zoom1_fig@xmin,zoom1_fig@xmin,zoom1_fig@xmax,zoom1_fig@xmax,zoom1_fig@xmin)
  zoom1_df[[i]]$y = c(zoom1_fig@ymax,zoom1_fig@ymin,zoom1_fig@ymin,zoom1_fig@ymax,zoom1_fig@ymax)
  zoom1_df[[i]] = zoom1_df[[i]][,3:4]
}


# locator data prep #


# ecoregion
CBR = readOGR("Data/ecoregion/", layer = "CBR", verbose = FALSE)
CBR@data$id = rownames(CBR@data)
CBR.points = fortify(CBR, region="id")
CBR.df = join(CBR.points, CBR@data, by="id")


#continental US
CUS = readOGR("Data/CUS/", layer = "CUS_5070", verbose = FALSE)
CUS@data$id = rownames(CUS@data)
CUS.points = fortify(CUS, region="id")
CUS.df = join(CUS.points, CUS@data, by="id")


loc_ext = extent(plots_nad83)
loc_ext@xmin = -2030000
loc_ext@xmax = -1300000
loc_ext@ymin = 1710000
loc_ext@ymax = 2320000

lextent_df = as.data.frame(matrix(nrow = 5, ncol =2), col.names = c(x,y))
lextent_df$x = c(loc_ext@xmin,loc_ext@xmin,loc_ext@xmax,loc_ext@xmax,loc_ext@xmin)
lextent_df$y = c(loc_ext@ymax,loc_ext@ymin,loc_ext@ymin,loc_ext@ymax,loc_ext@ymax)
lextent_df = lextent_df[,3:4]



# creating the map figure ------------------------------------------------------

map1 <- ggplot(ff.df, aes(easting, northing)) +
  geom_raster(aes(fill = fire_frequency)) +
  scale_fill_manual(values = my_colors, name = "Fire Frequency") +
  theme_void() +
  geom_path(data = roads.df, aes(long, lat, group=group), color = "grey40") +
  geom_point(data = plotpoints, aes(x = plotpoints[,1], y = plotpoints[,2]), color = "black") +
  labs(color = "Plot Locations") +
  theme(legend.position = "none")  +
  geom_path(data = zoom1_df[[1]], aes(x=x,y=y), color = "grey50", size = 1) +
  geom_path(data = zoom1_df[[2]], aes(x=x,y=y), color = "grey50", size = 1) +
  geom_path(data = zoom1_df[[3]], aes(x=x,y=y), color = "grey50", size = 1) +
  geom_path(data = zoom1_df[[4]], aes(x=x,y=y), color = "grey50", size = 1) +
  geom_path(data = zoom1_df[[5]], aes(x=x,y=y), color = "grey50", size = 1) +
  geom_path(data = zoom1_df[[7]], aes(x=x,y=y), color = "grey50", size = 1) +
  geom_path(data = zoom1_df[[6]], aes(x=x,y=y), color = "black", size = 1) +
  geom_path(data = extent_df, aes(x=x,y=y)) 
# +
# geom_text(x = -1795942, y = 2200081, label = "Winnemucca")
# 

loc <- ggplot(CBR.df) +
  geom_path(data = CBR.df, aes(long, lat, group=group), color = "grey40") +
  geom_path(data = CUS.df, aes(long, lat, group=group), color = "grey60") +
  geom_path(data = extent_df, aes(x=x,y=y)) +
  coord_cartesian(xlim = c(-2030000, -1300000), ylim = c(1710000,2320000)) +
  theme(panel.border = element_rect()) +
  #geom_path(data = lextent_df, aes(x=x, y=y)) +
  theme_void()

zoom <- ggplot(ffz.df, aes(easting, northing)) +
  geom_raster(aes(fill = fire_frequency)) +
  scale_fill_manual(values = my_colors, name = "Fire Frequency") +
  theme_void() +
  #geom_path(data = rz.df, aes(long, lat, group=group), color = "grey40") +
  geom_point(data = pzpoints, aes(x = pzpoints[,1], y = pzpoints[,2]), color = "black") +
  labs(color = "Plot Locations") +
  theme(legend.position = "none") +
  geom_path(data = zoom_df, aes(x=x,y=y), color = "black", size = 1)# +
# ggtitle("Block 6")

p_map<-ggdraw() +
  draw_plot(map1, 0,0,0.7,1) +
  draw_plot(loc, 0.7, 0.5, 0.3, 0.4) +
  draw_plot(zoom, 0.7, 0.0375, 0.28, 0.4) +
  draw_plot_label(label = c("A", "B", "C"), x = c(0.615,0.92,0.92), y = c(0.95, 0.95, 0.42)) +
  draw_plot_label(label = c("Winnemucca", "Inset (C)", "Study Area (A)"), x = c(0.125,0.12, 0.73), y = c(0.3,0.57,0.79), size = 8) +
  draw_plot_label(label = "NV", x = 0.8, y = 0.7, size = 10, col = "grey50")

ggsave(p_map, filename="figures/figure_2_map.pdf", 
       limitsize = FALSE,bg="white",
       width = 6, height = 4)

# species accumulation ---------------------------------------------------------
# many thanks to:
# http://kembellab.ca/r-workshop/biodivR/SK_Biodiversity_R.html

sac_x <- read.csv("Data/FF_All_plots.csv") %>%
  dplyr::select(-SS1, -SS2) %>%
  mutate(Plot = substr(Image,1,9)) %>%
  mutate(UNK9 = UNK2 + UNK9) %>% # UNK2 = UNK9!
  mutate(LEPE2 = LEPE2 + LEPE2.1) %>%
  dplyr::select(-UNK2, -Image, -LEPE2.1) %>%
  dplyr::rename(CYMOP2 = CYMO,
                BAPR5 = BAPR,
                ERNA10 = ERNA,
                CRYPT = CRYPTANTHASP,
                ELCA13 = ELCA,
                ACMI2 = ACMI) %>%
  group_by(Plot) %>%
  summarise_all(mean, na.rm=TRUE) %>%
  mutate(Seedling_presence = ifelse(Seedlings > 0 ,1,0))%>%
  dplyr::select(-ssMean, -Seedlings,
                -Seedling_presence, -LITTER, -BARE, -SHADOW,
                -UNCLEAR, -WOOD, -SCAT, -DUNG, -ROCK,-CRUST,
                plot=Plot)

sac_y <-read.csv("Data/FF_plot_level.csv") %>%
  dplyr::rename(plot = X) %>%
  dplyr::select(plot, FF) %>%
  left_join(sac_x)

sact <- list()
for(i in 0:3){
  sac_y %>% filter(FF == i) %>% 
    dplyr::select(-FF,-plot) %>% 
    specaccum(method = "exact") -> zero
  sact[[i+1]] <- tibble(sites = zero$sites, 
                 richness = zero$richness, 
                 sd = zero$sd,
                 fire_frequency = as.factor(i))
}
sact <- do.call("rbind", sact)

pd <- position_dodge(0.1)

p_sac<-ggplot(sact, aes(x=sites, y=richness, colour = fire_frequency)) + 
  geom_errorbar(aes(ymin=richness-sd, ymax=richness+sd), width =0, position = pd) +
  #geom_ribbon(aes(ymin=richness-sd, ymax=richness+sd, fill = fire_frequency),alpha = 0.2) +
  geom_line(aes(color = fire_frequency), position=pd) +
  scale_color_discrete(name="Fire\nFrequency") +
  theme_classic()+
  xlab("Sites") +
  ylab("Richness") +
  theme(legend.position = c(0.02,1),
        legend.justification = c(0,1))

ggsave(p_sac,filename = "figures/figure_6_sac_plot.png",
       limitsize = FALSE,
       width = 6,bg="white",
       height = 5)

# NMDS figure ------------------------------------------------------------------

ENV <- read.csv("Data/FF_All_plots.csv") %>%
  dplyr::select(Image, ssMean, Seedlings, Seedling_presence, LITTER, BARE, ROCK,
                CRUST, SHADOW, UNCLEAR, WOOD, SCAT, DUNG) %>%
  mutate(plot = substr(Image, 1, 9)) %>%
  dplyr::select(-Image) %>%
  group_by(plot) %>%
  summarise_all(mean, na.rm=TRUE) %>%
  mutate(Seedling_presence = as.factor(ifelse(Seedling_presence >0, "present", "absent")))


metadata <- read.csv("Data/FF_plot_level.csv") %>%
  dplyr::rename(plot = X)%>%
  mutate(burned = as.factor(ifelse(FF == 0, "unburned", "burned")),
         soil_carbon_gm2 = Soil_BD * Soil_OM * 100,
         soil_nitrogen_gm2 = Soil_BD * Soil_TN * 100,
         BRTE_carbon_gm2 = BRTE_mass/.9 * BRTE_TC/100,
         BRTE_nitrogen_gm2 = BRTE_mass/.9 * BRTE_TN/100,
         Soil_ag_stab = ENV$ssMean,
         Seedling_presence = as.factor(Seedling_presence),
         exotic_cover = EAF + EAG + EPF + EPG,
         native_cover = NAF + NPF + NPG + shrub_cover_PQ,
         adj_TSF = ifelse(is.na(TSF)==TRUE, 100, TSF),
         block = substr(plot, 1, 3),
         block = as.factor(ifelse(block == "F03", "F05", block))
  ) %>%
  tibble::column_to_rownames("plot") %>%
  dplyr::select(FF,BRTE_TN, BRTE_TC, Soil_BD, Soil_OM, Soil_TN,
                netMineralization, Elevation, folded_aspect, AUM_acre,
                soil_carbon_gm2, soil_nitrogen_gm2, Soil_ag_stab,burned
  )

SPs <- read.csv("Data/FF_All_plots.csv") %>%
  dplyr::select(-SS1, -SS2) %>%
  mutate(Plot = substr(Image,1,9)) %>%
  mutate(UNK9 = UNK2 + UNK9) %>% # UNK2 = UNK9!
  mutate(LEPE2 = LEPE2 + LEPE2.1) %>%
  dplyr::select(-UNK2, -Image, -LEPE2.1) %>%
  dplyr::rename(CYMOP2 = CYMO,
                BAPR5 = BAPR,
                ERNA10 = ERNA,
                CRYPT = CRYPTANTHASP,
                ELCA13 = ELCA,
                ACMI2 = ACMI) %>%
  group_by(Plot) %>%
  summarise_all(mean, na.rm=TRUE) %>%
  mutate(Seedling_presence = ifelse(Seedlings > 0 ,1,0))%>%
  dplyr::select(-ssMean, -Seedlings,-Seedling_presence, -LITTER, -BARE, -SHADOW,
                -UNCLEAR, -WOOD, -SCAT, -DUNG, -ROCK,-CRUST) %>%
  tibble::column_to_rownames("Plot")

# Turn percent cover to relative abundance by dividing each value by sample
# total abundance
comm <- decostand(SPs, method = "total")
comm.k.mds <- metaMDS(comm, distance = "kul", trace = 0)

# ggplot it 

scores <- as.data.frame(scores(comm.k.mds))
scores$site <- rownames(scores)
scores$FF <-as.factor(metadata$FF)
scores$burned <-as.factor(metadata$burned)

#https://stackoverflow.com/questions/14711470/plotting-envfit-vectors-vegan-package-in-ggplot2
ef <- envfit(comm.k.mds, SPs, na.rm = T, strata = metadata$block, permutations = 9999)
sp <-as.data.frame(ef$vectors$arrows*sqrt(ef$vectors$r))
species <- as.data.frame(cbind(sp, ef$vectors$pvals)) %>%
  tibble::rownames_to_column("species") %>%
  rename(p =  "ef$vectors$pvals") %>%
  filter(p < 0.05)

bu_txt <- data.frame(x = c(-1.2,0.7), y = c(0.6,-0.7), text = c("unburned", "burned"))
sp_txt <- data.frame(species = c("BRTE", "POSE", "ARTRW8", "ELEL5", "SIAL2", "ERCI6"),
                           x = c(0.95,    -0.5,    -0.7,    -0.7,      0.0,    0.25),
                           y = c(0.07,     0.8,    -0.8,     0.15,     0.8,   -0.55)
                     )

p_nmds<- ggplot(data=scores) +
  stat_ellipse(aes(x=NMDS1, y=NMDS2, color = burned),type = "t", level = 0.95) +
  geom_point(aes(x=NMDS1, y=NMDS2, shape=FF, size= 4)) +
  scale_shape_manual(values = c(0:3), name = "Fire\nFrequency") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  geom_segment(data=species,aes(x=0,xend=NMDS1,y=0,yend=NMDS2),
             arrow = arrow(length = unit(0.5, "cm")),colour="darkgrey",inherit_aes=FALSE) + 
  geom_text(data=sp_txt,aes(x=x,y=y,label=species),size=4, alpha = 0.7) +
  geom_text(data=bu_txt, aes(x=x,y=y,label=text, color = text), fontface = "bold")+
  coord_fixed() +
  guides(size = F, color = F) +
  theme(legend.justification=c(0,0), legend.position=c(0,0),
        legend.background = element_rect(fill = 'transparent'))


ggsave(p_nmds, filename = "figures/figure_3_nmds.png", limitsize = FALSE, width = 7, height = 6,
       dpi=600, bg="white")

# tukey plots ------------------------------------------------------------------
SPs <- read.csv("Data/FF_All_plots.csv") %>%
  dplyr::select(-SS1, -SS2) %>%
  mutate(Plot = substr(Image,1,9)) %>%
  mutate(UNK9 = UNK2 + UNK9) %>% # UNK2 = UNK9!
  mutate(LEPE2 = LEPE2 + LEPE2.1) %>%
  dplyr::select(-UNK2, -Image, -LEPE2.1) %>%
  dplyr::rename(CYMOP2 = CYMO,
                BAPR5 = BAPR,
                ERNA10 = ERNA,
                CRYPT = CRYPTANTHASP,
                ELCA13 = ELCA,
                ACMI2 = ACMI) %>%
  group_by(Plot) %>%
  summarise_all(mean, na.rm=TRUE) %>%
  mutate(Seedling_presence = ifelse(Seedlings > 0 ,1,0))%>%
  dplyr::select(-ssMean, -Seedlings,-Seedling_presence, -LITTER, -BARE, -SHADOW,
                -UNCLEAR, -WOOD, -SCAT, -DUNG, -ROCK,-CRUST) %>%
  tibble::column_to_rownames("Plot")

lspn <- log(specnumber(SPs)) %>% 
  as.data.frame() %>%
  tibble::rownames_to_column("plot") %>%
  dplyr::rename(lspn = ".")

shannon <- diversity(SPs, index = "shannon") %>%
  as.data.frame() %>%
  tibble::rownames_to_column("plot") %>%
  dplyr::rename(shannon_weaver = ".")

simpson <- diversity(SPs, index = "simpson") %>%
  as.data.frame() %>%
  tibble::rownames_to_column("plot") %>%
  dplyr::rename(simpson = ".")

invsimpson <- diversity(SPs, index = "invsimpson") %>%
  as.data.frame() %>%
  tibble::rownames_to_column("plot") %>%
  dplyr::rename(inverse_simpson = ".")

adiv <- left_join(shannon, simpson)%>%
  left_join(invsimpson)%>%
  left_join(lspn) %>%
  mutate(pielou_evenness = shannon_weaver/lspn) %>%
  dplyr::select(-lspn)

# Turn percent cover to relative abundance by dividing each value by sample
# total abundance


metadata_evn <- read.csv("Data/FF_plot_level.csv") %>%
  dplyr::rename(plot = X)%>%
  mutate(exotic_cover = EAF + EAG + EPF + EPG,
         native_cover = NAF + NPF + NPG + shrub_cover_PQ,
         block = substr(plot, 1, 3),
         block = as.factor(ifelse(block == "F03", "F05", block))
  ) %>%
  left_join(adiv)%>%
  dplyr::select(FF,exotic_cover,native_cover,sum_spp,plot, 
                shannon_weaver,
                #simpson, 
                #inverse_simpson, 
                pielou_evenness)%>%
  mutate(FF = as.factor(FF))


comm <- decostand(SPs, method = "total")
comm.k.mds <- metaMDS(comm, distance = "kul", trace = 0)

bdmod_sim <- betadisper(betadiver(comm, 1), metadata_evn$FF)
s <- data.frame(bdmod_sim$distances, bdmod_sim$group) %>%
  dplyr::rename(bd = bdmod_sim.distances, FF = bdmod_sim.group) %>%
  mutate(plot = rownames(.)) %>%
  left_join(metadata_evn)

md_tuk <- dplyr::select(s,-plot)
vars <- names(dplyr::select(md_tuk,-FF))
tuk_grps <- data.frame("FF" = rep(0:3, length(names(md_tuk))-1), "variable" = NA, grp = NA) %>% as_tibble
counter <- 1
for(i in 1:length(vars)){
  f <- formula(paste(vars[i], "~ FF"))
  mod <- aov(f, data = md_tuk)
  grps <- HSD.test(mod,"FF")$groups %>%
    tibble::rownames_to_column("FF")
  for(j in 0:3){
    thing <- grps %>% filter(FF==j)
    add_row(tuk_grps)
    tuk_grps$FF[counter] <- j
    tuk_grps$variable[counter] <- vars[i]
    tuk_grps$grp[counter] <- as.character(thing$groups)
    counter=counter+1
  }
}

tk <- mutate(tuk_grps, v_ff = paste0(FF,variable)) %>%
  dplyr::select(grp,v_ff)

lut1 <- c("shannon_weaver"="C. Shannon-Weaver Dominance",
          #"inverse_simpson"="A. Inverse Simpson",
          "pielou_evenness"="D. Pielou Evenness" ,
          "sum_spp"="E. Number of Species",
          #"simpson"="A. Alpha Diversity (Simpson)",
          "bd" = "F. Beta Diversity (Whittaker)",
          "native_cover" ="A. Native Cover",
          "exotic_cover"="B. Exotic Cover")

evn <- gather(dplyr::select(s,-plot), key = Variable, value = Value, -FF) %>%
  mutate(FF = as.factor(FF)) %>%
  mutate(v_ff = paste0(FF,Variable)) %>%
  left_join(tk) %>%
  as_tibble() %>%
  dplyr::select(-v_ff,
                "Fire Frequency" = FF,
                "TukeyHSD" = grp) %>%
  mutate(Variable = lut1[Variable]) %>%
  arrange(Variable)

p_tuk <- ggplot(evn, aes(x = `Fire Frequency`, y = Value, fill = TukeyHSD)) +
  geom_boxplot() +
  facet_wrap(~Variable, scales = "free", dir = "h", nrow = 2) +
  theme_bw() +
  theme(panel.background = element_rect(fill='white', colour='black')) + 
  scale_fill_manual(values=c("white","grey"))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(legend.position = "None") +
  theme(panel.border = element_rect(fill = NA, colour = "black"))+
  ylab(NULL)

ggsave(p_tuk, filename = "figures/figure_4_tukey_plots.png", limitsize = FALSE, width = 7, height = 5, bg="white")
# functional group figure ------------------------------------------------------

nbp = read.csv("Data/FF_plot_level.csv") %>%
  rename(plot = X) %>%
  dplyr::select(FF,
                "Exotic Annual Forb" = "EAF",
                "Exotic Annual Grass" ="EAG",
                "Exotic Perennial Forb" = "EPF",
                "Exotic Perennial Grass" = "EPG",
                "Native Annual Forb" = "NAF",
                "Native Perennial Forb" = "NPF",
                "Native Perennial Grass" = "NPG",
                "Native Shrub" = "shrub_cover_PQ"
  ) %>%
  mutate(FF = as.factor(FF)) %>%
  group_by(FF) %>%
  summarise_all(mean)%>%
  gather(key = "Origin and Life Form", value = "Percent Cover", -FF) %>%
  rename("Fire Frequency" = FF)

# colors are from RColorBrewer - brewer.pal(8, "Spectral")
p_olf<- ggplot(data = nbp, aes(x = `Fire Frequency`, y = `Percent Cover`, 
                       fill = `Origin and Life Form`)) +
  geom_bar(stat = 'identity', color = 'black') +
  scale_fill_manual(values=c("#D53E4F", "#F46D43", "#FDAE61", "#FEE08B", 
                             "#E6F598", "#ABDDA4", "#66C2A5", "#3288BD")) + 
  theme_bw() +
  theme(panel.background = element_rect(fill='white', colour='black')) +      
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        legend.position = c(1,1),
        legend.justification = c(1,1),
        legend.background = element_rect(fill = "transparent", color = "black")) 

ggsave(p_olf, filename = "figures/figure_5_origin_lf.png", limitsize = FALSE,
       width = 4, height = 7, dpi = 600, bg="white")



# cheatgrass vs alpha diversity ------------------------------------------------

SPs <- read.csv("Data/FF_All_plots.csv") %>%
  dplyr::select(-SS1, -SS2) %>%
  mutate(Plot = substr(Image,1,9)) %>%
  mutate(UNK9 = UNK2 + UNK9) %>% # UNK2 = UNK9!
  mutate(LEPE2 = LEPE2 + LEPE2.1) %>%
  dplyr::select(-UNK2, -Image, -LEPE2.1) %>%
  dplyr::rename(CYMOP2 = CYMO,
                BAPR5 = BAPR,
                ERNA10 = ERNA,
                CRYPT = CRYPTANTHASP,
                ELCA13 = ELCA,
                ACMI2 = ACMI) %>%
  group_by(Plot) %>%
  summarise_all(mean, na.rm=TRUE) %>%
  mutate(Seedling_presence = ifelse(Seedlings > 0 ,1,0))%>%
  dplyr::select(-ssMean, -Seedlings,-Seedling_presence, -LITTER, -BARE, -SHADOW,
                -UNCLEAR, -WOOD, -SCAT, -DUNG, -ROCK,-CRUST) %>%
  tibble::column_to_rownames("Plot")

n_sp <- read.csv("Data/FF_plot_level.csv") %>%
  dplyr::rename(plot = X)%>%
  mutate(FF = as.factor(FF),
         block = substr(plot, 1, 3),
         block = as.factor(ifelse(block == "F03", "F05", block))) %>%
  dplyr::select(plot, sum_spp,FF, block,Elevation)

lspn <- log(specnumber(SPs)) %>% 
  as.data.frame() %>%
  tibble::rownames_to_column("plot") %>%
  dplyr::rename(lspn = ".")

shannon <- diversity(SPs, index = "shannon") %>%
  as.data.frame() %>%
  tibble::rownames_to_column("plot") %>%
  dplyr::rename(shannon_weaver = ".")

simpson <- diversity(SPs, index = "simpson") %>%
  as.data.frame() %>%
  tibble::rownames_to_column("plot") %>%
  dplyr::rename(simpson = ".")

invsimpson <- diversity(SPs, index = "invsimpson") %>%
  as.data.frame() %>%
  tibble::rownames_to_column("plot") %>%
  dplyr::rename(inverse_simpson = ".")

adiv <- left_join(shannon, simpson)%>%
  left_join(invsimpson)%>%
  left_join(lspn) %>%
  mutate(pielou_evenness = shannon_weaver/lspn) %>%
  dplyr::select(-lspn)

comm <- decostand(SPs, method = "total")

# the name of the variable suggests simpson, but it's using whittaker (2nd arguement
# in the betadiver function is 1 for whittaker, would be 22 for simpson)
bdmod_sim <- betadisper(betadiver(comm, 1),n_sp$FF) 
s <- data.frame(bdmod_sim$distances, bdmod_sim$group) %>%
  dplyr::select(beta_diversity = bdmod_sim.distances) %>%
  mutate(plot = rownames(.)) 


# MJ: added a centered, scaled elevation column "scaled_elevation"
brte <- tibble::rownames_to_column(SPs,"plot")%>%
  dplyr::select(plot, BRTE) %>%
  left_join(adiv) %>%
  left_join(n_sp) %>%
  left_join(s) %>%
  dplyr::select(-simpson,-inverse_simpson) %>%
  mutate(scaled_elevation = c(scale(Elevation)))

# lmms, generating R2

get_mr2 <- function(x){
  return(as.numeric(MuMIn::r.squaredGLMM(x))[1] %>% round(2))
}

get_cr2 <- function(x){
  return(as.numeric(MuMIn::r.squaredGLMM(x))[2] %>% round(2))
}

# MJ: using lme4 instead, which is supported by remef
mods<-list()
mods$sw <- lmer(shannon_weaver ~ BRTE + scaled_elevation + (1 | block), data = brte)
mods$pe <- lmer(pielou_evenness ~ BRTE + scaled_elevation + (1| block), data = brte) # better is linear model wo elevation
mods$ns <- lmer(sum_spp ~ BRTE + scaled_elevation + (1 | block), data = brte)
mods$bd <- lmer(beta_diversity ~ BRTE + scaled_elevation + (1 | block), data = brte) # better is linear model w elevation

for (i in 1:length(mods)){print(car::vif(mods[[i]]))} # multicolinearity not an issue


mr2 <- sapply(mods, get_mr2, USE.NAMES = FALSE)
cr2 <- sapply(mods, get_cr2, USE.NAMES = FALSE)

# MJ: adjust shannon weaver for elevation effect
brte$sw_adj <- remef(mods$sw, fix = "scaled_elevation", keep.intercept = TRUE)
brte$pe_adj <- remef(mods$pe, fix = "scaled_elevation", keep.intercept = TRUE)
brte$ns_adj <- remef(mods$ns, fix = "scaled_elevation", keep.intercept = TRUE)
brte$bd_adj <- remef(mods$bd, fix = "scaled_elevation", keep.intercept = TRUE)

# MJ: get a data frame to make predictions, using the mean elevation value
brte_mean_elev <- brte %>%
  mutate(scaled_elevation = 0, 
         Elevation = mean(brte$Elevation))



# MJ: generate predictions from the fitted model using this new data frame
brte$predicted_sw_adj <- predict(mods$sw, newdata = brte_mean_elev)
brte$predicted_pe_adj <- predict(mods$pe, newdata = brte_mean_elev)
brte$predicted_ns_adj <- predict(mods$ns, newdata = brte_mean_elev)
brte$predicted_bd_adj <- predict(mods$bd, newdata = brte_mean_elev)


lut <- c("sw_adj" = "A. Shannon-Weaver Dominance",
         "pe_adj" = "B. Pielou Evenness",
         "ns_adj" = "C. Number of Species",
         "bd_adj" = "D. Beta Diversity (Whittaker)"
)

brte_g1 <- dplyr::select(brte, -shannon_weaver, -Elevation,-scaled_elevation,-pielou_evenness, -sum_spp, 
                         -beta_diversity, -starts_with("predicted")) %>%
  gather(value=value, key=key, -block, -FF,-plot,-BRTE)

brte_ga <- dplyr::select(brte, plot, starts_with("predicted")) %>%
  gather(value = fit, key=key, -plot) %>%
  mutate(key = substr(key, 11,16)) %>%
  left_join(brte_g1) %>%
  mutate(key = lut[key])


# plot it!

p0 <- ggplot(brte_ga, aes(x=BRTE, y=value, group=block)) +
  theme_bw() +
  geom_line(aes(y=fit, color=block), size=0.8, alpha = 0.75, show.legend=F) +
  #geom_smooth(method = "lm", color = "grey40",se=F) +
  geom_point(aes(shape=block,color=block), size=1, stroke =1.5) +
  facet_wrap(~key, scales = "free", dir = "h", nrow = 2) +
  scale_shape_manual(values = c(0:7), name = "Block") +
  scale_color_brewer(palette = "Dark2",name = "Block") +
  theme(panel.border = element_rect(fill = NA, colour = "black"))+
  theme(panel.background = element_rect(fill='white', colour='black')) + 
  ylab("Elevation-Adjusted Values") +
  xlab("Cheatgrass Cover")

ggsave(p0,filename = "figures/figure_7_brte_div.png",width=6, height = 5, limitsize = FALSE, bg="white")

