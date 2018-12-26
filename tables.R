libs <- c("dplyr", "vegan", "tibble")
lapply(libs, library, character.only = TRUE, verbose = FALSE)

climate_vars <- read.csv("Data/climate_vars.csv") %>%
  dplyr::select(-X)

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
         AUM_ha = AUM_acre/0.404,
         Seedling_presence = as.factor(Seedling_presence),
         exotic_cover = EAF + EAG + EPF + EPG,
         native_cover = NAF + NPF + NPG + shrub_cover_PQ,
         adj_TSF = ifelse(is.na(TSF)==TRUE, 100, TSF),
         block = substr(plot, 1, 3),
         block = as.factor(ifelse(block == "F03", "F05", block))) %>%
  dplyr::select(
    -mean_interval,
    -Seedling_presence,
    -burned,
    -AUM_acre,
    -allotment,
    -nearest_sage.m.,
    -nearest_shrub.m.,
    -soil_carbon_gm2,
    -soil_nitrogen_gm2,
    -BRTE_carbon_gm2,
    -BRTE_nitrogen_gm2,
    -Soil_ag_stab,
    -exotic_cover,
    -native_cover,
    -adj_TSF) %>%
  filter(FF > 0) %>%
  left_join(climate_vars)


SPs_b <- read.csv("Data/FF_All_plots.csv") %>%
  dplyr::select(-SS1, -SS2) %>%
  mutate(Plot = substr(Image,1,9)) %>%
  mutate(UNK9 = UNK2 + UNK9) %>% # UNK2 = UNK9!
  mutate(LEPE2 = LEPE2 + LEPE2.1) %>%
  dplyr::select(-UNK2, -Image, -LEPE2.1,-LUP1) %>%
  dplyr::rename(CYMOP2 = CYMO,
                BAPR5 = BAPR,
                ERNA10 = ERNA,
                CRYPT = CRYPTANTHASP,
                ELCA13 = ELCA,
                ACMI2 = ACMI) %>%
  filter(as.numeric(substr(Plot,5,5))>0) %>%
  group_by(Plot) %>%
  summarise_all(mean, na.rm=TRUE) %>%
  mutate(Seedling_presence = ifelse(Seedlings > 0 ,1,0))%>%
  dplyr::select(-ssMean, -Seedlings,-Seedling_presence, -LITTER, -BARE, -SHADOW,
                -UNCLEAR, -WOOD, -SCAT, -DUNG, -ROCK,-CRUST) %>%
  tibble::column_to_rownames("Plot")

SPs <- read.csv("Data/FF_All_plots.csv") %>%
  dplyr::select(-SS1, -SS2) %>%
  mutate(Plot = substr(Image,1,9)) %>%
  mutate(UNK9 = UNK2 + UNK9) %>% # UNK2 = UNK9!
  mutate(LEPE2 = LEPE2 + LEPE2.1) %>%
  dplyr::select(-UNK2, -Image, -LEPE2.1,-LUP1) %>%
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

#MC = multi.collinear(metadata_mc); MC
#MC = multi.collinear(climate_vars[,-13]); MC
# these were the collinear ones
# [1] "soil_carbon_gm2"   "soil_nitrogen_gm2" "BRTE_carbon_gm2"   "BRTE_nitrogen_gm2" "Soil_ag_stab"      "exotic_cover"     
# [7] "native_cover"      "adj_TSF"     

SPs_b <- decostand(SPs_b, "total")
comm.k.dist.burnt <- vegdist(SPs_b, method = "kul")
bd.burnt <- betadiver(SPs_b, 1)
bd_z.burnt <- betadiver(SPs_b, 24)



perm <- how(nperm = 9999)
setBlocks(perm) = with(metadata, block)

admod_bd_w = adonis(bd.burnt ~ 
                    FF +
                    ppt_1pre *
                    tmax_after + 
                    ppt_2pre +
                    TSF +
                    folded_aspect+ 
                    Elevation,
                  data = metadata,
                  permutations = perm, 
                  method = "kul")

write.csv(round(admod_bd_w$aov.tab,3), "writing/beta_diversity_permanova.csv", na="")

admod_bd_z = adonis(bd_z.burnt ~ 
                        FF +
                        ppt_1pre *
                        tmax_after + 
                        ppt_2pre +
                        TSF +
                        folded_aspect+ 
                        Elevation,
                      data = metadata,
                      permutations = perm, 
                      method = "kul")
write.csv(round(admod_bd_z$aov.tab,3), "writing/beta_diversity_z_permanova.csv", na="")


admod_hc = adonis(SPs_b ~
                    TSF * 
                    FF +
                    vpdmax_during  +
                    tmax_during +
                    tmax_pre +
                    AUM_ha, 
                  data = metadata,
                  permutations = perm, 
                  method = "kul")

write.csv(round(admod_hc$aov.tab,3), "writing/community_comp_permanova.csv", na = "")

# matplot(metadata$vpdmax_during, SPs_b, pch=colnames(SPs_b))
# for(i in 1:ncol(SPs_b)){
#   abline(a=admod_hc$coefficients[,i],b=admod_hc$coefficients[1,i], lty=i, col=i)
# }
# 
# matplot(metadata$TSF, SPs_b)
# for(i in 1:ncol(SPs_b)){
#   abline(a=admod_hc$coefficients[2,i],b=admod_hc$coefficients[1,i], lty=i, col=i)
# }
# 
# matplot(metadata$ppt_1pre, bd.burnt)
# for(i in 1:ncol(admod_bd_w$coef.sites)){
#   abline(a=admod_bd_w$coef.sites[3,i],b=admod_hc$coefficients[1,i], lty=i, col=i)
# }
# 
# matplot(metadata$tmax_after, bd.burnt)
# for(i in 1:ncol(admod_bd_w$coef.sites)){
#   abline(a=admod_bd_w$coef.sites[4,i],b=admod_hc$coefficients[1,i], lty=i, col=i)
# }


