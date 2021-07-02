# convert shapefiles to point csv
library(sf)
library(tidyverse)

rm(list = ls())

islands = c("HAW", "KAH", "KAL", "KAU", "LAN", "MAI", "MOL", "NII", "OAH")

for (ils in 1:length(islands)) {
  
  # ils = 1

  bottom_type = st_read(paste0("N:/GIS/Projects/CommonMaps/01_Preprocess/MHI/", islands[ils], "/hardsoft/biogeo/", tolower(islands[ils]), "_hs_biogeo_shp.shp"))
  
  p = bottom_type %>% 
    st_as_sf() %>% 
    # filter(!HardSoft %in% c("Unknown")) %>%
    filter(!HardSoft %in% c("Unknown", "Land", "Other")) %>%
    ggplot() + 
    geom_sf(aes(group = HardSoft, fill = HardSoft), color = "NA", show.legend = T) + 
    ggdark::dark_theme_minimal() + 
    # scale_fill_viridis_d() + 
    ggtitle(paste0(islands[ils], "_Bottom_Substrate"))
  
  pdf(paste0("outputs/", islands[ils], "_Bottom_Substrate.pdf"), height = 5, width = 6)
  print(p)
  dev.off()
  
  bottom_type = st_cast(bottom_type,"POINT")
  
  bottom_type <- bottom_type %>%
    mutate(lon = unlist(map(bottom_type$geometry,1)),
           lat = unlist(map(bottom_type$geometry,2))) %>% 
    as.data.frame()
  
  # st_write(bottom_type, paste0("data/", ils_lower, "_hs_biogeo_shp.csv"), layer_options = "GEOMETRY=AS_XY")
  # bottom_type = readr::read_csv(paste0("data/", ils_lower, "_hs_biogeo_shp.csv"))
  
  save(bottom_type, file = paste0("data/", tolower(islands[ils]), "_hs_biogeo_shp.RData"))
  
  print(paste0(islands[ils], "...done..."))
  
}


