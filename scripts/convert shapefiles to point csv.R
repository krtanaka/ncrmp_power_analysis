# convert shapefiles to point csv
library(sf)
bottom_type = st_read("data/oah_hs_biogeo/oah_hs_biogeo_shp.shp")

bottom_type %>% 
  st_as_sf() %>% 
  filter(!HardSoft %in% c("Unknown", "Land", "Other")) %>% 
  ggplot() + 
  geom_sf(aes(group = HardSoft, fill = HardSoft), color = "NA", show.legend = T) + 
  ggdark::dark_theme_minimal()

bottom_type_points = st_cast(bottom_type,"POINT")
st_write(bottom_type, "data/oah_hs_biogeo/oah_hs_biogeo_shp.csv", layer_options = "GEOMETRY=AS_XY")
bottom_type = read_csv("data/oah_hs_biogeo/oah_hs_biogeo_shp.csv")
save(bottom_type, file = 'data/oah_hs_biogeo/oah_hs_biogeo_shp.RData')
