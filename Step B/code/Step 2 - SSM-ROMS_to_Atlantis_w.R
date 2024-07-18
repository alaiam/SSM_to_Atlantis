#################################
# Step 2: From SSM regular grid NetCDF outputs to ww .csv file with the good format to before processing to go in atlantis
#################################
# w sign = positive downwards (downwelling), negative upwards (upwelling)

start_time <- Sys.time()

# Packages and settings
library(raster)#
library(ncdf4)
library(tabularaster)
library(rbgm)
library(sf) 
library(tidync)
library(tidyverse)
library(reshape2)



options(dplyr.summarise.inform=FALSE)

#### Var to change to choose the file to read and to write #
year = 2011
velma = T
path        <- "/home/atlantis/SSM_to_Atlantis/step_B"
dir_code    <- "/code"

###########################################################################
# Path and names definition


setwd(path)
output_path <- paste0("/outputWW_",year,"/")
if (!file.exists(output_path)){dir.create(output_path)}

if (velma = T){
  filename <- paste0("regular_grid_variables_",year,"_uvw_velma.nc")
}else{
  filename <- paste0("regular_grid_variables_",year,"_W.nc")  # careful --> this might change in the future
}



###########################################################################
# Read data ROMS data
roms <- tidync(paste(path,"data/",filename, sep = ""))
box_composition <- read.csv(paste0("/box_composition_ww.csv"))


# get list of ROMS variables
roms_vars <- tidync::hyper_grids(roms) %>% # all available grids in the ROMS ncdf
  pluck("grid") %>% # for each grid, pull out all the variables associated with that grid and make a reference table
  purrr::map_df(function(x){
    roms %>% tidync::activate(x) %>% tidync::hyper_vars() %>% 
      dplyr::mutate(grd=x)
  })



box_composition <- read.csv("code/box_composition_ww.csv")
layer_max <- read.csv("code/layer_max_ww.csv")
names(layer_max)[1] <- "Polygon #"
layer_depth <-box_composition %>%
  select("atlantis_layer", "dz.x")%>% 
  mutate(Layer = atlantis_layer-1, Layer_dest = atlantis_layer-1)%>%
  distinct()%>%
  select("dz.x", Layer, Layer_dest)
 
############################################################################################
############################################################################################
############################################################################################
step_file <- c(seq(0,730-30,30),730) #Days to divide the total files
# step_file = seq(0,60,30)

ww_dim <- roms_vars %>% dplyr::filter(name==c("ww")) %>% pluck('grd')

cores=detectCores()
cl <- cores -1 #not to overload your computer
cl <- 2 #not to overload your computer
registerDoParallel(cl)


foreach(month = 2:length(step_file)) %dopar%{
# for (month in 2:length(step_file)){



variable_before_Atlantis <- roms %>%
  tidync::activate(ww_dim) %>%
  tidync::hyper_tibble() %>%
  dplyr::select(ww, longitude, latitude, sigma_layer, time)%>%
  dplyr::rename(
    ww=ww, 
    longitude = longitude,  
    latitude = latitude,  
    roms_layer = sigma_layer, time = time)



variable_before_Atlantis<- variable_before_Atlantis %>% filter(time<=step_file[month]&time>step_file[month-1])



merge_test <- merge(box_composition, variable_before_Atlantis, by = c("latitude", "longitude", "roms_layer"))
rm(variable_before_Atlantis)

###################################################################
time = sort(unique(merge_test$time))          # To adapt after time intergration
box = 89
layer = 6

atlantis_input_ww <- array(rep(NA,box*(layer+1)*length(time)), dim = c((layer+1),box,length(time)))

for (i in 0:(box-1)){
  
  for (t in 1:length(time)){
    all.layers_ww = rep(NA,7)   # define an empty vector to receive the values of the 6 layers for ww
    # Calculate the layer
    for (j in 1:layer){ 
      subset <-merge_test %>%
        filter(.bx0 == i, atlantis_layer == j, time == time[t])
      
      
          if (dim(subset)[1] == 0){
            all.layers_ww[j] = NA
          }else{
            all.layers_ww[j] <- mean(subset$ww, na.rm = T)
          }
    }

    atlantis_input_ww[,i+1,t] <- all.layers_ww
      
    }
}
atlantis_input_ww <- atlantis_input_ww*12*60*60



table_flux_ww <- melt(atlantis_input_ww, varnames = c("Layer", "Polygon #", "time"), value.name = "water.exchange") %>%
  mutate(`Polygon #` = `Polygon #` - 1,
         `adjacent box` = `Polygon #`,
         Layer = Layer - 1,
         Layer_dest = Layer + 1, 
         time = min(merge_test$time) - 1 + time)

table_flux_ww <- table_flux_ww %>% left_join(layer_max, by = "Polygon #") 

table_flux_ww<- table_flux_ww %>%
  filter(Layer_dest <= layer_max) %>%
  select(Layer, `Polygon #` ,time, water.exchange, `adjacent box`, Layer_dest)



################## Hyperdiffusion correction ##############
################## /thickness layer (dz) #######################
# Without these line, the correction is /area of the polygon 

# atlantis_bgm <- read_bgm(paste(path,"PugetSound_89b_070116.bgm", sep = ""))
# ssm <- tidync(paste(path,"ssm_00117_copy_test.nc", sep = ""))
# hyperdiffusion_correction_area <- atlantis_sf %>% ungroup %>% select(area,box_id )
# hyperdiffusion_correction_area <- as.data.frame(hyperdiffusion_correction_area)[,c(1,2)]%>%
#   dplyr::rename(area= area,
#     "Polygon #" = box_id)
# 
# table_flux_ww2 <- table_flux_ww %>% left_join(hyperdiffusion_correction_area, by = "Polygon #") %>%
#   left_join(layer_depth[-3], by = "Layer") %>%
#   left_join(layer_depth[-2], by = "Layer_dest")
# table_flux_ww2 <- table_flux_ww2 %>%
#   mutate(dz = ifelse(water.exchange>0, dz.x.y,dz.x.x))
# 
# table_flux_ww_hyperdiff_corr <- table_flux_ww2 %>%
#   mutate(water.exchange = water.exchange*area/dz)%>%
#   select(Layer, `Polygon #` ,time, water.exchange, `adjacent box`, Layer_dest)

  
output_filename <- paste0("flux_ww_", step_file[month],".csv")
csv_filename <- paste0(output_path, output_filename)
write.csv(table_flux_ww, csv_filename, row.names = F)
}

