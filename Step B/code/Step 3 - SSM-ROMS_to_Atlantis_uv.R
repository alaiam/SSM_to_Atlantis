#################################
# Step 3: from SSM regular grid u and v correct format to create Atlantis outputs
# in a daily format
#################################

start_time <- Sys.time()

# Packages and settings
library(raster)
library(tidyverse)
library(sf) 
library(here)
library(tidync)
library(rbgm)
library(tabularaster)
library(ncdf4)
library(doParallel)

select <- dplyr::select
here <- here::here
map <- purrr::map
options(dplyr.summarise.inform=FALSE)

#### Var to change to choose the file to read and to write #
year = 2011
velma = T
path        <- "/home/atlantis/SSM_to_Atlantis/step_B"


###########################################################################
# Path and names definition


setwd(path)
output_path <- paste0("/outputUV_",year,"/")
if (!file.exists(output_path)){dir.create(output_path)}

if (velma = T){
  filename_u <- paste0("regular_grid_variables_",year,"_uvw_velma.nc")
  filename_v <- paste0("regular_grid_variables_",year,"_uvw_velma.nc")
}else{
  filename_u <- paste0("regular_grid_variables_",year,"_uv_velma.nc")
  filename_v <- paste0("regular_grid_variables_",year,"_uv_velma.nc")  # careful --> this might change in the future
}

###########################################################################
# Load grid information

atlantis_bgm <- read_bgm(paste(path,"PugetSound_89b_070116.bgm", sep = ""))
atlantis_sf <- atlantis_bgm %>% box_sf()
###########################################################################


full_face_composition <- read.csv("code/face_composition_uv.csv")
names(full_face_composition)[c(9,10,11)] <- c("Polygon #","Face #", "adjacent box")
uxy2 <- read.csv("code/uxy2_uv_code.csv")

list <- list()
for (i in 1:length(full_face_composition$uvec)){
  chaine <- gsub("^c\\(", "", full_face_composition$uvec[i])
  chaine <- gsub("\\)$", "", chaine)
  chaine <- unlist(strsplit(chaine, ","))
  if(length((grep(":", chaine)==1))==1){
    chaine <- unlist(strsplit(chaine, ":"))
    chaine <- as.numeric(chaine)
    chaine <- chaine[1]:chaine[2]
  }
  chaine <- as.numeric(chaine)
  list <- c(list, list(chaine))
}
full_face_composition$uvec <- list

###########################################################################
# Read data ROMS data
roms_u <- tidync(paste0(path,"data/",filename_u))

# get list of ROMS variables
roms_vars_u <- tidync::hyper_grids(roms_u) %>% # all available grids in the ROMS ncdf
  pluck("grid") %>% # for each grid, pull out all the variables associated with that grid and make a reference table
  purrr::map_df(function(x){
    roms_u %>% tidync::activate(x) %>% tidync::hyper_vars() %>% 
      dplyr::mutate(grd=x)
  })

roms_v <- tidync(paste(path,"data/",filename_v, sep = ""))

# get list of ROMS variables
roms_vars_v <- tidync::hyper_grids(roms_v) %>% # all available grids in the ROMS ncdf
  pluck("grid") %>% # for each grid, pull out all the variables associated with that grid and make a reference table
  purrr::map_df(function(x){
    roms_v %>% tidync::activate(x) %>% tidync::hyper_vars() %>% 
      dplyr::mutate(grd=x)
  })



####################################################################################


u_dim <- roms_vars_u %>% dplyr::filter(name==c("u")) %>% pluck('grd')
u_values <- roms_u %>%
  tidync::activate(u_dim) %>%
  tidync::hyper_tibble() %>%
  dplyr::select(u, longitude, latitude, sigma_layer, time)%>%
  dplyr::rename(
    u = u,
    longitude = longitude,  
    latitude = latitude,  
    roms_layer = sigma_layer,
     time = time)#%>%
  # filter(time >=360)

v_dim <- roms_vars_v %>% dplyr::filter(name==c("v")) %>% pluck('grd')
v_values <- roms_v %>%
  tidync::activate(v_dim) %>%
  tidync::hyper_tibble() %>%
  dplyr::select(v, longitude, latitude, sigma_layer, time)%>%
  dplyr::rename(
    v = v,
    longitude = longitude,   
    latitude = latitude,  
    roms_layer = sigma_layer,
    time = time)#%>%
  # filter(time >=360)


pdt <- sort(unique(u_values$time))
files <- sub("uv_",year,"_", "", list.files("output_",year,"/uv"))
files <- sort(as.numeric(sub(".csv", "", files)))
out <- (1:730)[!1:730 %in% files]
pdt <- out


cores=detectCores()
cl <- cores -1 #not to overload your computer
cl <- 4 #not to overload your computer
registerDoParallel(cl)


foreach(days = pdt) %dopar%{
# for (day in 1:max(pdt)){
  library(raster)#
  library(ncdf4)
  library(tabularaster)
  library(rbgm)
  library(sf) 
  library(tidync)
  library(tidyverse)
  print("Adrien")
sub_u_values <- u_values %>% filter(time == days)
sub_v_values <- v_values %>% filter(time == days)


uv_values2<-  merge(sub_u_values,sub_v_values , by = c("longitude", "latitude", "roms_layer"))
uv_values2<-  merge(uv_values2,uxy2 , by = c("latitude", "longitude"))



######################### Process the value for the first part of the day 


for (i in 1:dim(full_face_composition)[1]){
  extract <- full_face_composition$uvec[i][[1]]
  subset <- uv_values2[uv_values2$roms_layer==full_face_composition$roms_layer[i],]
  subset <- subset[subset$uidx%in%extract,]
  full_face_composition$u[i] <- mean(subset$u)
  full_face_composition$v[i] <- mean(subset$v)
  full_face_composition$lat[i] <- mean(subset$latitude, na.rm = T)
  full_face_composition$lon[i] <- mean(subset$longitude, na.rm = T)
  
}

face_vector <- full_face_composition %>%
  group_by(minz,FaceID) %>%
  mutate(new.u =mean(u, na.rm =T)) %>% 
  mutate(new.v =mean(v, na.rm =T)) %>% 
  select("FaceID",  "cosine" , "sine" ,"atlantis_layer", "face_area" , "maxz" , "minz", "Polygon #", "Face #", "adjacent box" ,"lat","lon","new.u","new.v") %>%
  distinct() %>%
  ungroup()

face_adj_box <- face_vector %>%select("FaceID", "adjacent box", "Polygon #")  %>%
  group_by(FaceID, pair_key = pmin(`Polygon #`, `adjacent box`) + pmax(`Polygon #`, `adjacent box`) * 1000) %>%
  slice(1) %>%
  ungroup() %>%
  select(-pair_key)


face_vector <- face_vector %>%
  select("FaceID",  "cosine" , "sine" ,"atlantis_layer", "face_area" , "maxz" , "minz","lat","lon","new.u","new.v") %>%
  distinct() %>%
  ungroup()

face_vector <- face_vector %>% left_join(face_adj_box,by=c('FaceID')) 

# flux résultant somme u + v est nommé z: z = (u² + v²)^1/2
face_vector$z = (face_vector$new.u^2 + face_vector$new.v^2)^(1/2)

# angle par rapport à u, la latitude, nommé a
# a = atan(U/v) --> to obtain the angle of th vector: atan2(v, u)*180/pi
face_vector$a = atan2(face_vector$new.v, face_vector$new.u)*180/pi



################ Calculate the speed perdilar to the face --> send in y
Polygons_to_dest <- face_vector %>%
  select("FaceID","atlantis_layer","cosine","face_area", "maxz" , "minz", "Polygon #", "adjacent box" ,"new.u","new.v") %>%
  mutate(z = (new.u^2 + new.v^2)^(1/2)) %>%
  mutate(alpha = (atan2(new.v, new.u)*180/pi)) %>%
  mutate(beta = acos(cosine)*180/pi) %>%
  # alpha second: angle between face and vector z
  mutate(alpha.second = 180 - (beta - alpha) ) %>% #TODO: check that this line converse the good direction through the face
  # y: speed perpendicular to the facce
  mutate(y = sin(alpha.second)*z) 

# From flux in m.s-1 -->  *area *60*60*12
Polygons_to_dest$water_flux <- Polygons_to_dest$y *Polygons_to_dest$face_area *60*60*12

# Sum of the flux with the same polygons in origin and desti
Flux_between_polygon <- Polygons_to_dest %>%
  group_by(minz,`Polygon #`, `adjacent box`) %>%
  mutate(water.transfert =sum(water_flux, na.rm =T)) %>% 
  select("atlantis_layer", "maxz" , "minz", "Polygon #", "adjacent box" ,"water.transfert") %>%
  mutate(time = days) %>%
  distinct() %>%
  ungroup() 


###########################################################################
# Hyperdiffusion correction: 

hyperdiffusion_correction_area <- atlantis_sf %>% ungroup %>% select(area,box_id )
hyperdiffusion_correction_area <- as.data.frame(hyperdiffusion_correction_area)[,c(1,2)]
names(hyperdiffusion_correction_area)[2] <- "destination"

Corr.Flux_between_polygon <- Flux_between_polygon %>% 
  mutate(destination = ifelse( water.transfert>0, `adjacent box`, `Polygon #`)) %>%
  left_join(hyperdiffusion_correction_area, by = c("destination")) %>%
  mutate(Corr.water.transfert = water.transfert/((area)))%>%   #Method correct hyperdiffusion
  mutate(atlantis_layer = atlantis_layer-1) %>%
  mutate(Layer_dest = atlantis_layer) %>%
  select(atlantis_layer, `Polygon #`,`adjacent box`, Corr.water.transfert, Layer_dest, time) 


write.csv(Corr.Flux_between_polygon, paste0("output_",year,"_uv/uv_", days,".csv"))
}


stopCluster(cl)




