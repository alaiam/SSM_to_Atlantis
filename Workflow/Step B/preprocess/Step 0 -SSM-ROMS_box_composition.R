#################################
# Step 0: create the correpondance between regular SSM and Atlantis 
# It will be register in box_compostion.csv
# If this step is already done (you already have box_composition), then go to step 1 (T° and salt)
#################################

start_time <- Sys.time()

# Packages and settings
library(raster)#
library(ncdf4)
library(tabularaster)
library(rbgm)
library(sf) 
library(tidync)
library(tidyverse)

select <- dplyr::select
map <- purrr::map
options(dplyr.summarise.inform=FALSE)


###########################################################################
# Path and names definition
path        <- "C:/Users/Alaia/Desktop/Postdoc/PS/From SSM to Atlantis/"
output_path <- "C:/Users/Alaia/Desktop/Postdoc/PS/From SSM to Atlantis/output/"
setwd("C:/Users/Alaia/Desktop/Postdoc/PS/From SSM to Atlantis/")
filename <- "regular_grid_variables_2011_TS.nc"


###########################################################################
# Read data
# read first ROMS data for depth and grid info - assumes this does not change between time steps and ROMS files
roms <- tidync(paste(path,"data/",filename, sep = ""))


# read Atlantis BGM
atlantis_bgm <- read_bgm(paste(path,"PugetSound_89b_070116.bgm", sep = ""))
# Read SSM data for position of the siglev and the siglay --> they were interpolated in Javier code and then 
# it is not the exact value anymore 
ssm <- tidync(paste(path,"ssm_00117_copy_test.nc", sep = ""))

###########################################################################
# Retrieve data

#Atlantis geometry as an sf shapefile
atlantis_sf <- atlantis_bgm %>% box_sf()


# get list of ssm variables
ssm_vars <- tidync::hyper_grids(ssm) %>% # all available grids in the ssm ncdf
  pluck("grid") %>% # for each grid, pull out all the variables associated with that grid and make a reference table
  purrr::map_df(function(x){
    ssm %>% tidync::activate(x) %>% tidync::hyper_vars() %>% 
      dplyr::mutate(grd=x)
  })

# # extract horizontal exact values

# Siglev = bottom of the cell vertically
siglev_dim <- ssm_vars %>% dplyr::filter(name==c("siglev")) %>% pluck('grd')
siglev_coord <- unname(unlist(ssm %>% tidync::activate(siglev_dim) %>% tidync::hyper_tibble() %>%
                                dplyr::select(siglev)))



###########################################################################

# get list of ROMS variables
roms_vars <- tidync::hyper_grids(roms) %>% # all available grids in the ROMS ncdf
  pluck("grid") %>% # for each grid, pull out all the variables associated with that grid and make a reference table
  purrr::map_df(function(x){
    roms %>% tidync::activate(x) %>% tidync::hyper_vars() %>% 
      dplyr::mutate(grd=x)
  })

# # extract time
time_grd <- roms_vars %>% filter(name=="time_vector") %>% pluck('grd')
roms_time <- roms %>% activate(time_grd) %>% hyper_tibble() %>% pull()

# ROMS grids
# Rho grid
# Horizontal coordinates: on of rho points

lat_dim <- roms_vars %>% dplyr::filter(name==c("latitude")) %>% pluck('grd')
lat_coord <- unname(unlist(roms %>% tidync::activate(lat_dim) %>% tidync::hyper_tibble() %>%
                             dplyr::select(latitude)))

#
lon_dim <- roms_vars %>% dplyr::filter(name==c("longitude")) %>% pluck('grd')
lon_coord <- unname(unlist(roms %>% tidync::activate(lon_dim) %>% tidync::hyper_tibble() %>%
                             dplyr::select(longitude)))

#
siglay_dim <- roms_vars %>% dplyr::filter(name==c("siglay")) %>% pluck('grd')
siglay_coord <- unname(unlist(roms %>% tidync::activate(siglay_dim) %>% tidync::hyper_tibble() %>%
                             dplyr::select(siglay)))
 
# 
grid <- expand.grid(longitude = lon_coord, latitude = lat_coord) 
roms_rho <- as_tibble(grid)
roms_rho <- roms_rho %>% dplyr::mutate(rhoidx=dplyr::row_number())

h_coord <- expand.grid(longitude = lat_coord, latitude = lon_coord, sigma_layer = siglay_coord)
h_coord <- as_tibble(h_coord)%>%  dplyr::rename()

# we create a variable with actually siglay in it (coordinate between cells)
# we use the coordinate of the bottom of the cell
h_coord$siglev <- h_coord$sigma_layer
for (i in 1:length(unique(h_coord$siglev))){
  h_coord$siglev[h_coord$siglev==sort(unique(h_coord$siglev))[i]] <- sort(siglev_coord)[i]
}

############################################################################################
############################################################################################
############################################################################################

# append coordinates to the native lat lon from ROMS
append_xy_coords <- function(lonlatdat,xyproj=atlantis_bgm$extra$projection,lon_col="lon_rho",lat_col="lat_rho"){
  lonlatdat %>% 
    st_as_sf(coords=c(lon_col,lat_col),crs=4326,remove=F) %>%  # convert to spatial object
    st_transform(xyproj) %>%  # convert to Atlantis coords  --> this one is not working, return EMPTY
    mutate(x = st_coordinates(.)[,1],
           y = st_coordinates(.)[,2]) # grab x and y coordinates and add them as attributes
}

rhoxy<- append_xy_coords(roms_rho,lon_col="longitude",lat_col="latitude") %>% mutate(rhoidx=dplyr::row_number())
h_xy<- append_xy_coords(h_coord,lon_col="latitude",lat_col="longitude") %>% mutate(rhoidx=dplyr::row_number())


# In our case, velocity and physical variables measures are done at the same coordinate point 
# join ROMS grids with Atlantis boxes (rho grid) and faces (u and v grids)
# Boxes
boxes_rho_join <- atlantis_sf %>% st_join(rhoxy)

######################### Atlantis depth
# enter the depth breaks in the model
atlantis_z <-sort(unique(atlantis_sf$botz),decreasing = T)[-1] # Obtain the depth of atlantis layer (without surface value)


# small function to build enough depth layers for each box, starting with the defined layers above
# in this function, botz is the bottom depth given for each box, which is available in the .bgm file 
build_Atlantis_depths <- function(botz,lyrs){
  # bottom of each layer, starting from shallowest
  lyrbot<-lyrs
  # layers to use are all those that are shallower than the given botz
  lyr_vec <- lyrbot[lyrbot>botz]
  # the depth of the deepest layer is equal to botz
  lyr_vec <- c(lyr_vec,botz)
  # in Atlantis, each box has the same number of depth layers, but some layers have zero thickness
  # so we have to pad with zeroes to make all boxes have the same number of layers
  nzeroes <- length(lyrs)-length(lyr_vec)
  lyr_vec <- c(lyr_vec,rep(0,nzeroes))
  return(lyr_vec)
}

# function to estimate the overlap between a ROMS cell an Atlantis cell that are overlapping
# the function return the proportion of the ROMS cells that is in the Atlantis cell 
overlap <- function(x){
  top_atlantis <- x[1]
  bottom_atlantis <- x[2]
  top_ROMS <- x[3]
  bottom_ROMS <- x[4]
  if(findInterval(top_ROMS, c(bottom_atlantis,top_atlantis))==1){
    if(findInterval(bottom_ROMS, c(bottom_atlantis,top_atlantis))==1){p =1}
    else{p = (bottom_atlantis-top_ROMS)/(bottom_ROMS-top_ROMS)}
  }else{
    if(findInterval(bottom_ROMS, c(bottom_atlantis,top_atlantis))==1){
      p = (bottom_ROMS-top_atlantis)/(bottom_ROMS-top_ROMS)
    }
    else{p = (bottom_atlantis-top_atlantis)/(bottom_ROMS-top_ROMS)}
  }
  
  return(p)
}

# construct the depth profiles of each Atlantis box

atlantis_depths <- atlantis_bgm$boxes %>% select(.bx0,botz) %>% 
  # apply the function above to create the layers
  mutate(maxz=purrr::map(botz,~build_Atlantis_depths(.,lyrs=atlantis_z))) %>% 
  unnest(cols=c(maxz)) %>% 
  # add a minimum depth for each layer
  group_by(.bx0) %>% 
  mutate(minz=lag(maxz,1,default = 0),atlantis_layer=1:length(atlantis_z)) %>% 
  # add a layer thickness calculation
  mutate(dz=minz-maxz) %>% 
  # "dummy" layers (layers too deep for a given Atlantis box) should have minz and dz=0
  mutate(minz=ifelse(maxz==0,0,minz),dz=ifelse(maxz==0,0,dz)) %>% 
  ungroup() %>% 
  select(.bx0,atlantis_layer,minz,maxz,dz)

atlantis_depths_ROMS <- atlantis_bgm$boxes %>% select(.bx0,botz) %>% 
  # apply the function above to create the layers
  mutate(siglev=purrr::map(botz,~siglev_coord[-1])) %>% 
  unnest(cols=c(siglev)) %>%
  mutate(sigma = siglev)%>% # Keep siglev information in sigma
  mutate(siglev = -siglev * botz) %>% # calculate the depth of each layer
  # add a minimum depth for each layer
  group_by(.bx0) %>% 
  mutate(minz=lag(siglev,1,default = 0),roms_layer=1:length(siglev_coord[-1])) %>% 
  # add a layer thickness calculation
  mutate(dz=minz-siglev) %>% 
  # "dummy" layers (layers too deep for a given Atlantis box) should have minz and dz=0
  mutate(minz=ifelse(siglev==0,0,minz),dz=ifelse(siglev==0,0,dz)) %>% 
  ungroup() %>% 
  select(.bx0,roms_layer,minz,siglev,dz,sigma)

# Keep only polygone of sea 
atlantis_depths <- atlantis_depths %>%
  filter(dz > 0)

# Create 
atlantis_depth_overlap <- left_join(atlantis_depths, atlantis_depths_ROMS, by = ".bx0")%>%
  filter(siglev <= minz.x, minz.y >= maxz) # Keep only ROMS cell that are overlapping Atlantis cells vertically

# temporary data frame for the next step
x <- data.frame(top_atlantis = atlantis_depth_overlap$minz.x, bottom_atlantis = atlantis_depth_overlap$maxz, 
                  top_ROMS = atlantis_depth_overlap$minz.y, bottom_ROMS = atlantis_depth_overlap$siglev)

# Calculate proportion overlap between ROMS cells and atlantis cells vertically 
# with the overlap function
atlantis_depth_overlap$p_overlapping <- apply(x, FUN = overlap, MARGIN = 1)

atlantis_depth_overlap <- atlantis_depth_overlap %>%
  filter(p_overlapping == 0 |p_overlapping >= 0.45) # Keep overlap > 0.5, to avoid duplicate ROMS cells
# 0.45 avoid to have cells without value, keeping a value for box 5 and layer 1

# Merge with horizontal overlap
box_composition <- merge(atlantis_depth_overlap, boxes_rho_join, by=".bx0")
box_composition2 <- box_composition[,c(-11,-12,-13,-14,-15,-16,-17,-18,-19,-20,-21,-22,-23,-26,-27,-28)]
write.csv(box_composition2, "code_PS/box_composition.csv", row.names =F)





