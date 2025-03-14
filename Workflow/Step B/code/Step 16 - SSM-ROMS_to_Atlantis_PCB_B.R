


###########################################################################
# Path and names definition

path        <- paste0(here(), "/Workflow/Step B/")
input_path <- paste0(here(),"/Workflow/Step A/File_regular_grid/")

# Velma?
Velma = F
Nyear = 2011
if (Velma){
  filename <- paste0("VELMA/",Nyear,"/regular_grid_PCB_B_velma_",Nyear,".nc")
  output_path <- paste0(path, "intermediate output archive/output_VELMA_",Nyear,"_PCB_B/")
  
}else{
  filename <- paste0("No_VELMA/",Nyear,"/regular_grid_PCB_B_novelma_",Nyear,".nc")
  output_path <- paste0(path, "intermediate output archive/output_No_VELMA_",Nyear,"_PCB_B/")
}

if (!file.exists(output_path)){dir.create(output_path)}


###########################################################################
# Read data ROMS data
roms <- tidync(paste(input_path,filename, sep = ""))
box_composition <- read.csv("Step B/code/box_composition.csv")

###########################################################################

# get list of ROMS variables
roms_vars <- tidync::hyper_grids(roms) %>% # all available grids in the ROMS ncdf
  pluck("grid") %>% # for each grid, pull out all the variables associated with that grid and make a reference table
  purrr::map_df(function(x){
    roms %>% tidync::activate(x) %>% tidync::hyper_vars() %>% 
      dplyr::mutate(grd=x)
  })

#### 
atlantis_bgm <- read_bgm(paste(path,"PugetSound_89b_070116.bgm", sep = ""))
atlantis_sf <- atlantis_bgm %>% box_sf()
area  <- (atlantis_sf %>% ungroup %>%
            st_drop_geometry()%>% select(area,box_id ))

layer_thickness <- c(5,20,25,50,50,200)

############################################################################################
############################################################################################
############################################################################################
step_file <- 1:730 #Days to divide the total files

files <- sub("PCB_B_Atlantis_", "", list.files(output_path))
files <- sort(as.numeric(sub(".nc", "", files)))
out <- (1:730)[!1:730 %in% files]
step_file <- out

PCB1_dim <- roms_vars %>% dplyr::filter(name==c("PCB1")) %>% pluck('grd')


variable_before_Atlantis2 <- roms %>%
  tidync::activate(PCB1_dim) %>%
  tidync::hyper_tibble(force = TRUE) %>%
  dplyr::select(PCB1, PCB2, longitude, latitude, sigma_layer,time)%>%
  dplyr::rename(
    PCB1=PCB1, 
    PCB2=PCB2, 
    longitude = longitude,  
    latitude = latitude,  
    roms_layer = sigma_layer, time = time)



gc() #free unused memory before parallelization
cores=detectCores()
cl <- cores -1 #not to overload your computer
cl <- 4 #not to overload your computer
registerDoParallel(cl)

foreach(days = step_file) %dopar%{
  # for (days in 1:length(step_file)){
  
  
  
  
  variable_before_Atlantis<- variable_before_Atlantis2 %>% filter(time== days)
  
  
  
  variables_polygons <- merge(box_composition, variable_before_Atlantis, by = c("latitude", "longitude", "roms_layer"))
  
  ###################################################################
  time = sort(unique(variables_polygons$time)) 
  box = 89
  layer = 6
  N_var = 2
  
  atlantis_input_PCB1 <- array(rep(NA,box*(layer+1)*length(time)), dim = c((layer+1),box,length(time)))
  atlantis_input_PCB2 <- array(rep(NA,box*(layer+1)*length(time)), dim = c((layer+1),box,length(time)))
  
  for (i in 0:(box-1)){
    for (t in 1:length(time)){
      all.layers_PCB1 = rep(NA,6)   # define an empty vector to receive the values of the 6 layers for PCB1
      all.layers_PCB2 = rep(NA,6) # define an empty vector to receive the values of the 6 layers for PCB2
      # Calculate the layer
      for (j in 1:layer){ 
        subset <-variables_polygons %>%
          filter(.bx0 == i, atlantis_layer == j, time == time[t])
        
        
        if (dim(subset)[1] == 0){
          all.layers_PCB1[j] = NA
          all.layers_PCB2[j] = NA
        }else{
          all.layers_PCB1[j] <- (mean(subset$PCB1, na.rm = T)*1000)[[1]] #g to mg
          all.layers_PCB2[j] <- (mean(subset$PCB2, na.rm = T)*1000)[[1]] #g to mg
          # all.layers_PCB1[j] <- (mean(subset$PCB1, na.rm = T)*area[i+1,1]*layer_thickness[j]*1000)[[1]] #redfield ratio
          # all.layers_PCB2[j] <- (mean(subset$PCB2, na.rm = T)*area[i+1,1]*layer_thickness[j]*1000)[[1]] #redfield ratio
        }
      }
      
      keep <- all.layers_PCB1[is.na(all.layers_PCB1)]
      all.layers_PCB1 <- c(rev(all.layers_PCB1[!is.na(all.layers_PCB1)]),keep,NA)
      atlantis_input_PCB1[,i+1,t] <- all.layers_PCB1
      
      keep <- all.layers_PCB2[is.na(all.layers_PCB2)]
      all.layers_PCB2 <- c(rev(all.layers_PCB2[!is.na(all.layers_PCB2)]),keep,NA)
      atlantis_input_PCB2[,i+1,t] <- all.layers_PCB2
      
    }
  }
  
  
  ###################################################################################
  # Define nc file
  ###################################################################################
  # Define dimensions
  z_dim <- ncdim_def("z","layerNum", 1:(layer+1))
  b_dim <- ncdim_def("b","boxNum", 0:(box-1))
  t_dim <- ncdim_def("t","seconds since 2011-01-01", (time-1)*60*60)
  # Define variables
  z_var <- ncvar_def("z", "int", dim = list(z_dim), units = "depthBin", longname = "z")
  b_var <- ncvar_def("b", "int", dim = list(b_dim), units = "boxNum", longname = "b")
  t_var <- ncvar_def("t", "double", dim = list(t_dim), units = "seconds since 2011-01-01", longname = "t")
  PCB1 <- ncvar_def("PCB1", "double", dim = list( z_dim,b_dim, t_dim),
                   units = "mgN", missval = NA, longname = "PCB1")
  PCB2 <- ncvar_def("PCB2", "double", dim = list( z_dim,b_dim, t_dim),
                   units = "g.L-1", missval = NA, longname = "PCB2")
  output_filename = paste0("PCB_B__Atlantis_", days, ".nc")
  # Create a NetCDF file
  nc_filename <- paste0(output_path, output_filename)
  nc <- nc_create(nc_filename, vars = list(PCB1 = PCB1, PCB2 = PCB2))
  
  # Put dimensions and variables in the NetCDF file
  
  ncvar_put(nc, z_var, 1:(layer+1))
  ncvar_put(nc, b_var, 0:(box-1))
  ncvar_put(nc, t_var, (time-1)*60*60)
  ncvar_put(nc, PCB1, atlantis_input_PCB1, start = c(1,1,1),count = c( layer+1,box, length(time)))
  ncvar_put(nc, PCB2, atlantis_input_PCB2, start = c(1,1,1),count = c( layer+1,box, length(time)))
  
  # Add minimum and maximum values to PCB1 variable attributes
  ncatt_put(nc, "PCB1", "valid_min", -50)
  ncatt_put(nc, "PCB1", "valid_max", 200)
  
  # Add minimum and maximum values to PCB2 variable attributes
  ncatt_put(nc, "PCB2", "valid_min", 0)
  ncatt_put(nc, "PCB2", "valid_max", 2000)
  
  # Add dt attribute to t variable
  ncatt_put(nc, "t", "dt", 43200.0)
  
  # Global attributes
  ncatt_put(nc, 0, "title", "PSIMF Atlantis forcing")
  ncatt_put(nc, 0, "geometry", "PugetSound_89b_070116.bgm")
  ncatt_put(nc, 0, "parameters", "")
  
  # Close the NetCDF file
  nc_close(nc)
  
}


