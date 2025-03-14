


###########################################################################
# Path and names definition

path        <- paste0(here(), "/Workflow/Step B/")
input_path <- paste0(here(),"/Workflow/Step A/File_regular_grid/")

# Velma?
if (velma){
  filename <- paste0("VELMA/",Nyear,"/regular_grid_TS_velma_",Nyear,".nc")
  output_path <- paste0(path, "intermediate output archive/output_VELMA_",Nyear,"_TS/")
  
}else{
  filename <- paste0("No_VELMA/",Nyear,"/regular_grid_TS_novelma_",Nyear,".nc")
  output_path <- paste0(path, "intermediate output archive/output_No_VELMA_",Nyear,"_TS/")
}

if (!file.exists(output_path)){dir.create(output_path)}


###########################################################################
# Read data ROMS data
roms <- tidync(paste0(input_path,filename))
box_composition <- read.csv(paste0(path, "code/box_composition.csv"))

###########################################################################

# get list of ROMS variables
roms_vars <- tidync::hyper_grids(roms) %>% # all available grids in the ROMS ncdf
  pluck("grid") %>% # for each grid, pull out all the variables associated with that grid and make a reference table
  purrr::map_df(function(x){
    roms %>% tidync::activate(x) %>% tidync::hyper_vars() %>% 
      dplyr::mutate(grd=x)
  })


############################################################################################
############################################################################################
############################################################################################
step_file <- 1:730 

files <- sub("Physical_var_AtlantisTS_", "", list.files(output_path))
files <- sort(as.numeric(sub(".nc", "", files)))
out <- (1:730)[!1:730 %in% files]
step_file <- out

temperature_dim <- roms_vars %>% dplyr::filter(name==c("temp")) %>% pluck('grd')

# load entire dataset #
variable_before_Atlantis2 <- roms %>%
  tidync::activate(temperature_dim) %>%
  tidync::hyper_tibble(force = TRUE) %>%
  dplyr::select(temp, salinity, longitude, latitude, sigma_layer,time)%>%
  dplyr::rename(
    temp=temp, 
    salinity=salinity, 
    longitude = longitude,  
    latitude = latitude,  
    roms_layer = sigma_layer, time = time)
gc()  # free unsused memory
cores=detectCores()
cl <- cores -1 #not to overload your computer
cl <- 4 #not to overload your computer
registerDoParallel(cl)


foreach(days = step_file) %dopar%{

  variable_before_Atlantis<- variable_before_Atlantis2 %>% filter(time== days)
  variables_polygons <- merge(box_composition, variable_before_Atlantis, by = c("latitude", "longitude", "roms_layer"))

  ###################################################################
  time = sort(unique(variables_polygons$time)) 
  box = 89
  layer = 6
  N_var = 2
  
  atlantis_input_Temp <- array(rep(NA,box*(layer+1)*length(time)), dim = c((layer+1),box,length(time)))
  atlantis_input_salinity <- array(rep(NA,box*(layer+1)*length(time)), dim = c((layer+1),box,length(time)))
  
  for (i in 0:(box-1)){
    for (t in 1:length(time)){
      all.layers_T = rep(NA,6)   # define an empty vector to receive the values of the 6 layers for temperature
      all.layers_sal = rep(NA,6) # define an empty vector to receive the values of the 6 layers for salinity
      
      # Calculate the layer
      for (j in 1:layer){ 
        subset <-variables_polygons %>%
          filter(.bx0 == i, atlantis_layer == j, time == time[t])
        
        
        if (dim(subset)[1] == 0){
          all.layers_T[j] = NA
          all.layers_sal[j] = NA
        }else{
          all.layers_T[j] <- mean(subset$temp, na.rm = T)
          all.layers_sal[j] = mean(subset$salinity, na.rm = T)
        }
      }
      
      keep <- all.layers_T[is.na(all.layers_T)]
      all.layers_T <- c(rev(all.layers_T[!is.na(all.layers_T)]),keep,NA)
      atlantis_input_Temp[,i+1,t] <- all.layers_T
      
      keep <- all.layers_sal[is.na(all.layers_sal)]
      all.layers_sal <- c(rev(all.layers_sal[!is.na(all.layers_sal)]),keep,NA)
      atlantis_input_salinity[,i+1,t] <- all.layers_sal
      
    }
  }
  
  
  ###################################################################################
  # Define nc file
  ###################################################################################
  # Define dimensions
  z_dim <- ncdim_def("z","layerNum", 1:(layer+1))
  b_dim <- ncdim_def("b","boxNum", 0:(box-1))
  t_dim <- ncdim_def("t",paste0("seconds since ",Nyear,"-01-01"), (time-1)*60*60)
  # Define variables
  z_var <- ncvar_def("z", "int", dim = list(z_dim), units = "depthBin", longname = "z")
  b_var <- ncvar_def("b", "int", dim = list(b_dim), units = "boxNum", longname = "b")
  t_var <- ncvar_def("t", "double", dim = list(t_dim), units = paste0("seconds since ",Nyear,"-01-01"), longname = "t")
  temperature <- ncvar_def("temperature", "double", dim = list( z_dim,b_dim, t_dim),
                           units = "°C", missval = NA, longname = "Temperature")
  salinity <- ncvar_def("salinity", "double", dim = list( z_dim,b_dim, t_dim),
                        units = "g.L-1", missval = NA, longname = "Salinity")
  output_filename = paste0("Physical_var_AtlantisTS_", days, ".nc")
  # Create a NetCDF file
  nc_filename <- paste0(output_path, output_filename)
  nc <- nc_create(nc_filename, vars = list(temperature = temperature, salinity = salinity))
  
  # Put dimensions and variables in the NetCDF file
  
  ncvar_put(nc, z_var, 1:(layer+1))
  ncvar_put(nc, b_var, 0:(box-1))
  ncvar_put(nc, t_var, (time-1)*60*60)
  ncvar_put(nc, temperature, atlantis_input_Temp, start = c(1,1,1),count = c( layer+1,box, length(time)))
  ncvar_put(nc, salinity, atlantis_input_salinity, start = c(1,1,1),count = c( layer+1,box, length(time)))
  
  # Add minimum and maximum values to temperature variable attributes
  ncatt_put(nc, "temperature", "valid_min", -50)
  ncatt_put(nc, "temperature", "valid_max", 200)
  
  # Add minimum and maximum values to salinity variable attributes
  ncatt_put(nc, "salinity", "valid_min", 0)
  ncatt_put(nc, "salinity", "valid_max", 2000)
  
  # Add dt attribute to t variable
  ncatt_put(nc, "t", "dt", 43200.0)
  
  # Global attributes
  ncatt_put(nc, 0, "title", "PSIMF Atlantis forcing")
  ncatt_put(nc, 0, "geometry", "PugetSound_89b_070116.bgm")
  ncatt_put(nc, 0, "parameters", "")
  
  # Close the NetCDF file
  nc_close(nc)
  
}
end_time <- Sys.time()
setwd(here())
