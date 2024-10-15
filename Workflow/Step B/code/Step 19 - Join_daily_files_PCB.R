
# Set path
if (velma){
  path <- paste0(here(), "/Workflow/Step B/intermediate output archive/output_VELMA_",Nyear,"_PCB")
}else{
  path <- paste0(here(), "/Workflow/Step B/intermediate output archive/output_No_VELMA_",Nyear,"_PCB")
}
setwd(path)


list.file <- sort(list.files(path))


time = seq(0,730*12*60*60-1, 12*60*60) 


Ndt = 1:length(time)
box = 89
layer = 6
N_var = 2

atlantis_input_WC <- array(rep(NA,box*(layer+1)*length(time)), dim = c((layer+1),box,length(time)))
atlantis_input_RPOC <- array(rep(NA,box*(layer+1)*length(time)), dim = c((layer+1),box,length(time)))
atlantis_input_LPOC <- array(rep(NA,box*(layer+1)*length(time)), dim = c((layer+1),box,length(time)))
atlantis_input_DOC <- array(rep(NA,box*(layer+1)*length(time)), dim = c((layer+1),box,length(time)))


liste <- sort(list.file)
for (i in 1:length(list.file)){
  nc <- nc_open(paste0("PCB_Atlantis_",i,".nc"))
  pdt <- ncvar_get(nc, varid = "t")/60/60+1

  atlantis_input_WC <- ncvar_get(nc, varid = "WC")
  # Division by 2 is a rough way to split the PCB between refractory and labile 
  # In can be improve using the [] of each group as a proxy
  atlantis_input_RPOC <- ncvar_get(nc, varid = "POC")/2/5.7 # Redfield ratio
  atlantis_input_LPOC <- ncvar_get(nc, varid = "POC")/2/5.7 # Redfield ratio
  atlantis_input_DOC <- ncvar_get(nc, varid = "DOC")/5.7 # Redfield ratio
  nc_close(nc)
}



###################################################################################
# WC file
###################################################################################
# Define dimensions
z_dim <- ncdim_def("z","layerNum", 1:(layer+1))
b_dim <- ncdim_def("b","boxNum", 0:(box-1))
t_dim <- ncdim_def("t","seconds since 2011-01-01", time, unlim = T)
# Define variables
z_var <- ncvar_def("z", "int", dim = list(z_dim), units = "depthBin", longname = "z")
b_var <- ncvar_def("b", "int", dim = list(b_dim), units = "boxNum", longname = "b")
t_var <- ncvar_def("t", "double", dim = list(t_dim), units = "seconds since 2011-01-01", longname = "t")
WC <- ncvar_def("PCB_habitat", "double", dim = list( z_dim,b_dim, t_dim),
                units = "mg N m-3", missval = 0, longname = "PCB in the WC")


# Create a NetCDF file
if (velma){
  nc_filename <- paste0(here(), "/Workflow/Step B/Final outputs/VELMA/",Nyear,"/pugetsound_SSM_Atlantis_PCB_velma_",Nyear,".nc")
}else{
  nc_filename <- paste0(here(), "/Workflow/Step B/Final outputs/No_VELMA/",Nyear,"/pugetsound_SSM_Atlantis_PCB_","2011.nc")
}

nc <- nc_create(nc_filename, vars = list(WC = WC))

# Put dimensions and variables in the NetCDF file

ncvar_put(nc, z_var, 1:(layer+1))
ncvar_put(nc, b_var, 0:(box-1))
ncvar_put(nc, t_var, time)
ncvar_put(nc, WC, atlantis_input_WC, start = c(1,1,1),count = c( layer+1,box, length(time)))

# Add minimum and maximum values to PCB_habitat variable attributes
ncatt_put(nc, "PCB_habitat", "valid_min", -50)
ncatt_put(nc, "PCB_habitat", "valid_max", 2000)

# Add dt attribute to t variable
ncatt_put(nc, "t", "dt", 43200.0)

# Global attributes
ncatt_put(nc, 0, "title", "PSIMF Atlantis forcing")
ncatt_put(nc, 0, "geometry", "PugetSound_89b_070116.bgm")

# Close the NetCDF file
nc_close(nc)



###################################################################################
# LPON file
###################################################################################
# Define dimensions
t_dim <- ncdim_def("t","seconds since 2011-01-01", time, unlim = T)
z_dim <- ncdim_def("z","layerNum", 1:(layer+1))
b_dim <- ncdim_def("b","boxNum", 0:(box-1))

# Define variables
z_var <- ncvar_def("z", "int", dim = list(z_dim), units = "depthBin", longname = "z")
b_var <- ncvar_def("b", "int", dim = list(b_dim), units = "boxNum", longname = "b")
t_var <- ncvar_def("t", "double", dim = list(t_dim), units = "seconds since 2011-01-01", longname = "t")
LPON <- ncvar_def("PCB_LPON", "double", dim = list( z_dim,b_dim, t_dim),
                units = "mg N m-3", missval = 0, longname = "PCB in LPON")


# Create a NetCDF file
if (velma){
  nc_filename <- paste0(here(), "/Workflow/Step B/Final outputs/VELMA/",Nyear,"/pugetsound_SSM_Atlantis_PCB_LPON_velma_",Nyear,".nc")
}else{
  nc_filename <- paste0(here(), "/Workflow/Step B/Final outputs/No_VELMA/",Nyear,"/pugetsound_SSM_Atlantis_PCB_LPON_","2011.nc")
}



nc <- nc_create(nc_filename, vars = list(LPON = LPON))

# Put dimensions and variables in the NetCDF file

ncvar_put(nc, z_var, 1:(layer+1))
ncvar_put(nc, b_var, 0:(box-1))
ncvar_put(nc, t_var, time)
ncvar_put(nc, LPON, atlantis_input_LPOC, start = c(1,1,1),count = c( layer+1,box, length(time)))

# Add minimum and maximum values to LPON variable attributes
ncatt_put(nc, "PCB_LPON", "valid_min", -1)
ncatt_put(nc, "PCB_LPON", "valid_max", 2000)

# Add dt attribute to t variable
ncatt_put(nc, "t", "dt", 43200.0)

# Global attributes
ncatt_put(nc, 0, "title", "PSIMF Atlantis forcing")
ncatt_put(nc, 0, "geometry", "PugetSound_89b_070116.bgm")





###################################################################################
# RPON file
###################################################################################
# Define dimensions
z_dim <- ncdim_def("z","layerNum", 1:(layer+1))
b_dim <- ncdim_def("b","boxNum", 0:(box-1))
t_dim <- ncdim_def("t","seconds since 2011-01-01", time, unlim = T)
# Define variables
z_var <- ncvar_def("z", "int", dim = list(z_dim), units = "depthBin", longname = "z")
b_var <- ncvar_def("b", "int", dim = list(b_dim), units = "boxNum", longname = "b")
t_var <- ncvar_def("t", "double", dim = list(t_dim), units = "seconds since 2011-01-01", longname = "t")
RPON <- ncvar_def("PCB_RPON", "double", dim = list( z_dim,b_dim, t_dim),
                units = "mg N m-3", missval = 0, longname = "PCB in RPON")


# Create a NetCDF file
if (velma){
  nc_filename <- paste0(here(), "/Workflow/Step B/Final outputs/VELMA/",Nyear,"/pugetsound_SSM_Atlantis_PCB_RPON_velma_",Nyear,".nc")
}else{
  nc_filename <- paste0(here(), "/Workflow/Step B/Final outputs/No_VELMA/",Nyear,"/pugetsound_SSM_Atlantis_PCB_RPON_","2011.nc")
}

nc <- nc_create(nc_filename, vars = list(RPON = RPON))

# Put dimensions and variables in the NetCDF file

ncvar_put(nc, z_var, 1:(layer+1))
ncvar_put(nc, b_var, 0:(box-1))
ncvar_put(nc, t_var, time)
ncvar_put(nc, RPON, atlantis_input_RPOC, start = c(1,1,1),count = c( layer+1,box, length(time)))

# Add minimum and maximum values to RPON variable attributes
ncatt_put(nc, "PCB_RPON", "valid_min", -50)
ncatt_put(nc, "PCB_RPON", "valid_max", 2000)

# Add dt attribute to t variable
ncatt_put(nc, "t", "dt", 43200.0)

# Global attributes
ncatt_put(nc, 0, "title", "PSIMF Atlantis forcing")
ncatt_put(nc, 0, "geometry", "PugetSound_89b_070116.bgm")

# Close the NetCDF file
nc_close(nc)



###################################################################################
# DON file
###################################################################################
# Define dimensions
t_dim <- ncdim_def("t","seconds since 2011-01-01", time, unlim = T)
z_dim <- ncdim_def("z","layerNum", 1:(layer+1))
b_dim <- ncdim_def("b","boxNum", 0:(box-1))

# Define variables
z_var <- ncvar_def("z", "int", dim = list(z_dim), units = "depthBin", longname = "z")
b_var <- ncvar_def("b", "int", dim = list(b_dim), units = "boxNum", longname = "b")
t_var <- ncvar_def("t", "double", dim = list(t_dim), units = "seconds since 2011-01-01", longname = "t")
DON <- ncvar_def("PCB_DON", "double", dim = list( z_dim,b_dim, t_dim),
                units = "mg N m-3", missval = 0, longname = "PCB_DON")


# Create a NetCDF file
if (velma){
  nc_filename <- paste0(here(), "/Workflow/Step B/Final outputs/VELMA/",Nyear,"/pugetsound_SSM_Atlantis_PCB_DON_velma_",Nyear,".nc")
}else{
  nc_filename <- paste0(here(), "/Workflow/Step B/Final outputs/No_VELMA/",Nyear,"/pugetsound_SSM_Atlantis_PCB_DON_","2011.nc")
}



nc <- nc_create(nc_filename, vars = list(DON = DON))

# Put dimensions and variables in the NetCDF file

ncvar_put(nc, z_var, 1:(layer+1))
ncvar_put(nc, b_var, 0:(box-1))
ncvar_put(nc, t_var, time)
ncvar_put(nc, DON, atlantis_input_DOC, start = c(1,1,1),count = c( layer+1,box, length(time)))

# Add minimum and maximum values to DON variable attributes
ncatt_put(nc, "PCB_DON", "valid_min", -1)
ncatt_put(nc, "PCB_DON", "valid_max", 2000)

# Add dt attribute to t variable
ncatt_put(nc, "t", "dt", 43200.0)

# Global attributes
ncatt_put(nc, 0, "title", "PSIMF Atlantis forcing")
ncatt_put(nc, 0, "geometry", "PugetSound_89b_070116.bgm")

# Close the NetCDF file
nc_close(nc)
setwd(here())
