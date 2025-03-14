library(ncdf4)
library(here)

# Variables to change #
year = 2011
velma = F

# Set path
if (velma){
  path <- paste0(here(), "/Step B/intermediate output archive/output_VELMA_",year,"_PON")
}else{
  path <- paste0(here(), "/Step B/intermediate output archive/output_No_VELMA_",year,"_PON")
}
# setwd(path)

list.file <- sort(list.files(path))
time = seq(0,730*12*60*60-1, 12*60*60) 


Ndt = 1:length(time)
box = 89
layer = 6
N_var = 2

atlantis_input_LPON <- array(rep(NA,box*(layer+1)*length(time)), dim = c((layer+1),box,length(time)))
atlantis_input_RPON <- array(rep(NA,box*(layer+1)*length(time)), dim = c((layer+1),box,length(time)))
liste <- sort(list.file)
for (i in 1:length(list.file)){
  nc <- nc_open(paste0(path, "/PON_Atlantis_",i,".nc"))
  pdt <- ncvar_get(nc, varid = "t")/60/60+1
  atlantis_input_LPON[,,i]      <- ncvar_get(nc, varid = "LPON")
  atlantis_input_RPON[,,i]      <- ncvar_get(nc, varid = "RPON")
  nc_close(nc)
}

apply(X = is.na(atlantis_input_RPON),  FUN = sum, MARGIN = c(3))
apply(X = is.na(atlantis_input_LPON),  FUN = sum, MARGIN = c(3))
###################################################################################
# LPON file
###################################################################################
# Define dimensions
z_dim <- ncdim_def("z","layerNum", 1:(layer+1))
b_dim <- ncdim_def("b","boxNum", 0:(box-1))
t_dim <- ncdim_def("t","seconds since 2011-01-01", time, unlim = T)
# Define variables
z_var <- ncvar_def("z", "int", dim = list(z_dim), units = "depthBin", longname = "z")
b_var <- ncvar_def("b", "int", dim = list(b_dim), units = "boxNum", longname = "b")
t_var <- ncvar_def("t", "double", dim = list(t_dim), units = "seconds since 2011-01-01", longname = "t")
LPON <- ncvar_def("Lab_Det_N", "double", dim = list( z_dim,b_dim, t_dim),
                         units = "mg N m-3", missval = 0, longname = "Labile particulate organic nitrogen")


# Create a NetCDF file
if (velma){
  nc_filename <- paste0(here(), "/pugetsound_SSM_Atlantis_LPON_velma_",year,".nc")
}else{
  nc_filename <- paste0(here(), "/pugetsound_SSM_Atlantis_LPON_","2011.nc")
}
nc <- nc_create(nc_filename, vars = list(LPON = LPON))

# Put dimensions and variables in the NetCDF file

ncvar_put(nc, z_var, 1:(layer+1))
ncvar_put(nc, b_var, 0:(box-1))
ncvar_put(nc, t_var, time)
ncvar_put(nc, LPON, atlantis_input_LPON, start = c(1,1,1),count = c( layer+1,box, length(time)))

# Add minimum and maximum values to LPON variable attributes
ncatt_put(nc, "Lab_Det_N", "valid_min", -50)
ncatt_put(nc, "Lab_Det_N", "valid_max", 2000)

# Add dt attribute to t variable
ncatt_put(nc, "t", "dt", 43200.0)

# Global attributes
ncatt_put(nc, 0, "title", "PSIMF Atlantis forcing")
ncatt_put(nc, 0, "geometry", "PugetSound_89b_070116.bgm")

# Close the NetCDF file
nc_close(nc)



###################################################################################
# RPON file
###################################################################################
# Define dimensions
t_dim <- ncdim_def("t","seconds since 2011-01-01", time, unlim = T)
z_dim <- ncdim_def("z","layerNum", 1:(layer+1))
b_dim <- ncdim_def("b","boxNum", 0:(box-1))

# Define variables
z_var <- ncvar_def("z", "int", dim = list(z_dim), units = "depthBin", longname = "z")
b_var <- ncvar_def("b", "int", dim = list(b_dim), units = "boxNum", longname = "b")
t_var <- ncvar_def("t", "double", dim = list(t_dim), units = "seconds since 2011-01-01", longname = "t")
RPON <- ncvar_def("Ref_Det_N", "double", dim = list( z_dim,b_dim, t_dim),
                      units = "mg N m-3", missval = 0, longname = "Refractory particulate organic nitrogen")


# Create a NetCDF file
if (velma){
  nc_filename <- paste0(here(), "/pugetsound_SSM_Atlantis_RPON_velma_",year,".nc")
}else{
  nc_filename <- paste0(here(), "/pugetsound_SSM_Atlantis_RPON_","2011.nc")
}
nc <- nc_create(nc_filename, vars = list(RPON = RPON))

# Put dimensions and variables in the NetCDF file

ncvar_put(nc, z_var, 1:(layer+1))
ncvar_put(nc, b_var, 0:(box-1))
ncvar_put(nc, t_var, time)
ncvar_put(nc, RPON, atlantis_input_RPON, start = c(1,1,1),count = c( layer+1,box, length(time)))

# Add minimum and maximum values to RPON variable attributes
ncatt_put(nc, "Ref_Det_N", "valid_min", -1)
ncatt_put(nc, "Ref_Det_N", "valid_max", 2000)

# Add dt attribute to t variable
ncatt_put(nc, "t", "dt", 43200.0)

# Global attributes
ncatt_put(nc, 0, "title", "PSIMF Atlantis forcing")
ncatt_put(nc, 0, "geometry", "PugetSound_89b_070116.bgm")

# Close the NetCDF file
nc_close(nc)

# meanNA <- function(x){return(mean(x,na.rm = T))}
# plot(1:730, apply(atlantis_input_RPON, FUN = meanNA, MARGIN = 3))
# plot(1:730, apply(atlantis_input_LPON, FUN = meanNA, MARGIN = 3))

