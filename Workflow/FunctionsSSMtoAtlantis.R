# year = 2095
# velma = T

library(ncdf4)
library(here)
library(raster)#
library(ncdf4)
library(tabularaster)
library(rbgm)
library(sf) 
library(tidync)
library(tidyverse)
library(doParallel)
library(reticulate)
library(reshape2)
source("Workflow/Step A/Start_Salish_Sea_env.R")

select <- dplyr::select
map <- purrr::map
options(dplyr.summarise.inform=FALSE)


preprocess <- function(){
  print("Hello, before to use the functions SSMtoAtlantis, var.SSMtoAtlantis, Step A or Step B, be sure to:")
  print("1) Have created your python environment Salish_Sea_env") 
  print("2) Have run the preprocess codes for step B, that create the files box/face _composition, layer_max, and uxy_uv_code")
}


SSMtoAtlantis <- function(year, velma, filename.hyd, filename.wq){
  .year <<- year
  .velma <<- velma
  .filename.hyd <<- filename.hyd
  .filename.wq <<- filename.wq
  
  start.time <- Sys.time()
  write.csv(start.time, "starttime.csv")
  step.time <- c()
   # hyd
  # print("TS")
  # var.SSMtoAtlantis(year = .year, variable = "temperature", velma = .velma, filename = .filename.hyd)
  # print("UVW")
  # step.time <- c(step.time,Sys.time())
  # write.csv(step.time, "steptime.csv")
  # 
  # var.SSMtoAtlantis(year = .year, variable = "U", velma = .velma, filename = .filename.hyd)
  # step.time <- c(step.time,Sys.time())
  # write.csv(step.time, "steptime.csv")
  # 
  # #wq
  # print("Z")
  # var.SSMtoAtlantis(year = .year, variable = "SZ", velma = .velma, filename = .filename.wq)
  # step.time <- c(step.time,Sys.time())
  # write.csv(step.time, "steptime.csv")
  # 
  # print("P")
  # var.SSMtoAtlantis(year = .year, variable = "SP", velma = .velma, filename = .filename.wq)
  # step.time <- c(step.time,Sys.time())
  # write.csv(step.time, "steptime.csv")
  # 
  # print("N")
  # var.SSMtoAtlantis(year = .year, variable = "NO3", velma = .velma, filename = .filename.wq)
  # step.time <- c(step.time,Sys.time())
  # write.csv(step.time, "steptime.csv")
  # 
  # print("DON")
  # var.SSMtoAtlantis(year = .year, variable = "RDON", velma = .velma, filename = .filename.wq)
  # step.time <- c(step.time,Sys.time())
  # write.csv(step.time, "steptime.csv")
  
  print("O2")
  var.SSMtoAtlantis(year = .year, variable = "Oxygen", velma = .velma, filename = .filename.wq)
  step.time <- c(step.time,Sys.time())
  write.csv(step.time, "steptime.csv")
  
  print("PON")
  var.SSMtoAtlantis(year = .year, variable = "LPON", velma = .velma, filename = .filename.wq)
  step.time <- c(step.time,Sys.time())
  write.csv(step.time, "steptime.csv")
  
  end.time <- Sys.time()
  write.csv(end.time, "endtime.csv")
  
}
var.SSMtoAtlantis <- function(year, variable, velma, filename){
  
  list.var = c("salinity","temperature",
               "U", "V", "W",
               "NO3", "NH4", 
               "SZ", "LZ", "MZ", "SP", "LP",
               "Oxygen","LPON","RPON",
               "RDON")
  
  if(!variable%in%list.var) stop("The variable is not in SSM, please try:
salinity, temperature, U, V, W, NO3, NH4, SZ, LZ, MZ, SP, LP, Oxygen, LPON, RPON, RDON")
  # print("var1")
  StepA(year = year, variable = variable, velma = velma, filename = filename)
  print("var2")
  StepB(year = year, variable = variable, velma = velma)
  print("var3")
  Joindailyfile(year = year, variable = variable, velma = velma)

}

# Step A 
StepA <- function(year, variable, velma, filename){
  file_name_input <<- filename
  print(filename)
  list.var = c("salinity","temperature",
               "U", "V", "W",
               "NO3", "NH4", 
               "SZ", "LZ", "MZ", "SP", "LP",
               "Oxygen","LPON","RPON",
               "RDON")
  
  if(!variable%in%list.var) stop("The variable is not in SSM, please try:
salinity, temperature, U, V, W, NO3, NH4, SZ, LZ, MZ, SP, LP, Oxygen, LPON, RPON, RDON")
  if(velma == F) stop("Option not available yet --> TODO --> add the condition to have on the output name file")
 
  if (velma){
    path =   paste0('/home/atlantis/amps_hydrodynamics/Workflow/Step A/File_regular_grid/VELMA/', year, "/")
  }else{
    path =   paste0('/home/atlantis/amps_hydrodynamics/Workflow/Step A/File_regular_grid/No_VELMA/', year, "/")
  }
  
  if (!file.exists(path)){dir.create(path)}
  
  if (variable == "salinity"|| variable == "temperature"){
    file_name_output <<- paste0(path, 'regular_grid_TS_velma_',year,'.nc') 
    print(file_name_output)
    py_run_file("Workflow/Step A/StepA_TS.py")
  }
  
  if (variable == "U"|| variable == "V"|| variable == "W"){
    file_name_output <<- paste0(path, 'regular_grid_UVW_velma_',year,'.nc') 
    print(file_name_output)
    py_run_file("Workflow/Step A/StepA_UVW.py")
  }
  
  if (variable == "NO3"|| variable == "NH4"){
    file_name_output <<- paste0(path, 'regular_grid_N_velma_',year,'.nc') 
    print(file_name_output)
    py_run_file("Workflow/Step A/StepA_N.py")
  }
  
  if (variable == "SZ"|| variable == "LZ"|| variable == "MZ"){
    file_name_output <<- paste0(path, 'regular_grid_Z_velma_',year,'.nc') 
    print(file_name_output)
    py_run_file("Workflow/Step A/StepA_Z.py")
    }
  
  if (variable == "SP"|| variable == "LP"){
    file_name_output <<- paste0(path, 'regular_grid_B_velma_',year,'.nc') 
    print(file_name_output)
    py_run_file("Workflow/Step A/StepA_B.py")
    }
  
  if (variable == "Oxygen"){
    file_name_output <<- paste0(path, 'regular_grid_Oxygen_velma_',year,'.nc') 
    print(file_name_output)
    py_run_file("Workflow/Step A/StepA_O2.py")
    }
  
  if (variable == "LPON"|| variable == "RPON"){
    file_name_output <<- paste0(path, 'regular_grid_PON_velma_',year,'.nc') 
    print(file_name_output)
    py_run_file("Workflow/Step A/StepA_PON.py")
    }
  
  if (variable == "RDON"){
    file_name_output <<- paste0(path, 'regular_grid_DON_velma_',year,'.nc') 
    print(file_name_output)
    py_run_file("Workflow/Step A/StepA_DON.py")
    }
}
# Step B
StepB <- function(year, variable, velma){
  Nyear <<- year
  velma <<- velma
  variable <<- variable
  list.var = c("salinity","temperature",
               "U", "V", "W",
               "NO3", "NH4", 
               "SZ", "LZ", "MZ", "SP", "LP",
               "Oxygen","LPON","RPON",
               "RDON")
  
  if(!variable%in%list.var) stop("The variable is not in SSM, please try:
salinity, temperature, U, V, W, NO3, NH4, SZ, LZ, MZ, SP, LP, Oxygen, LPON, RPON, RDON. 
For a pdf sum up, try all ")
  

  if (variable == "salinity"|| variable =="temperature")   {source("Workflow/Step B/code/Step 1 - SSM-ROMS_to_Atlantis_T_S.R")}
  if (variable == "U"|| variable =="V"|| variable =="W")   {source("Workflow/Step B/code/Step 2 - SSM-ROMS_to_Atlantis_w.R")
                                                            source("Workflow/Step B/code/Step 3 - SSM-ROMS_to_Atlantis_uv.R")}
  if (variable == "SP"|| variable =="LP")                  {source("Workflow/Step B/code/Step 6 - SSM-ROMS_to_Atlantis_B.R")}
  if (variable == "NO3"|| variable =="NH4")                {source("Workflow/Step B/code/Step 7 - SSM-ROMS_to_Atlantis_N.R")}
  if (variable == "SZ"|| variable =="LZ"|| variable =="MZ"){source("Workflow/Step B/code/Step 8 - SSM-ROMS_to_Atlantis_Z.R")}
  if (variable == "Oxygen")                                {source("Workflow/Step B/code/Step 10 - SSM-ROMS_to_Atlantis_O2.R")}
  if (variable == "LPON"|| variable =="RPON")              {source("Workflow/Step B/code/Step 12 - SSM-ROMS_to_Atlantis_PON.R")}
  if (variable == "RDON"|| variable =="DON")               {source("Workflow/Step B/code/Step 14 - SSM-ROMS_to_Atlantis_DON.R")}
  
  
}
# Join daily files
Joindailyfile <- function(year, variable, velma){
  
  Nyear <<- year
  velma <<- velma
  variable <<- variable 
  
  list.var = c("salinity","temperature",
               "U", "V", "W",
               "NO3", "NH4", 
               "SZ", "LZ", "MZ", "SP", "LP",
               "Oxygen","LPON","RPON",
               "RDON")
  
  
  if(!variable%in%list.var) stop("The variable is not in SSM, please try:
salinity, temperature, U, V, W, NO3, NH4, SZ, LZ, MZ, SP, LP, Oxygen, LPON, RPON, RDON")
  
  list.folder.names = c("TS","TS",
                        "uv", "uv", "uv",
                        "N", "N", 
                        "Z", "Z", "Z", "B", "B",
                        "O2","PON","PON",
                        "DON")
  
  pos <- match(variable, list.var)
  if (velma){
    output_path = paste0("Workflow/Step B/intermediate output archive/output_VELMA_",year, "_", list.folder.names[pos])
    final_path = paste0("Workflow/Step B/Final outputs/VELMA/",year)
  }else{
    output_path = paste0("Workflow/Step B/intermediate output archive/output_No_VELMA_",year, "_", list.folder.names[pos])
    final_path = paste0("Workflow/Step B/Final outputs/No_VELMA/",year)
    
  }
  print(output_path)
  if (!file.exists(final_path)){dir.create(final_path)}
  

  if (length(list.files(output_path))<730){
    print("toto")
    StepB(year, variable, velma)
    Joindailyfile(year, variable, velma)
  }else{
    if (variable == "salinity"|| variable =="temperature"){
      source("Workflow/Step B/code/Step 4 - Join TS daily files.R")
    }
    if (variable == "U"|| variable =="V"|| variable =="W"){
      output_path = paste0("Workflow/Step B/intermediate output archive/output_VELMA_",year, "_ww")
      if (length(list.files(output_path))<730){
        print("totou")
        StepB(year, variable, velma)
        Joindailyfile(year, variable, velma)
      }else{
      source("Workflow/Step B/code/Step 5 - Join_daily_files_uv.R")}
    }
    if (variable == "NO3"|| variable =="NH4"){
      print("oy")
      source("Workflow/Step B/code/Step 9 - Join_daily_files_N.R")
    }
    if (variable == "SZ"|| variable =="LZ"|| variable =="MZ"){
      source("Workflow/Step B/code/Step 9 - Join_daily_files_Z.R")
    }
    if (variable == "SP"|| variable =="LP"){
      source("Workflow/Step B/code/Step 9 - Join_daily_files_B.R")
    }
    if (variable == "Oxygen"){
      source("Workflow/Step B/code/Step 11 - Join_daily_files_O2.R")
    }
    if (variable == "LPON"|| variable =="RPON"){
      source("Workflow/Step B/code/Step 13 - Join_daily_files_PON.R")
    }
    if (variable == "RDON"|| variable =="DON"){
      source("Workflow/Step B/code/Step 15 - Join_daily_files_DON.R")
    }
    }}
  


