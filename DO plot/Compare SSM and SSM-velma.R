# Charger la bibliothèque nécessaire
library(ncdf4)
library(rbgm)
library(tidyverse)
library(here)

# Définir le répertoire de travail
setwd(here())


plot_SSM_physics <- function(varname, unit = "", varid = varname, title = title){
  
  if (varname == "salinity"){
    ylab = "Salinity (g.L-1)"
    file_ROMS <- "data/pugetsound_roms_salt.nc"
  }
  if (varname == "temperature"){
    ylab = "Temperature (°C)"
    file_ROMS <- "data/pugetsound_roms_temp.nc" 
  }
  
  file_SSM <- paste0("Workflow/Step B/Final outputs/No_VELMA/2011/pugetsound_SSM_Atlantis_",varname,"_novelma_2011.nc")
  file_SSM_velma <- paste0("Workflow/Step B/Final outputs/VELMA/2011/pugetsound_SSM_Atlantis_",varname,"_velma_2011.nc")
  
  # volume data for weighted average
  atlantis_bgm <- read_bgm("Workflow/Step B/PugetSound_89b_070116.bgm")
  atlantis_sf <- atlantis_bgm %>% box_sf()
  area <- atlantis_sf$area
  
  
  # Obtenir les données de la variable
  nc <- nc_open(file_SSM)
  var <- ncvar_get(nc, varid = varname)
  
  nc <- nc_open(file_SSM_velma)
  var_velma <- ncvar_get(nc, varid = varname)
  
  nc <- nc_open(file_ROMS)
  var_roms <- ncvar_get(nc, varid = varname)
  
  
  # Transforamtion des données #
  tickness <-  sort((sort(unique(atlantis_sf$botz))[-7] - sort(unique(atlantis_sf$botz))[-1])*-1)
  
  st <- apply(X = var[,,1], FUN = is.na, MARGIN = c(1,2))*-1+1
  layer <-  st
  for (i in 1:89){
    a <- layer[,i]
    a[1:sum(a)] <- sum(a):1
    layer[,i] <- a
  }
  volume <- st
  for (i in 1:6){
    for(j in 1:89){
      if (layer[i,j]>0){
        volume[i,j] <- volume[i,j]*area[j]*tickness[layer[i,j]]
      }
    }
  }
  
  mean_var_overall <- mean_var_1 <- mean_var_2 <- mean_var_3 <- mean_var_4 <- mean_var_5 <- mean_var_6 <- c()
  min_var_overall <- min_var_1 <- min_var_2 <- min_var_3 <- min_var_4 <- min_var_5 <- min_var_6 <- c()
  max_var_overall <- max_var_1 <- max_var_2 <- max_var_3 <- max_var_4 <- max_var_5 <- max_var_6 <- c()
  
  mean_var_velma_overall <- mean_var_velma_1 <- mean_var_velma_2 <- mean_var_velma_3 <- mean_var_velma_4 <- mean_var_velma_5 <- mean_var_velma_6 <- c()
  min_var_velma_overall <- min_var_velma_1 <- min_var_velma_2 <- min_var_velma_3 <- min_var_velma_4 <- min_var_velma_5 <- min_var_velma_6 <- c()
  max_var_velma_overall <- max_var_velma_1 <- max_var_velma_2 <- max_var_velma_3 <- max_var_velma_4 <- max_var_velma_5 <- max_var_velma_6 <- c()
  
  mean_var_roms_overall <- mean_var_roms_1 <- mean_var_roms_2 <- mean_var_roms_3 <- mean_var_roms_4 <- mean_var_roms_5 <- mean_var_roms_6 <- c()
  min_var_roms_overall <- min_var_roms_1 <- min_var_roms_2 <- min_var_roms_3 <- min_var_roms_4 <- min_var_roms_5 <- min_var_roms_6 <- c()
  max_var_roms_overall <- max_var_roms_1 <- max_var_roms_2 <- max_var_roms_3 <- max_var_roms_4 <- max_var_roms_5 <- max_var_roms_6 <- c()
  
  for (i in 1:730){
    mean_var_overall <- c(mean_var_overall, sum(volume*var[,,i]/sum(volume)))
    mean_var_1 <- c(mean_var_1, sum(volume[layer == 1]*var[,,i][layer == 1]/sum(volume[layer == 1])))
    mean_var_2 <- c(mean_var_2, sum(volume[layer == 2]*var[,,i][layer == 2]/sum(volume[layer == 2])))
    mean_var_3 <- c(mean_var_3, sum(volume[layer == 3]*var[,,i][layer == 3]/sum(volume[layer == 3])))
    mean_var_4 <- c(mean_var_4, sum(volume[layer == 4]*var[,,i][layer == 4]/sum(volume[layer == 4])))
    mean_var_5 <- c(mean_var_5, sum(volume[layer == 5]*var[,,i][layer == 5]/sum(volume[layer == 5])))
    mean_var_6 <- c(mean_var_6, sum(volume[layer == 6]*var[,,i][layer == 6]/sum(volume[layer == 6])))
    
    max_var_overall <- c(max_var_overall, max(var[,,i]))
    max_var_1 <- c(max_var_1, max(var[,,i][layer == 1]))
    max_var_2 <- c(max_var_2,max(var[,,i][layer == 2]))
    max_var_3 <- c(max_var_3,max(var[,,i][layer == 3]))
    max_var_4 <- c(max_var_4,max(var[,,i][layer == 4]))
    max_var_5 <- c(max_var_5,max(var[,,i][layer == 5]))
    max_var_6 <- c(max_var_6,max(var[,,i][layer == 6]))
    
    min_var_overall <- c(min_var_overall, min(var[,,i]))
    min_var_1 <- c(min_var_1, min(var[,,i][layer == 1]))
    min_var_2 <- c(min_var_2,min(var[,,i][layer == 2]))
    min_var_3 <- c(min_var_3,min(var[,,i][layer == 3]))
    min_var_4 <- c(min_var_4,min(var[,,i][layer == 4]))
    min_var_5 <- c(min_var_5,min(var[,,i][layer == 5]))
    min_var_6 <- c(min_var_6,min(var[,,i][layer == 6]))
    
    mean_var_velma_overall <- c(mean_var_velma_overall, sum(volume*var_velma[,,i]/sum(volume)))
    mean_var_velma_1 <- c(mean_var_velma_1, sum(volume[layer == 1]*var_velma[,,i][layer == 1]/sum(volume[layer == 1])))
    mean_var_velma_2 <- c(mean_var_velma_2, sum(volume[layer == 2]*var_velma[,,i][layer == 2]/sum(volume[layer == 2])))
    mean_var_velma_3 <- c(mean_var_velma_3, sum(volume[layer == 3]*var_velma[,,i][layer == 3]/sum(volume[layer == 3])))
    mean_var_velma_4 <- c(mean_var_velma_4, sum(volume[layer == 4]*var_velma[,,i][layer == 4]/sum(volume[layer == 4])))
    mean_var_velma_5 <- c(mean_var_velma_5, sum(volume[layer == 5]*var_velma[,,i][layer == 5]/sum(volume[layer == 5])))
    mean_var_velma_6 <- c(mean_var_velma_6, sum(volume[layer == 6]*var_velma[,,i][layer == 6]/sum(volume[layer == 6])))
    
    max_var_velma_overall <- c(max_var_velma_overall, max(var_velma[,,i]))
    max_var_velma_1 <- c(max_var_velma_1, max(var_velma[,,i][layer == 1]))
    max_var_velma_2 <- c(max_var_velma_2,max(var_velma[,,i][layer == 2]))
    max_var_velma_3 <- c(max_var_velma_3,max(var_velma[,,i][layer == 3]))
    max_var_velma_4 <- c(max_var_velma_4,max(var_velma[,,i][layer == 4]))
    max_var_velma_5 <- c(max_var_velma_5,max(var_velma[,,i][layer == 5]))
    max_var_velma_6 <- c(max_var_velma_6,max(var_velma[,,i][layer == 6]))
    
    min_var_velma_overall <- c(min_var_velma_overall, min(var_velma[,,i]))
    min_var_velma_1 <- c(min_var_velma_1, min(var_velma[,,i][layer == 1]))
    min_var_velma_2 <- c(min_var_velma_2,min(var_velma[,,i][layer == 2]))
    min_var_velma_3 <- c(min_var_velma_3,min(var_velma[,,i][layer == 3]))
    min_var_velma_4 <- c(min_var_velma_4,min(var_velma[,,i][layer == 4]))
    min_var_velma_5 <- c(min_var_velma_5,min(var_velma[,,i][layer == 5]))
    min_var_velma_6 <- c(min_var_velma_6,min(var_velma[,,i][layer == 6]))
    
    mean_var_roms_overall <- c(mean_var_roms_overall, sum(volume*var_roms[,,i]/sum(volume)))
    mean_var_roms_1 <- c(mean_var_roms_1, sum(volume[layer == 1]*var_roms[,,i][layer == 1]/sum(volume[layer == 1])))
    mean_var_roms_2 <- c(mean_var_roms_2, sum(volume[layer == 2]*var_roms[,,i][layer == 2]/sum(volume[layer == 2])))
    mean_var_roms_3 <- c(mean_var_roms_3, sum(volume[layer == 3]*var_roms[,,i][layer == 3]/sum(volume[layer == 3])))
    mean_var_roms_4 <- c(mean_var_roms_4, sum(volume[layer == 4]*var_roms[,,i][layer == 4]/sum(volume[layer == 4])))
    mean_var_roms_5 <- c(mean_var_roms_5, sum(volume[layer == 5]*var_roms[,,i][layer == 5]/sum(volume[layer == 5])))
    mean_var_roms_6 <- c(mean_var_roms_6, sum(volume[layer == 6]*var_roms[,,i][layer == 6]/sum(volume[layer == 6])))
    
    max_var_roms_overall <- c(max_var_roms_overall, max(var_roms[,,i]))
    max_var_roms_1 <- c(max_var_roms_1, max(var_roms[,,i][layer == 1]))
    max_var_roms_2 <- c(max_var_roms_2,max(var_roms[,,i][layer == 2]))
    max_var_roms_3 <- c(max_var_roms_3,max(var_roms[,,i][layer == 3]))
    max_var_roms_4 <- c(max_var_roms_4,max(var_roms[,,i][layer == 4]))
    max_var_roms_5 <- c(max_var_roms_5,max(var_roms[,,i][layer == 5]))
    max_var_roms_6 <- c(max_var_roms_6,max(var_roms[,,i][layer == 6]))
    
    min_var_roms_overall <- c(min_var_roms_overall, min(var_roms[,,i]))
    min_var_roms_1 <- c(min_var_roms_1, min(var_roms[,,i][layer == 1]))
    min_var_roms_2 <- c(min_var_roms_2,min(var_roms[,,i][layer == 2]))
    min_var_roms_3 <- c(min_var_roms_3,min(var_roms[,,i][layer == 3]))
    min_var_roms_4 <- c(min_var_roms_4,min(var_roms[,,i][layer == 4]))
    min_var_roms_5 <- c(min_var_roms_5,min(var_roms[,,i][layer == 5]))
    min_var_roms_6 <- c(min_var_roms_6,min(var_roms[,,i][layer == 6]))
    
  }
  
  #######################################################
  ######################### MEAN ########################
  #######################################################
  
  
  # pdf(width=9,height = 4.5)
  # Configuration de la disposition des graphiques et ajustement des marges globales
  par(mfrow = c(3, 2), oma = c(2.5, 2.5, 2, 1), mar = c(2, 2, 2, 1))
  
  # Définir les limites de l'axe y
  abs_min <- min(c(min_var_roms_1, min_var_1, min_var_velma_1, min_var_roms_5, min_var_5, min_var_velma_5))
  abs_max <- max(c(max_var_roms_1, max_var_1, max_var_velma_1, max_var_roms_5, max_var_5, max_var_velma_5))*1.2
  
  # Surface layer
  plot(1:730/2, mean_var_roms_1, type = "l", col = "#FF4500", lwd = 1.5, xlab = "",
       ylab = "", main = "Layer 1 (0-5m)", ylim = c(abs_min, abs_max), axes = FALSE)
  axis(1, at = c(-100,0,100,200,300,400))
  axis(2, at = seq(round(abs_min)/2,round(abs_max)*2, (round(abs_max)-round(abs_min))/8))
  # box(col = )
  lines(1:730/2, mean_var_1, col = "#1E3A8A", lwd = 1.5)
  lines(1:730/2, mean_var_velma_1, col = "#ADFF2F", lwd = 1.5)
  grid(nx = NULL, ny = NULL, col = "lightgray", lty = "dotted")
  legend("topleft", legend = c("ROMS", "SSM", "SSM-VELMA"), col = c("#FF4500", "#1E3A8A", "#ADFF2F"), lty = 1, lwd = 1.5, bty = "n")
  legend("topright", legend = c("AVERAGE"), col = c( "white"), lty = 1, lwd = 1.5, bty = "n")
  
  # Depth layer
  plot(1:730/2, mean_var_roms_5, type = "l", col = "#FF4500", lwd = 1.5, xlab = "",
       ylab = "", main = "Layer 5 (100-150 m)", ylim = c(abs_min, abs_max), axes = FALSE)
  axis(1, at = c(-100,0,100,200,300,400))
  axis(2, at = seq(round(abs_min)/2,round(abs_max)*2, (round(abs_max)-round(abs_min))/8), labels = FALSE)
  
  # box(col = )
  lines(1:730/2, mean_var_5, col = "#1E3A8A", lwd = 1.5)
  lines(1:730/2, mean_var_velma_5, col = "#ADFF2F", lwd = 1.5)
  grid(nx = NULL, ny = NULL, col = "lightgray", lty = "dotted")
  legend("topleft", legend = c("ROMS", "SSM", "SSM-VELMA"), col = c("#FF4500", "#1E3A8A", "#ADFF2F"), lty = 1, lwd = 1.5, bty = "n")
  
  # Ajouter les étiquettes des axes globaux
  mtext("Days", side = 1, outer = TRUE, line = 1, col = , cex = 1.2)
  mtext(ylab, side = 2, outer = TRUE, line = 1, col = , cex = 1.2)
  mtext(title, side = 3, outer = TRUE, line = 0, col = , cex = 1.2)
  legend("topright", legend = c("AVERAGE"), col = c( "white"), lty = 1, lwd = 1.5, bty = "n")
  
  
  #######################################################
  ######################### MAX ########################
  #######################################################
  
  # Surface layer
  plot(1:730/2, max_var_roms_1, type = "l", col = "#FF4500", lwd = 1.5, xlab = "",
       ylab = "", main = "Layer 1 (0-5m)", ylim = c(abs_min, abs_max), axes = FALSE)
  axis(1, at = c(-100, 0,100,200,300,400))
  axis(2, at = seq(round(abs_min)/2,round(abs_max)*2, (round(abs_max)-round(abs_min))/8))
  # box(col = )
  lines(1:730/2, max_var_1, col = "#1E3A8A", lwd = 1.5)
  lines(1:730/2, max_var_velma_1, col = "#ADFF2F", lwd = 1.5)
  grid(nx = NULL, ny = NULL, col = "lightgray", lty = "dotted")
  legend("topleft", legend = c("ROMS", "SSM", "SSM-VELMA"), col = c("#FF4500", "#1E3A8A", "#ADFF2F"), lty = 1, lwd = 1.5, bty = "n")
  legend("topright", legend = c("MAXIMUM"), col = c( "white"), lty = 1, lwd = 1.5, bty = "n")
  
  # Depth layer
  plot(1:730/2, max_var_roms_5, type = "l", col = "#FF4500", lwd = 1.5, xlab = "",
       ylab = "", main = "Layer 5 (100-150 m)", ylim = c(abs_min, abs_max), axes = FALSE)
  axis(1, at = c(-100,0,100,200,300,400))
  axis(2, at = seq(round(abs_min)/2,round(abs_max)*2, (round(abs_max)-round(abs_min))/8), labels = FALSE)
  
  # box(col = )
  lines(1:730/2, max_var_5, col = "#1E3A8A", lwd = 1.5)
  lines(1:730/2, max_var_velma_5, col = "#ADFF2F", lwd = 1.5)
  grid(nx = NULL, ny = NULL, col = "lightgray", lty = "dotted")
  legend("topleft", legend = c("ROMS", "SSM", "SSM-VELMA"), col = c("#FF4500", "#1E3A8A", "#ADFF2F"), lty = 1, lwd = 1.5, bty = "n")
  legend("topright", legend = c("MAXIMUM"), col = c( "white"), lty = 1, lwd = 1.5, bty = "n")

  
  #######################################################
  ######################### MIN ########################
  #######################################################
  
  # Surface layer
  plot(1:730/2, min_var_roms_1, type = "l", col = "#FF4500", lwd = 1.5, xlab = "",
       ylab = "", main = "Layer 1 (0-5m)", ylim = c(abs_min, abs_max), axes = FALSE)
  axis(1, at = c(-100, 0,100,200,300,400))
  axis(2, at = seq(round(abs_min)/2,round(abs_max)*2, (round(abs_max)-round(abs_min))/8))
  # box(col = )
  lines(1:730/2, min_var_1, col = "#1E3A8A", lwd = 1.5)
  lines(1:730/2, min_var_velma_1, col = "#ADFF2F", lwd = 1.5)
  grid(nx = NULL, ny = NULL, col = "lightgray", lty = "dotted")
  legend("topleft", legend = c("ROMS", "SSM", "SSM-VELMA"), col = c("#FF4500", "#1E3A8A", "#ADFF2F"), lty = 1, lwd = 1.5, bty = "n")
  legend("topright", legend = c("MINIMUM"), col = c( "white"), lty = 1, lwd = 1.5, bty = "n")
  
  # Depth layer
  plot(1:730/2, min_var_roms_5, type = "l", col = "#FF4500", lwd = 1.5, xlab = "",
       ylab = "", main = "Layer 5 (100-150 m)", ylim = c(abs_min, abs_max), axes = FALSE)
  axis(1, at = c(-100,0,100,200,300,400))
  axis(2, at = seq(round(abs_min)/2,round(abs_max)*2, (round(abs_max)-round(abs_min))/8), labels = FALSE)
  
  # box(col = )
  lines(1:730/2, min_var_5, col = "#1E3A8A", lwd = 1.5)
  lines(1:730/2, min_var_velma_5, col = "#ADFF2F", lwd = 1.5)
  grid(nx = NULL, ny = NULL, col = "lightgray", lty = "dotted")
  legend("topleft", legend = c("ROMS", "SSM", "SSM-VELMA"), col = c("#FF4500", "#1E3A8A", "#ADFF2F"), lty = 1, lwd = 1.5, bty = "n")
  legend("topright", legend = c("MINIMUM"), col = c( "white"), lty = 1, lwd = 1.5, bty = "n")
  
  
}
plot_SSM_wq <- function(varname, unit = "", varid = varname, title = title){
  
  ylab = paste0(varname, " ", unit)
  file_SSM <- paste0("Workflow/Step B/Final outputs/No_VELMA/2011/pugetsound_SSM_Atlantis_",varname,"_novelma_2011.nc")
  file_SSM_velma <- paste0("Workflow/Step B/Final outputs/VELMA/2011/pugetsound_SSM_Atlantis_",varname,"_velma_2011.nc")
  
  
  # volume data for weighted average
  atlantis_bgm <- read_bgm("Workflow/Step B/PugetSound_89b_070116.bgm")
  atlantis_sf <- atlantis_bgm %>% box_sf()
  area <- atlantis_sf$area
  
  
  # Obtenir les données de la variable
  nc <- nc_open(file_SSM)
  var <- ncvar_get(nc, varid = varid)
  
  nc <- nc_open(file_SSM_velma)
  var_velma <- ncvar_get(nc, varid = varid)
  
  
  # Transforamtion des données #
  tickness <-  sort((sort(unique(atlantis_sf$botz))[-7] - sort(unique(atlantis_sf$botz))[-1])*-1)
  
  st <- apply(X = var[,,1], FUN = is.na, MARGIN = c(1,2))*-1+1
  layer <-  st
  for (i in 1:89){
    a <- layer[,i]
    a[1:sum(a)] <- sum(a):1
    layer[,i] <- a
  }
  volume <- st
  for (i in 1:6){
    for(j in 1:89){
      if (layer[i,j]>0){
        volume[i,j] <- volume[i,j]*area[j]*tickness[layer[i,j]]
      }
    }
  }
  
  mean_var_overall <- mean_var_1 <- mean_var_2 <- mean_var_3 <- mean_var_4 <- mean_var_5 <- mean_var_6 <- c()
  min_var_overall <- min_var_1 <- min_var_2 <- min_var_3 <- min_var_4 <- min_var_5 <- min_var_6 <- c()
  max_var_overall <- max_var_1 <- max_var_2 <- max_var_3 <- max_var_4 <- max_var_5 <- max_var_6 <- c()
  
  mean_var_velma_overall <- mean_var_velma_1 <- mean_var_velma_2 <- mean_var_velma_3 <- mean_var_velma_4 <- mean_var_velma_5 <- mean_var_velma_6 <- c()
  min_var_velma_overall <- min_var_velma_1 <- min_var_velma_2 <- min_var_velma_3 <- min_var_velma_4 <- min_var_velma_5 <- min_var_velma_6 <- c()
  max_var_velma_overall <- max_var_velma_1 <- max_var_velma_2 <- max_var_velma_3 <- max_var_velma_4 <- max_var_velma_5 <- max_var_velma_6 <- c()
  
  for (i in 1:730){
    mean_var_overall <- c(mean_var_overall, sum(volume*var[,,i]/sum(volume)))
    mean_var_1 <- c(mean_var_1, sum(volume[layer == 1]*var[,,i][layer == 1]/sum(volume[layer == 1])))
    mean_var_2 <- c(mean_var_2, sum(volume[layer == 2]*var[,,i][layer == 2]/sum(volume[layer == 2])))
    mean_var_3 <- c(mean_var_3, sum(volume[layer == 3]*var[,,i][layer == 3]/sum(volume[layer == 3])))
    mean_var_4 <- c(mean_var_4, sum(volume[layer == 4]*var[,,i][layer == 4]/sum(volume[layer == 4])))
    mean_var_5 <- c(mean_var_5, sum(volume[layer == 5]*var[,,i][layer == 5]/sum(volume[layer == 5])))
    mean_var_6 <- c(mean_var_6, sum(volume[layer == 6]*var[,,i][layer == 6]/sum(volume[layer == 6])))
    
    max_var_overall <- c(max_var_overall, max(var[,,i]))
    max_var_1 <- c(max_var_1, max(var[,,i][layer == 1]))
    max_var_2 <- c(max_var_2,max(var[,,i][layer == 2]))
    max_var_3 <- c(max_var_3,max(var[,,i][layer == 3]))
    max_var_4 <- c(max_var_4,max(var[,,i][layer == 4]))
    max_var_5 <- c(max_var_5,max(var[,,i][layer == 5]))
    max_var_6 <- c(max_var_6,max(var[,,i][layer == 6]))
    
    min_var_overall <- c(min_var_overall, min(var[,,i]))
    min_var_1 <- c(min_var_1, min(var[,,i][layer == 1]))
    min_var_2 <- c(min_var_2,min(var[,,i][layer == 2]))
    min_var_3 <- c(min_var_3,min(var[,,i][layer == 3]))
    min_var_4 <- c(min_var_4,min(var[,,i][layer == 4]))
    min_var_5 <- c(min_var_5,min(var[,,i][layer == 5]))
    min_var_6 <- c(min_var_6,min(var[,,i][layer == 6]))
    
    mean_var_velma_overall <- c(mean_var_velma_overall, sum(volume*var_velma[,,i]/sum(volume)))
    mean_var_velma_1 <- c(mean_var_velma_1, sum(volume[layer == 1]*var_velma[,,i][layer == 1]/sum(volume[layer == 1])))
    mean_var_velma_2 <- c(mean_var_velma_2, sum(volume[layer == 2]*var_velma[,,i][layer == 2]/sum(volume[layer == 2])))
    mean_var_velma_3 <- c(mean_var_velma_3, sum(volume[layer == 3]*var_velma[,,i][layer == 3]/sum(volume[layer == 3])))
    mean_var_velma_4 <- c(mean_var_velma_4, sum(volume[layer == 4]*var_velma[,,i][layer == 4]/sum(volume[layer == 4])))
    mean_var_velma_5 <- c(mean_var_velma_5, sum(volume[layer == 5]*var_velma[,,i][layer == 5]/sum(volume[layer == 5])))
    mean_var_velma_6 <- c(mean_var_velma_6, sum(volume[layer == 6]*var_velma[,,i][layer == 6]/sum(volume[layer == 6])))
    
    max_var_velma_overall <- c(max_var_velma_overall, max(var_velma[,,i]))
    max_var_velma_1 <- c(max_var_velma_1, max(var_velma[,,i][layer == 1]))
    max_var_velma_2 <- c(max_var_velma_2,max(var_velma[,,i][layer == 2]))
    max_var_velma_3 <- c(max_var_velma_3,max(var_velma[,,i][layer == 3]))
    max_var_velma_4 <- c(max_var_velma_4,max(var_velma[,,i][layer == 4]))
    max_var_velma_5 <- c(max_var_velma_5,max(var_velma[,,i][layer == 5]))
    max_var_velma_6 <- c(max_var_velma_6,max(var_velma[,,i][layer == 6]))
    
    min_var_velma_overall <- c(min_var_velma_overall, min(var_velma[,,i]))
    min_var_velma_1 <- c(min_var_velma_1, min(var_velma[,,i][layer == 1]))
    min_var_velma_2 <- c(min_var_velma_2,min(var_velma[,,i][layer == 2]))
    min_var_velma_3 <- c(min_var_velma_3,min(var_velma[,,i][layer == 3]))
    min_var_velma_4 <- c(min_var_velma_4,min(var_velma[,,i][layer == 4]))
    min_var_velma_5 <- c(min_var_velma_5,min(var_velma[,,i][layer == 5]))
    min_var_velma_6 <- c(min_var_velma_6,min(var_velma[,,i][layer == 6]))
    
  }
  
  #######################################################
  ######################### MEAN ########################
  #######################################################
  
  
  # pdf(width=9,height = 4.5)
  # Configuration de la disposition des graphiques et ajustement des marges globales
  par(mfrow = c(3, 2), oma = c(2.5, 2.5, 2, 1), mar = c(2, 2, 2, 1))
  
  # Définir les limites de l'axe y
  abs_min <- min(c(min_var_1, min_var_velma_1, min_var_5, min_var_velma_5))
  abs_max <- max(c(max_var_1, max_var_velma_1, max_var_5, max_var_velma_5))*1.2
  
  # Surface layer
  plot(1:730/2, mean_var_1, type = "l", col = "#1E3A8A", lwd = 1.5, xlab = "",
       ylab = "", main = "Layer 1 (0-5m)", ylim = c(abs_min, abs_max), axes = FALSE)
  axis(1, at = c(-100,0,100,200,300,400))
  axis(2, at = seq(round(abs_min)/2,round(abs_max)*2, (round(abs_max)-round(abs_min))/8))
  # box(col = )
  lines(1:730/2, mean_var_velma_1, col = "#ADFF2F", lwd = 1.5)
  grid(nx = NULL, ny = NULL, col = "lightgray", lty = "dotted")
  legend("topleft", legend = c("SSM", "SSM-VELMA"), col = c("#1E3A8A", "#ADFF2F"), lty = 1, lwd = 1.5, bty = "n")
  legend("topright", legend = c("AVERAGE"), col = c( "white"), lty = 1, lwd = 1.5, bty = "n")
  
  # Depth layer
  plot(1:730/2, mean_var_5, type = "l", col = "#1E3A8A", lwd = 1.5, xlab = "",
       ylab = "", main = "Layer 5 (100-150 m)", ylim = c(abs_min, abs_max), axes = FALSE)
  axis(1, at = c(-100,0,100,200,300,400))
  axis(2, at = seq(round(abs_min)/2,round(abs_max)*2, (round(abs_max)-round(abs_min))/8), labels = FALSE)
  
  # box(col = )
  lines(1:730/2, mean_var_velma_5, col = "#ADFF2F", lwd = 1.5)
  grid(nx = NULL, ny = NULL, col = "lightgray", lty = "dotted")
  legend("topleft", legend = c("SSM", "SSM-VELMA"), col = c( "#1E3A8A", "#ADFF2F"), lty = 1, lwd = 1.5, bty = "n")
  legend("topright", legend = c("AVERAGE"), col = c( "white"), lty = 1, lwd = 1.5, bty = "n")
  
  # Ajouter les étiquettes des axes globaux
  mtext("Days", side = 1, outer = TRUE, line = 1, col = , cex = 1.2)
  mtext(ylab, side = 2, outer = TRUE, line = 1, col = , cex = 1.2)
  mtext(title, side = 3, outer = TRUE, line = 0, col = , cex = 1.2)
  
  
  
  #######################################################
  ######################### MAX ########################
  #######################################################
  
  # Configuration de la disposition des graphiques et ajustement des marges globales
  # par(mfrow = c(1, 2), oma = c(2.5, 2.5, 2, 1), mar = c(2, 2, 2, 1))
  
  # Surface layer
  plot(1:730/2, max_var_1, type = "l", col = "#1E3A8A", lwd = 1.5, xlab = "",
       ylab = "", main = "Layer 1 (0-5m)", ylim = c(abs_min, abs_max), axes = FALSE)
  axis(1, at = c(-100,0,100,200,300,400))
  axis(2, at = seq(round(abs_min)/2,round(abs_max)*2, (round(abs_max)-round(abs_min))/8))
  # box(col = )
  lines(1:730/2, max_var_velma_1, col = "#ADFF2F", lwd = 1.5)
  grid(nx = NULL, ny = NULL, col = "lightgray", lty = "dotted")
  legend("topleft", legend = c("SSM", "SSM-VELMA"), col = c("#1E3A8A", "#ADFF2F"), lty = 1, lwd = 1.5, bty = "n")
  legend("topright", legend = c("MAXIMUM"), col = c( "white"), lty = 1, lwd = 1.5, bty = "n")
  
  # Depth layer
  plot(1:730/2, max_var_5, type = "l", col = "#1E3A8A", lwd = 1.5, xlab = "",
       ylab = "", main = "Layer 5 (100-150 m)", ylim = c(abs_min, abs_max), axes = FALSE)
  axis(1, at = c(-100,0,100,200,300,400))
  axis(2, at = seq(round(abs_min)/2,round(abs_max)*2, (round(abs_max)-round(abs_min))/8), labels = FALSE)
  
  # box(col = )
  lines(1:730/2, max_var_velma_5, col = "#ADFF2F", lwd = 1.5)
  grid(nx = NULL, ny = NULL, col = "lightgray", lty = "dotted")
  legend("topleft", legend = c("SSM", "SSM-VELMA"), col = c( "#1E3A8A", "#ADFF2F"), lty = 1, lwd = 1.5, bty = "n")
  legend("topright", legend = c("MAXIMUM"), col = c( "white"), lty = 1, lwd = 1.5, bty = "n")
  
  #######################################################
  ######################### MIN ########################
  #######################################################
  
  # Configuration de la disposition des graphiques et ajustement des marges globales
  # par(mfrow = c(1, 2), oma = c(2.5, 2.5, 2, 1), mar = c(2, 2, 2, 1))
  
  # Surface layer
  plot(1:730/2, min_var_1, type = "l", col = "#1E3A8A", lwd = 1.5, xlab = "",
       ylab = "", main = "Layer 1 (0-5m)", ylim = c(abs_min, abs_max), axes = FALSE)
  axis(1, at = c(-100,0,100,200,300,400))
  axis(2, at = seq(round(abs_min)/2,round(abs_max)*2, (round(abs_max)-round(abs_min))/8))
  # box(col = )
  lines(1:730/2, min_var_velma_1, col = "#ADFF2F", lwd = 1.5)
  grid(nx = NULL, ny = NULL, col = "lightgray", lty = "dotted")
  legend("topleft", legend = c("SSM", "SSM-VELMA"), col = c("#1E3A8A", "#ADFF2F"), lty = 1, lwd = 1.5, bty = "n")
  legend("topright", legend = c("MINIMUM"), col = c( "white"), lty = 1, lwd = 1.5, bty = "n")
  
  # Depth layer
  plot(1:730/2, min_var_5, type = "l", col = "#1E3A8A", lwd = 1.5, xlab = "",
       ylab = "", main = "Layer 5 (100-150 m)", ylim = c(abs_min, abs_max), axes = FALSE)
  axis(1, at = c(-100,0,100,200,300,400))
  axis(2, at = seq(round(abs_min)/2,round(abs_max)*2, (round(abs_max)-round(abs_min))/8), labels = FALSE)
  
  # box(col = )
  lines(1:730/2, min_var_velma_5, col = "#ADFF2F", lwd = 1.5)
  grid(nx = NULL, ny = NULL, col = "lightgray", lty = "dotted")
  legend("topleft", legend = c("SSM", "SSM-VELMA"), col = c( "#1E3A8A", "#ADFF2F"), lty = 1, lwd = 1.5, bty = "n")
  legend("topright", legend = c("MINIMUM"), col = c( "white"), lty = 1, lwd = 1.5, bty = "n")
  
}
plot_SSM <- function(varname){
  
  list.var = c("salinity","temperature",
               "NO3", "NH4", 
               "SZ", "LZ", "MZ", "SP", "LP",
               "Oxygen","LPON","RPON",
               "RDON")
  
  list.title = c("Salinity","Temperature",
               "Nitrate", "Ammonium", 
               "Small Zooplankton", "Large Zooplankton", "Meso Zooplankton", 
               "Small Phytoplankton", "Large Phytoplankton",
               "Oxygen","Labile Particulate Organic Nitrogen","Refractory Particulate Organic Nitrogen",
               "Dissolved Organic Nitrogen")
  
  list.unit <- c("(g.L-1)","(°C)",
                 "(mg N.m-3)", "(mg N.m-3)", 
                 "(mg N.m-3)", "(mg N.m-3)", "(mg N.m-3)", "(mg N.m-3)", "(mg N.m-3)",
                 "(mg O2.m-3)","(mg N.m-3)","(mg N.m-3)",
                 "(mg N.m-3)")
  
  atlantis.name <- c("salinity","temperature",
                     "NO3", "NH3", 
                     "Micro_Zoo_N", "Lrg_Zoo_N", "Meso_Zoo_N", "Sm_Phyto_N", "Lrg_Phyto_N",
                     "Oxygen","Lab_Det_N","Ref_Det_N",
                     "DON")
  
  if (varname == "all"){
    pdf(file = "SSM-VELMA.pdf",   # The directory you want to save the file in
        width = 10, # The width of the plot in inches
        height = 9)
    plot_SSM(varname = "salinity")
    plot_SSM(varname = "temperature")
    plot_SSM(varname = "Oxygen")

    plot_SSM(varname = "NH4")
    plot_SSM(varname = "NO3")
    plot_SSM(varname = "SP")
    plot_SSM(varname = "LP")
    plot_SSM(varname = "SZ")
    plot_SSM(varname = "MZ")
    plot_SSM(varname = "LZ")
    plot_SSM(varname = "LPON")
    plot_SSM(varname = "RPON")
    plot_SSM(varname = "RDON")
    dev.off()
  }else{
    
    
    if(!varname%in%list.var) stop("The variable is not in SSM, please try:
salinity, temperature, NO3, NH4, SZ, LZ, MZ, SP, LP, Oxygen, LPON, RPON, RDON. 
For a pdf sum up, try all ")
    
    pos = match(varname, list.var)  
    
    if (varname=="salinity"||varname=="temperature"){
      plot_SSM_physics(varname = varname, unit = list.unit[pos], 
                       varid = atlantis.name[pos], title = list.title[pos])
    }else{
      plot_SSM_wq(varname = varname, unit = list.unit[pos], 
                  varid = atlantis.name[pos], title = list.title[pos])
    }
    
    
    
  }
  
  
  
}
plot_SSM("all")
plot_SSM(varname = "salinity")
plot_SSM(varname = "temperature")
plot_SSM(varname = "Oxygen")
plot_SSM(varname = "NH4")
plot_SSM(varname = "NO3")
plot_SSM(varname = "SP")
plot_SSM(varname = "LP")
plot_SSM(varname = "SZ")
plot_SSM(varname = "MZ")
plot_SSM(varname = "LZ")
plot_SSM(varname = "LPON")
plot_SSM(varname = "RPON")
plot_SSM(varname = "RDON")




# plot_SSM(varname = "LPON")
# plot_SSM(varname = "RPON")

