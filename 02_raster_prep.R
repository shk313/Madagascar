#################################################################################################
########################       Script for 02_raster_prep                 #########################
##################################################################################################

### used in Arambepola et al Spatiotemporal mapping of malaria prevalence in Madgascar using routine surveillance and health survey data


library(raster)
library(rgdal)
library(malariaAtlas)
library(doParallel)

## moving to Z to free up disc space on the blade
#hf.rast.files <-cbind(c(list.files("/home/suzanne/MDG_time_travel_rasters", pattern ="tif$", full.names = TRUE)))
hf.rast.files <-cbind(c(list.files("/home/suzanne/zmount/Madagascar-NMCP/MDG_time_travel_rasters", pattern ="tif$", full.names = TRUE)))

new_folder <- ("/home/suzanne/Madagascar_time_travel_rasters_cropped/")
dir.create(new_folder)

MDG <- getShp(ISO = "MDG")  ## shapefile to crop to

######################
##### Step 1: crop and re-write files to remove water from the rasters.
######################

cores <- detectCores() 
cl <- makeCluster(10)   ##if running on personal computer this number will change
registerDoParallel(cl)
start_time <- Sys.time()
b_crop <- foreach(i = 1:nrow(hf.rast.files)) %dopar% {
  
  b <- raster::raster(hf.rast.files[i])
  b_crop <- raster::crop(b, MDG) 
  
  raster::writeRaster(b_crop, filename = file.path(paste0(new_folder, names(b), ".tif"))) 
  
}

end_time <- Sys.time()
end_time - start_time  ## 10mins
stopCluster(cl)

##alternative to above
test <- lapply(hf.rast.files, function(x) raster(x))
test_crop <- lapply(test, function(x) crop(x, MDG))
lapply(test_crop, function(x) writeRaster(x, filename=file.path(paste0(new_folder, names(x), ".tif"))))

# setwd("/home/suzanne/Madagascar_time_travel_rasters_cropped/")
# mapply(writeRaster, test_crop, names(test_crop), 'GTiff')


## make sure that all of our masked rasters have the same number of cells
rasters_check <- sapply(seq_along(cropped), function(x) ncell(cropped[[x]]))
if(length(unique(rasters_check)) > 1){
  getmode <- function(v) {
    uniqv <- unique(v)
    uniqv[which.max(tabulate(match(v, uniqv)))]
  }
  
  warning("Check the following files: ", 
          paste(filenames[which(is.na(match(rasters_check, getmode(rasters_check))))], 
                collapse = "\n"), 
          sep = "")
  
}
rm(rasters_check)

#####################################
## Step 2: 1/travel time squared.
#####################################

#hf.rast.files<-cbind(c(paste("/home/suzanne/Madagascar_time_travel_rasters/X", c(1:3247), "HF.access.tif", sep ="")))
hf.rast.files <-cbind(c(list.files("/home/suzanne/Madagascar_time_travel_rasters_cropped", pattern ="tif$", full.names = TRUE)))


nrow <- ncell(raster(hf.rast.files[1]))
ncol <- length(hf.rast.files)
mat <- matrix(nrow = nrow, ncol = ncol)  ## pre defining should speed things up a little.

### save travel time distances as a matrix prior to inverse
start_time <- Sys.time()
for(i in 1:nrow(hf.rast.files)){
  a <- raster(hf.rast.files[[i]])
  x <- as.matrix(a)
  
  ind <- i*1
  mat[, ind] <- unlist(x, use.names = TRUE)
  
}
end_time <- Sys.time()
end_time - start_time 

save(mat, file = "/home/suzanne/hf_matrix.RData")

### inverse travel time distances
start_time <- Sys.time()
for(i in 1:nrow(hf.rast.files)){
  
  a <- raster(hf.rast.files[[i]])
  x <- as.matrix(a)
  
  x <- 1/(x^2)  
  x <- list(x)
  
  ind <- i*1
  mat_inv[, ind] <- unlist(x, use.names = TRUE)
  
  
  print(paste("written", names(a)))
  
}  

end_time <- Sys.time()
end_time - start_time ## 16 mins


### parralelise, check if any quicker.
cores <- detectCores() 
cl <- makeCluster(10)   ##if running on personal computer this number will change
registerDoParallel(cl)
start_time <- Sys.time()
mat_inv <- foreach(i = 1:nrow(hf.rast.files)) %dopar% {
  
  a <- raster(hf.rast.files[[i]])
  x <- as.matrix(a)
  
  x <- 1/(x^2)  
  x <- list(x)
  
  ind <- i*1
  mat[, ind] <- unlist(x, use.names = TRUE)
  
  
  print(paste("written", names(a)))
  
}  

end_time <- Sys.time()
end_time - start_time ## 
stopCluster(cl)

save(mat, file = "/home/suzanne/hf_inv_matrix.RData")
## latest save 25/06/2019

###############
## Step 3: Set cut offs
################

mat <- get(load("/home/suzanne/hf_matrix.RData"))
mat_inv <- get(load("/home/suzanne/hf_inv_matrix.RData"))

mat_inv[mat < 1] <- 1
mat_inv_old <- mat_inv
mat_inv[mat > 200] <- 0

## if any pixel has no hf within 200 mins, assign it to the nearest hf.
which_zero <- which(rowSums(mat_inv) == 0)
for(pixel_i in which_zero){
  min_hf <- which.min(mat_inv_old[pixel_i, ])
  mat_inv[pixel_i, min_hf] <- mat_inv_old[pixel_i, min_hf]
}

catchments <- mat_inv^2


#######
# Step 4: Derive the proportion of the pixel population going to each health facility
#######

rm(proportion)
proportion <- matrix(nrow = 1434125, ncol = 2801)

## separate column names
column_names <- colnames(mat_inv)
colnames(mat) <- NULL
rownames(mat) <- NULL ## slowed down the next bit having these

start_time <- Sys.time()
proportion <- t(apply(mat_inv, 1, function(x) x/sum(x)))  ## na.rm = TRUE, should break if there are NA's
end_time <- Sys.time()
end_time - start_time   ###9 mins used about 70% of mem in bld1

### do any pixels have no people going to a hf? i.e. NA in all cols.
missing_any <- proportion[rowSums(is.na(proportion)) !=ncol(proportion), ]  # 55%


###############
# Step 5: Prepare population file
###############

pop <- raster("/home/suzanne/Utility_files/FB_HRSL_Pop_AfricaMostly_MGMatched.2015.Annual.Data.1km.sum.tif")
mdg_hf <- raster("/home/suzanne/Madagascar_time_travel_rasters_cropped/X3278HF.access.tif")
pop <- crop(pop, mdg_hf)
#pop_matrix <- matrix(pop, nrow = 1434125, ncol = 1)  
#pop_rast <- values(pop)


#############
## Step 6: Combine TS work
#############

ts <- raster("/home/suzanne/Utility_files/ts_logistic.tif")
ts <- crop(ts, mdg_hf)

pop_ts <- pop*ts
pop_ts_matrix <- matrix(pop_ts, nrow= 1434125, ncol = 1)
pop_ts_matrix <- as.numeric(pop_ts_matrix)

proportion_2 <- matrix(nrow = 1434125, ncol = 2801)
proportion_2 <- proportion*pop_ts_matrix  

#can't just use colsums
population <- as.data.frame(t(apply(proportion_2, 2, function(x) sum(x, na.rm = TRUE))))
pop_vec <- t(population)


a <- read.csv("~/catchment_populations_old.csv")
a <- a[,2]

pop_vec <- as.data.frame(cbind(pop_vec, a))
colnames(pop_vec)[2] <- "FID"

pop_vec_coords <- left_join(pop_vec, lat, by = "FID")

write.csv(pop_vec_coords, "catchment_populations.csv", row.names=FALSE)




population <- read.csv("catchment_populations.csv")
save(proportion_2, file = "catchment_pop_matrix.RData")
