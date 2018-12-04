
# load hf travel time rasters
hf.rast.files<-c(paste("J:/Madagascar_time_travel_rasters/", c(1:3247), "HF.access.tif", sep =""))
#stack
hf.rast<-stack(hf.rast.files)
