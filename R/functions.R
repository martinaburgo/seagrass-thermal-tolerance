# CLIMATE VARIABLES CALCULATIONS

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

# because of the land sea mask in the temperature raster, some coastal cells where seagrass are present
# do not have temperature values
# fix this by assigning missing pixels as the mean of the surrounding 8 cells (i.e. queens move)

fill.na <- function(x, i=5) {
  if( is.na(x)[i] ) {
    return(mean(x, na.rm=TRUE))
  } else {
    return(x[i])
  }
} 

# https://gis.stackexchange.com/questions/181011/fill-the-gaps-using-nearest-neighbors

# check the function works i.e. fills in NA values but doesn't alter numeric values
rmat <- matrix(rnorm(64,0,1), nrow=8, ncol=8)
r <- raster(rmat)
r[1:2,1:2] <- NA
r2 <- focal(r, w = matrix(1,3,3), fun = fill.na,pad = TRUE, na.rm = FALSE)
plot(values(r),values(r2))
values(r)-values(r2)
mean(r[3,1:2])
values(r2)[which(is.na(values(r))==TRUE)]
