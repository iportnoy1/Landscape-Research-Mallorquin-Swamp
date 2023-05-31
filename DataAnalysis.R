##Setting working directory
setwd("C:/Users/idpdl/Desktop/Landscape Research")

##Calling package
library("raster"); library("SpatialPack"); library("tiff")

##Loading images generated in Google Earth Engine (GEE)
NDVI_2013_I <- raster("NDVI_2013-I.tif"); NDVI_2013_II <- raster("NDVI_2013-II.tif")
NDVI_2018_I <- raster("NDVI_2018-I.tif"); NDVI_2018_II <- raster("NDVI_2018-II.tif")
NDVI_2022_I <- raster("NDVI_2022-I.tif"); NDVI_2022_II <- raster("NDVI_2022-II.tif")

##Plotting
par(mfrow=c(3,2))
plot(NDVI_2013_I, xlab="Latitude", ylab="Longitude", main="NDVI,  2013-I")
plot(NDVI_2013_II, xlab="Latitude", ylab="Longitude", main="NDVI,  2013-II")
plot(NDVI_2018_I, xlab="Latitude", ylab="Longitude", main="NDVI,  2018-I")
plot(NDVI_2018_II, xlab="Latitude", ylab="Longitude", main="NDVI,  2018-II")
plot(NDVI_2022_I, xlab="Latitude", ylab="Longitude", main="NDVI,  2022-I")
plot(NDVI_2022_II, xlab="Latitude", ylab="Longitude", main="NDVI,  2022-II")

##Loading raw NDVI values
X_2013_I <- readTIFF("NDVI_2013-I.tif", native = FALSE)
X_2013_II <- readTIFF("NDVI_2013-II.tif", native = FALSE)
X_2018_I <- readTIFF("NDVI_2018-I.tif", native = FALSE)
X_2018_II <- readTIFF("NDVI_2018-II.tif", native = FALSE)
X_2022_I <- readTIFF("NDVI_2022-I.tif", native = FALSE)
X_2022_II <- readTIFF("NDVI_2022-II.tif", native = FALSE)

##Getting rid of the common non-vegetable-populated areas
cond<-(X_2013_I<0.2)&(X_2013_II<0.2)&(X_2018_I<0.2)&(X_2018_II<0.2)&(X_2022_I<0.2)&(X_2022_II<0.2)
X_2013_I[cond] <- NA
X_2013_I_na.rm <- X_2013_I; X_2013_I_na.rm[is.na(X_2013_I)] <- 0
X_2013_II[cond] <- NA
X_2013_II_na.rm <- X_2013_II; X_2013_II_na.rm[is.na(X_2013_II)] <- 0
X_2018_I[cond] <- NA
X_2018_I_na.rm <- X_2018_I; X_2018_I_na.rm[is.na(X_2018_I)] <- 0
X_2018_II[cond] <- NA
X_2018_II_na.rm <- X_2018_II; X_2018_II_na.rm[is.na(X_2018_II)] <- 0
X_2022_I[cond] <- NA
X_2022_I_na.rm <- X_2022_I; X_2022_I_na.rm[is.na(X_2022_I)] <- 0
X_2022_II[cond] <- NA
X_2022_II_na.rm <- X_2022_II; X_2022_II_na.rm[is.na(X_2022_II)] <- 0

##For Comparison with the NDVI mean values
ndvi.mean.vals <- c(mean(X_2013_I, na.rm = T), mean(X_2013_II, na.rm = T), mean(X_2018_I,na.rm = T), 
                    mean(X_2018_II, na.rm = T), mean(X_2022_I, na.rm = T), mean(X_2022_II, na.rm = T))
plot(ndvi.mean.vals, xlab = "Period", ylab = "Average NDVI", xaxt = "n", col = "blue", lty = 2, 
     type = 'b', ylim = c(0.28,0.4))
axis(1, at=c(1:6), labels=c("2013-I", "2013-II", "2018-I", "2018-II", "2022-I", "2022-II"))

##Assessing significance of the observed difference with bootstrapping
#Creating Bootstrapping function
statistic <- function(x1, x2) {
  mean(x1, na.rm=T) - mean(x2, na.rm=T)
}
statistic <- compiler:::cmpfun(statistic)
bS <- function(x1, x2, Bootstraps = 10000){
  n1 <- nrow(x1); n2 <- nrow(x2); m1 <- ncol(x1); m2 <- ncol(x2)
  X <- rbind(x1, x2)
  out <- lapply(1:Bootstraps, function(b){
    #Scrambling rows
    x1.samp <- X[sample(nrow(X), n1, replace = T), ]
    x2.samp <- X[sample(nrow(X), n2, replace = T), ]
    #Scrambling columns
    x.samp <- rbind(x1.samp, x2.samp)
    x1.samp2 <-  x.samp[,sample(ncol(x.samp), m1, replace = T)]
    x2.samp2 <-  x.samp[,sample(ncol(x.samp), m2, replace = T)]
    out <-statistic(x1.samp2, x2.samp2)
  })
  unlist(out)		
}
bS <- compiler:::cmpfun(bS)
comparison_test <- function(x1, x2, Bootstraps = 10000, plot.histogram = TRUE){
  stat <- statistic(x1, x2)
  res <- bS(x1, x2, Bootstraps)
  # if(stat < mean(res)){
  #   p <- mean(res<stat)
  # }else{
  #   p <- mean(res>stat)
  # }
  p <- mean(abs(res) >= abs(stat))
  if(plot.histogram == TRUE){
    hist(res, yaxt = "n", prob = TRUE, ylab = "", col = 1, border = 1, xlab = 'statistic', 
         main = '', xlim=c(min(c(res,stat))-0.001,max(c(res,stat))+0.001)); 
    abline(v = stat, col = "red", lwd = 2)
  }
  out = list(p,stat, res)
  names(out) <- c("p.value", "statistic", "bootstrapped_values")
  out
}
comparison_test <- compiler:::cmpfun(comparison_test)

Comp.Matrix <- matrix(data = rep(0,15), nrow = 3, ncol = 5)
colnames(Comp.Matrix) <- c("2013-II", "2018-I", "2018-II", "2022-I", "2022-II")
rownames(Comp.Matrix) <- c("SSIM", "Mean Difference Value", "Meand Difference p-value")

Comp.Matrix[1,1] <- SSIM(X_2013_II_na.rm, X_2013_I_na.rm, alpha=1, beta=1, gamma=1, 
                         eps=c(0.01, 0.03), L=2)$SSIM
Comp.Matrix[1,2] <- SSIM(X_2018_I_na.rm, X_2013_I_na.rm, alpha=1, beta=1, gamma=1, 
                         eps=c(0.01, 0.03), L=2)$SSIM
Comp.Matrix[1,3] <- SSIM(X_2018_II_na.rm, X_2013_I_na.rm, alpha=1, beta=1, gamma=1, 
                         eps=c(0.01, 0.03), L=2)$SSIM
Comp.Matrix[1,4] <- SSIM(X_2022_I_na.rm, X_2013_I_na.rm, alpha=1, beta=1, gamma=1, 
                         eps=c(0.01, 0.03), L=2)$SSIM
Comp.Matrix[1,5] <- SSIM(X_2022_II_na.rm, X_2013_I_na.rm, alpha=1, beta=1, gamma=1, 
                         eps=c(0.01, 0.03), L=2)$SSIM
temp <- comparison_test(X_2013_II, X_2013_I, Bootstraps=10000, plot.histogram=F)
Comp.Matrix[2,1] <- temp$statistic
Comp.Matrix[3,1] <- temp$p.value
temp <- comparison_test(X_2018_I, X_2013_I, Bootstraps=10000, plot.histogram=F)
Comp.Matrix[2,2] <- temp$statistic
Comp.Matrix[3,2] <- temp$p.value
temp <- comparison_test(X_2018_II, X_2013_I, Bootstraps=10000, plot.histogram=F)
Comp.Matrix[2,3] <- temp$statistic
Comp.Matrix[3,3] <- temp$p.value
temp <- comparison_test(X_2022_I, X_2013_I, Bootstraps=10000, plot.histogram=F)
Comp.Matrix[2,4] <- temp$statistic
Comp.Matrix[3,4] <- temp$p.value
temp <- comparison_test(X_2022_II, X_2013_I, Bootstraps=10000, plot.histogram=F)
Comp.Matrix[2,5] <- temp$statistic
Comp.Matrix[3,5] <- temp$p.value
