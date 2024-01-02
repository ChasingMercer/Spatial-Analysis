# Load packages using library() function
library("tmap")
library("sf")
library("spdep")
library("sp")
library("spatialreg")
library("sfdep")
library("spgwr")
library("dplyr")
library("car")

uk_msoa <- read_sf('data/infuse_msoa_lyr_2011/infuse_msoa_lyr_2011.shp', crs="OSGB36")
uk_msoa <- uk_msoa[c("geo_code", "geo_label", "geometry")]

df_features <- read.csv('data/merge_grouped_msoa.csv')
df_features <- df_features[c("msoa_geo_code", "mean_percent_over_65", "ru_cat",
              "mean_imd_decile", "mean_opioid_prescribing_rate", "mean_samhi_index")]

var_cols <- c("mean_percent_over_65", "ru_cat", "mean_imd_decile", 
                   "mean_opioid_prescribing_rate", "mean_samhi_index")

df_features[var_cols] <- scale(df_features[var_cols], center = TRUE, scale = TRUE)

# Set up the layout for the combined plot
par(mfrow = c(2, 3))

# Loop through columns
for (col in var_cols) {
  hist(df_features[[col]], main = " ", xlab = col, col = "lightblue", border = "black")
}
dev.off()

spatialdatafile <- merge(uk_msoa, df_features, by.x="geo_code", by.y="msoa_geo_code")
spatialdatafile <- st_as_sf(spatialdatafile, sf_column_name="geometry", crs="OSGB36")
spatialdatafile <- spatialdatafile %>% mutate(ROWNUM = row_number())
spatialdatafile <- st_set_crs(spatialdatafile, 27700)
spatialdatafile$ROWNUM <- 1:nrow(spatialdatafile)

# tm_shape(spatialdatafile) +
#   tm_fill("mean_opioid_prescribing_rate", style = "cont", n = 5, palette = "magma") +
#   tm_shape(spatialdatafile) + tm_polygons(alpha = 0, border.alpha = 1, border.col = "black", lwd=0.1) +
#   tm_compass(position = c("right", "top")) +
#   tm_scalebar(position = c("left", "bottom")) +
#   tm_layout(frame = FALSE, legend.title.size = 0.5, legend.text.size = 0.5)
# 
# tmap_save(filename = "plots/opioid_distributions.png")
# 
# # map for public transport accessibility categories (PTACAT)
# tm_shape(spatialdatafile) + tm_fill("mean_samhi_index", style = "cont", n=5, palette = "magma") +
#   tm_shape(spatialdatafile) + tm_polygons(alpha = 0, border.alpha = 1, border.col = "black",  lwd=0.1) +
#   tm_compass(position = c("right", "top")) +
#   tm_scalebar(position = c("left", "bottom")) +
#   tm_layout(frame = FALSE, legend.title.size = 0.5, legend.text.size = 0.5)
# 
# 
# tmap_save(filename = "plots/samhi_distributions.png")



Weights <- poly2nb(spatialdatafile, row.names = spatialdatafile$ROWNUM)

neighbor_counts <- card(Weights)

# Find indices of observations with zero neighbors
zero_neighbor_indices <- which(neighbor_counts == 0)

spatialdatafile <- subset(spatialdatafile,
                                subset = rownames(spatialdatafile) != zero_neighbor_indices)

Weights <- poly2nb(spatialdatafile, row.names = spatialdatafile$ROWNUM)

WeightsMatrix <- nb2mat(Weights, style='B', zero.policy=FALSE)
Residual_WeightMatrix <- mat2listw(WeightsMatrix , style='W', zero.policy=FALSE)

modelMLR1 <- lm(mean_opioid_prescribing_rate ~ mean_percent_over_65 + 
                 ru_cat + mean_imd_decile, data = spatialdatafile)

modelMLR <- lm(mean_opioid_prescribing_rate ~ mean_percent_over_65 + 
                 ru_cat + mean_imd_decile + mean_samhi_index, data = spatialdatafile)

spatialdatafile$RESIDUALS <- modelMLR$residuals

summary(modelMLR1)
summary(modelMLR)

summary(modelMLR)$coefficients[,4]

coef(modelMLR)

lm.morantest(modelMLR, Residual_WeightMatrix, alternative="two.sided")
vif_MLR <- vif(modelMLR)



modelSLY <- lagsarlm(mean_opioid_prescribing_rate ~ mean_percent_over_65 + 
                       ru_cat + mean_imd_decile + mean_samhi_index, data = spatialdatafile, Residual_WeightMatrix)

Weights_2.0 <- as(Residual_WeightMatrix, "CsparseMatrix")
trMC <- trW(Weights_2.0, type="MC")

IMPACTS_SLY <- impacts(modelSLY, tr = trMC, R=100)

sink("results/impacts_sly.txt")
print(IMPACTS_SLY)
sink()

spatialdatafile$RESID_SLY <- modelSLY$residuals

moran.mc(spatialdatafile$RESID_SLY, Residual_WeightMatrix, 1000, zero.policy = T)

spatialdatafile <- st_set_crs(spatialdatafile, 27700)

tm_shape(spatialdatafile) +
  tm_fill("RESIDUALS", title = "Residuals MLR", style = "cont", midpoint = 0, palette = "magma") +
  tm_shape(spatialdatafile) +
  tm_polygons(alpha = 0, border.alpha = 0.3, border.col = "black", lwd=0.1) +
  tm_compass(position = c("left", "top")) +
  tm_scalebar(position = c("right", "bottom")) +
  tm_layout(frame = FALSE, legend.title.size = 0.5, legend.text.size = 0.5)

tmap_save(filename = "plots/residuals_MLR.png")

tm_shape(spatialdatafile) +
  tm_fill("RESID_SLY", title = "Residuals SLY", style = "cont", midpoint = 0, palette = "magma") +
  tm_shape(spatialdatafile) +
  tm_polygons(alpha = 0, border.alpha = 0.3, border.col = "black", lwd=0.1) +
  tm_compass(position = c("left", "top")) +
  tm_scalebar(position = c("right", "bottom")) +
  tm_layout(frame = FALSE, legend.title.size = 0.5, legend.text.size = 0.5)

tmap_save(filename = "plots/residuals_SLY.png")

vif <- vif(modelMLR)
moran <- lm.morantest(modelMLR, Residual_WeightMatrix, alternative="two.sided")

spatialdatafile <- st_centroid(spatialdatafile)
spatialdatafile <- cbind(spatialdatafile, st_coordinates(spatialdatafile))

na_count <- colSums(is.na(spatialdatafile))
print(na_count)

BwG <- gwr.sel(mean_opioid_prescribing_rate ~ mean_percent_over_65 + 
                 ru_cat + mean_imd_decile + mean_samhi_index, data = spatialdatafile,
               coords = cbind(spatialdatafile$X, spatialdatafile$Y), adapt = TRUE)


# start timer to time how long it takes to run a gwr() on computer
start.timer <- proc.time()

# gwr() model. You need hatmatrix and se.fit specified as TRUE for testing statistical significance 
gwr.model <- gwr(mean_opioid_prescribing_rate ~ mean_percent_over_65 + ru_cat +
                   mean_imd_decile + mean_samhi_index, data = spatialdatafile, coords = cbind(spatialdatafile$X, spatialdatafile$Y), adapt=BwG, hatmatrix=TRUE, se.fit=TRUE)

gwr.data <- as.data.frame(gwr.model$SDF)

write.csv(gwr.data, 'results/gwr_data.csv')

gwr.data <- read.csv('results/gwr_data.csv')

gwr.data$pred - spatialdatafile$mean_opioid_prescribing_rate

r2 <- 0.876054
adjusted_rsquared <- 1 - ((1 - r2) * (6790 - 1) / (6790 - 4 - 1))
print(adjusted_rsquared)
hist(gwr.data[["localR2"]], main = "Local R-Squared", xlab = "Local R-Squared", col = "lightblue", border = "black")

#Create spatial dataframe containing model statistics
msoa_result <-as.data.frame(spatialdatafile[c("geo_code", "geo_label")])
msoa_result <- merge(msoa_result, uk_msoa, by.x="geo_code", by.y="geo_code")
msoa_result <- msoa_result[c("geo_code", "geo_label.x", "geometry.y")]
msoa_result <- rename(msoa_result, geometry = geometry.y, geo_label = geo_label.x)

msoa_result$CoefPOver65 <- gwr.data[,"mean_percent_over_65"]
msoa_result$CoefRuCat <- gwr.data[,"ru_cat"]
msoa_result$CoefMuSamhi <- gwr.data[,"mean_samhi_index"]
msoa_result$CoefMuImd <- gwr.data[,"mean_imd_decile"]
msoa_result$SePOver65 <- gwr.data[,"mean_percent_over_65_se"]
msoa_result$SeRuCat <- gwr.data[,"ru_cat_se"]
msoa_result$SeMuSamhi <- gwr.data[,"mean_samhi_index_se"]
msoa_result$SeMuImd <- gwr.data[,"mean_imd_decile_se"]
msoa_result$tstatPOver65  <- msoa_result$CoefPOver65  / msoa_result$SePOver65
msoa_result$tstatMuRuCat <- msoa_result$CoefRuCat / msoa_result$SeRuCat
msoa_result$tstatMuSamhi <- msoa_result$CoefMuSamhi / msoa_result$SeMuSamhi
msoa_result$tstatMuImd <- msoa_result$CoefMuImd / msoa_result$SeMuImd

msoa_result$tstatPOver65  <- msoa_result$CoefPOver65  / msoa_result$SePOver65
msoa_result$tstatRuCat <- msoa_result$CoefRuCat / msoa_result$SeRuCat
msoa_result$tstatMuSamhi <- msoa_result$CoefMuSamhi / msoa_result$SeMuSamhi
msoa_result$tstatMuImd <- msoa_result$CoefMuImd / msoa_result$SeMuImd
msoa_result$localR2 <- gwr.data[,"localR2"]
msoa_result$residuals <- gwr.model$lm$residuals

write.csv(as.data.frame(msoa_result), "results/msoa_result.csv")
write.table(BwG, 'results/bwg.txt')
write.table(vif, 'results/vif.txt')

print(moran)

msoa_result <- st_as_sf(msoa_result)
msoa_result <- st_set_crs(msoa_result, value=27700)

st_crs(msoa_result)

tm_shape(msoa_result) + 
  tm_fill("CoefPOver65", title = "Coefficient: Population Over 65 [%]", style = "cont", midpoint = 0, palette = "magma") +
  tm_polygons(alpha = 0, border.alpha = 0.5, border.col = "black",  lwd=0.1) +
  tm_compass(position = c("left", "top")) +
  tm_scalebar(position = c("right", "bottom")) +
  tm_layout(frame = FALSE, legend.title.size = 0.5, legend.text.size = 0.5,
            legend.position=c("right", "top"))

tmap_save(filename = "plots/over65.png")

# Map 2: Coefficient: Rural Urban Category
tm_shape(msoa_result) + 
  tm_fill("CoefRuCat", title = "Coefficient: Rural Urban Category [%]", style = "cont", 
          midpoint = 0, palette = "magma") +
  tm_polygons(alpha = 0, border.alpha = 0.5, border.col = "black",  lwd=0.1) +
  tm_compass(position = c("left", "top")) +
  tm_scalebar(position = c("right", "bottom")) +
  tm_layout(frame = FALSE, legend.title.size = 0.5, legend.text.size = 0.5,,
            legend.position=c("right", "top"))

tmap_save(filename = "plots/rural.png")

# Map 3: Coefficient: Small Area Mental Health Index
tm_shape(msoa_result) + 
  tm_fill("CoefMuSamhi", title = "Coefficient: Small Area Mental Health Index", style = "cont", midpoint = 0, palette = "magma") +
  tm_polygons(alpha = 0, border.alpha = 0.5, border.col = "black",  lwd=0.1) +
  tm_compass(position = c("left", "top")) +
  tm_scalebar(position = c("right", "bottom")) +
  tm_layout(frame = FALSE, legend.title.size = 0.5, legend.text.size = 0.5, ,
            legend.position=c("right", "top"))

tmap_save(filename = "plots/samhi.png")

# 
tm_shape(msoa_result) + 
  tm_fill("CoefMuImd", title = "Coefficient: Index Multiple Deprivation", style = "cont", midpoint = 0, palette = "magma") +
  tm_polygons(alpha = 0, border.alpha = 0.5, border.col = "black",  lwd=0.1) +
  tm_compass(position = c("left", "top")) +
  tm_scalebar(position = c("right", "bottom")) +
  tm_layout(frame = FALSE, legend.title.size = 0.5, legend.text.size = 0.5,
            legend.position=c("right", "top"))

tmap_save(filename = "plots/imd.png")

# Map 4: Coefficient: Index Multiple Deprivation
tm_shape(msoa_result) + 
  tm_fill("localR2", title = "Adaptive: Local R2", style = "cont", midpoint = NA, palette = "magma") +
  tm_polygons(alpha = 0, border.alpha = 0.5, border.col = "black",  lwd=0.1) +
  tm_compass(position = c("left", "top")) +
  tm_scalebar(position = c("right", "bottom")) +
  tm_layout(frame = FALSE, legend.title.size = 0.5, legend.text.size = 0.5, 
            legend.position=c("right", "top"))

tmap_save(filename = "plots/rsquared.png")

# Map 4: Coefficient: Index Multiple Deprivation
tm_shape(msoa_result) + 
  tm_fill("localR2", title = "Adaptive: Local R2", style = "cont", midpoint = NA, palette = "magma") +
  tm_polygons(alpha = 0, border.alpha = 0.5, border.col = "black",  lwd=0.1) +
  tm_compass(position = c("left", "top")) +
  tm_scalebar(position = c("right", "bottom")) +
  tm_layout(frame = FALSE, legend.title.size = 0.5, legend.text.size = 0.5, 
            legend.position=c("right", "top"))

tmap_save(filename = "plots/rsquared.png")


str(msoa_result)






