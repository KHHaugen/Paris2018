### Global Moran's I paper
# gwr model based on code from http://rspatial.org/analysis/rst/6-local_regression.html

# Some issue with identical points. Identify where and perhapps make an offset.

## Oslo kommune Språksenteret, orgnr 813993422, have clear issues regarding class size

library(pacman)
p_load(sp, spdep)

data <- read.csv('file:///C:/Users/kha168/Desktop/RURED/data/Grunnskoler/grunnstruktur med koordinater.csv', encoding = 'UTF-8')

# Removing schools located outside the country
data$lon[data$lat < 57.5] <- NA
data$lat[data$lat < 57.5] <- NA

data$lat[(data$lat < 66 | data$lat > 72) & data$lon > 15] <- NA
data$lon[is.na(data$lat)] <- NA

data <- data[!is.na(data$lat), ]

# Removing Språksenteret and schools with reported 0 students
data <- data[data$Org.nr != 813993422,]
data <- data[data$grpst2 != 0,]

# Aggregating the data
data <- aggregate(data[,c('grpoeng', 'grpst2', 'utd', 'pop', 'eleverN', 'lon', 'lat')], list(data$Org.nr, data$Navn), mean, na.rm=T)

data <- data[!is.na(data$utd),]
data <- data[!is.na(data$grpst2),]
#data <- data[!is.na(data$andel_syss),]
normsize <- 20

data$id <- 1:nrow(data)
data$norm <- 0
data$norm[data$grpst2 > normsize] <- 1

colnames(data)[colnames(data) == 'utd'] <- 'edu'

# Creating a variable for expected change in class size after norm


data$delta <- 0
data$delta[data$grpst2 > normsize] <- data$grpst2[data$grpst2 > normsize] - normsize

# Defining coordinates and CRS
sp_point <- cbind(data$lon, data$lat)
proj <- CRS("+init=epsg:4326")
#data.sp <- SpatialPointsDataFrame(coords=sp_point,data,proj4string=proj)
#map_crd <- coordinates(data.sp)

# Finding duplicated geolocations:
#zerodist(data.sp)

# Setting an offset to one in each pair
sp_point[,2][duplicated(sp_point)] <- sp_point[,2][duplicated(sp_point)] + 0.00001

data.sp <- SpatialPointsDataFrame(coords=sp_point,data,proj4string=proj)
data.sp <- spTransform(data.sp, CRS("+init=epsg:25833"))
map_crd <- coordinates(data.sp)
rownames(map_crd) <- data.sp$id

## Defining number of neighbours for 

k <- 1

# Creating a neighbor-list
W_knn1 <- knn2nb(knearneigh(map_crd, k=k))
W_knn1_mat <- nb2listw(W_knn1, style = 'W', zero.policy = TRUE)
# k = 3 
W_knn3 <- knn2nb(knearneigh(map_crd, k=3))
W_knn3_mat <- nb2listw(W_knn3, style = 'W', zero.policy = TRUE)
# k = 5
W_knn5 <- knn2nb(knearneigh(map_crd, k=5))
W_knn5_mat <- nb2listw(W_knn5, style = 'W', zero.policy = TRUE)


# Creating connectivity matrix based on inverted distance
dist <- unlist(nbdists(W_knn1, map_crd)) 
distance <- dnearneigh(map_crd, d1=0, d2=max(dist))
distance.neigh.list <- nbdists(distance, coordinates(map_crd))
W_inv.dist <- lapply(distance.neigh.list, function(x) 1/x)

W_inv.distance <- nb2listw(distance, glist=W_inv.dist, style="W") 

# Displaying the connected schools
plot(W_knn1_mat,coords=map_crd,pch=19, cex=0.1, col="gray")

# Testing the global Moran's I
moran.test(data$grpoeng, listw = W_knn1_mat)

moran_test <- (moran.test(data$grpoeng, listw = W_knn1_mat))

moran_test_mod <- moran_test$estimate
moran_test_mod$'standard deviate' <- moran_test$statistic
moran_test_mod$p.value <- moran_test$p.value
moran_test_mod <- data.frame(moran_test_mod)
rownames(moran_test_mod) <- ""

print(xtable::xtable(moran_test_mod, label = 'moran', caption= "Global Moran's I test"), file = 'Moran_test.tex', table.placement = 'H')

# Plotting GPA against the spatially lagged GPA
moran.plot(data$grpoeng, listw = W_knn1_mat)

# Cook's distance: (influence of each point)
summary(moran.plot(data$grpoeng, listw = W_knn1_mat)$infmat[,c('cook.d')])

## Creating Models
mod1 <- lm(grpoeng ~ grpst2 + I(log(pop)) + edu, data = data) # OLS
mod2 <- lm(grpoeng ~ grpst2 + I(log(pop)) + edu, data = data, weights = eleverN) # OLS with weights
mod3 <- errorsarlm(grpoeng ~ grpst2 + I(log(pop)) + edu, data = data, listw = W_knn1_mat) # simultaneous autoregressive model
mod4 <- errorsarlm(grpoeng ~ grpst2 + I(log(pop)) + edu, data = data, listw = W_knn1_mat, weights = eleverN) # simultaneous autoregressive model with weights
mod5 <- lagsarlm(grpoeng ~ grpst2 + I(log(pop)) + edu, data = data, listw = W_knn1_mat) # Spatially lagged dependent variable
mod6 <-  errorsarlm(grpoeng ~ grpst2 + I(log(pop)) + edu, data = data, listw=W_inv.distance, zero.policy=T, weights = eleverN) # Geostatistical simultaneous model
mod7 <- lagsarlm(grpoeng ~ grpst2 + I(log(pop)) + edu, data = data.sp, listw=W_inv.distance, zero.policy=T) # Geostatistical spatially lagged dependent

## Models with k = 3
mod1.3 <- lm(grpoeng ~ grpst2 + I(log(pop)) + edu, data = data) # OLS
mod2.3 <- lm(grpoeng ~ grpst2 + I(log(pop)) + edu, data = data, weights = eleverN) # OLS with weights
mod3.3 <- errorsarlm(grpoeng ~ grpst2 + I(log(pop)) + edu, data = data, listw = W_knn3_mat) # simultaneous autoregressive model
mod4.3 <- errorsarlm(grpoeng ~ grpst2 + I(log(pop)) + edu, data = data, listw = W_knn3_mat, weights = eleverN) # simultaneous autoregressive model with weights
mod5.3 <- lagsarlm(grpoeng ~ grpst2 + I(log(pop)) + edu, data = data, listw = W_knn3_mat) # Spatially lagged dependent variable

moran_test3 <- (moran.test(data$grpoeng, listw = W_knn3_mat))

moran_test_mod3 <- moran_test3$estimate
moran_test_mod3$'standard deviate' <- moran_test3$statistic
moran_test_mod3$p.value <- moran_test3$p.value
moran_test_mod3 <- data.frame(moran_test_mod3)
rownames(moran_test_mod3) <- ""

print(xtable::xtable(moran_test_mod3, label = 'moran3', caption= "Global Moran's I test"), file = 'Moran_test3.tex', table.placement = 'H')

## Models with k = 5
mod1.5 <- lm(grpoeng ~ grpst2 + I(log(pop)) + edu, data = data) # OLS
mod2.5 <- lm(grpoeng ~ grpst2 + I(log(pop)) + edu, data = data, weights = eleverN) # OLS with weights
mod3.5 <- errorsarlm(grpoeng ~ grpst2 + I(log(pop)) + edu, data = data, listw = W_knn5_mat) # simultaneous autoregressive model
mod4.5 <- errorsarlm(grpoeng ~ grpst2 + I(log(pop)) + edu, data = data, listw = W_knn5_mat, weights = eleverN) # simultaneous autoregressive model with weights
mod5.5 <- lagsarlm(grpoeng ~ grpst2 + I(log(pop)) + edu, data = data, listw = W_knn5_mat) # Spatially lagged dependent variable

moran_test5 <- (moran.test(data$grpoeng, listw = W_knn5_mat))

moran_test_mod5 <- moran_test5$estimate
moran_test_mod5$'standard deviate' <- moran_test5$statistic
moran_test_mod5$p.value <- moran_test5$p.value
moran_test_mod5 <- data.frame(moran_test_mod5)
rownames(moran_test_mod5) <- ""

print(xtable::xtable(moran_test_mod5, label = 'moran5', caption= "Global Moran's I test"), file = 'Moran_test5.tex', table.placement = 'H')

## Moran test geostatistical models
moran_test_geo <- (moran.test(data$grpoeng, listw = W_inv.distance))

moran_test_mod_geo <- moran_test_geo$estimate
moran_test_mod_geo$'standard deviate' <- moran_test_geo$statistic
moran_test_mod_geo$p.value <- moran_test_geo$p.value
moran_test_mod_geo <- data.frame(moran_test_mod_geo)
rownames(moran_test_mod_geo) <- ""

print(xtable::xtable(moran_test_mod_geo, label = 'moranGeo', caption= "Global Moran's I test"), file = 'Moran_test_geo.tex', table.placement = 'H')

## Combining Moran I's from the model specifications

morancomb <- rbind(moran_test_mod, moran_test_mod3, moran_test_mod5, moran_test_mod_geo)
rownames(morancomb) <- c('k = 1', 'k = 3', 'k = 5', 'Inverse distance')

print(xtable::xtable(morancomb, label = 'morancomb', caption= "Global Moran's I test",digits=c(1,2,3,4,2,4)), file = 'Moran_test_comb.tex', table.placement = 'H')

## Testing for spatial autocorrelation in OLS models
lm.morantest(mod1, W_knn1_mat)
lm.morantest(mod2, W_knn1_mat)

## Displaying Models
summary(mod1)
summary(mod2)
summary(mod3)
summary(mod4)
summary(mod5)
summary(mod6)
summary(mod7)


# Local regression
library(spgwr)

#bw <- gwr.sel(grpoeng ~grpst2, data=data.sp)
bw <- gwr.sel(grpoeng ~ grpst2 + I(log(pop)) + edu, data = data.sp, weights = eleverN)

mod9 <- gwr(grpoeng ~ grpst2 + I(log(pop)) + edu, data = data.sp, weights = eleverN, bandwidth = bw)
mod9

gwrMod <- rbind(intercept = as.numeric(summary(mod9$SDF$'(Intercept)')), grst2 = as.numeric(summary(mod9$SDF$grpst2)),
                     logpop = as.numeric(summary(mod9$SDF$'I(log(pop))')), edu = as.numeric(summary(mod9$SDF$edu)))
colnames(gwrMod) <- names(summary(mod9$SDF$grpst2))


print(xtable::xtable(gwrMod, caption = 'Summary table of the GWR parameters', label = 'gwrparams'), file = 'gwrMod.tex', table.placement = 'H')

gwrRes <- data.frame(lon = data.sp$lon, lat = data.sp$lat, grpoeng = data.sp$grpoeng, delta = data.sp$delta, skole = data.sp$Group.2)
gwrRes$b.grpst <- mod9$SDF$grpst2
gwrRes$b.int <- mod9$SDF$'(Intercept)'


gwrRes$exp.change <- gwrRes$delta*gwrRes$b.grpst
gwrRes$y <- mod9$lm$y

summary(gwrRes$exp.change)

# Moran's I after norm using GWR
gwrRes$postnorm <- gwrRes$y + (gwrRes$delta*gwrRes$b.grpst)

moran_test_post_gwr <- (moran.test(gwrRes$postnorm, listw = W_knn1_mat))

moran_test_mod_post_gwr <- moran_test_post_gwr$estimate
moran_test_mod_post_gwr$'standard deviate' <- moran_test_post_gwr$statistic
moran_test_mod_post_gwr$p.value <- moran_test_post_gwr$p.value
moran_test_mod_post_gwr <- data.frame(moran_test_mod_post_gwr)
rownames(moran_test_mod_post_gwr) <- ""

moran_test_post3_gwr <- (moran.test(gwrRes$postnorm, listw = W_knn3_mat))

moran_test_mod_post3_gwr <- moran_test_post3_gwr$estimate
moran_test_mod_post3_gwr$'standard deviate' <- moran_test_post3_gwr$statistic
moran_test_mod_post3_gwr$p.value <- moran_test_post3_gwr$p.value
moran_test_mod_post3_gwr <- data.frame(moran_test_mod_post3_gwr)
rownames(moran_test_mod_post3_gwr) <- ""

moran_test_post5_gwr <- (moran.test(gwrRes$postnorm, listw = W_knn5_mat))

moran_test_mod_post5_gwr <- moran_test_post5_gwr$estimate
moran_test_mod_post5_gwr$'standard deviate' <- moran_test_post5_gwr$statistic
moran_test_mod_post5_gwr$p.value <- moran_test_post5_gwr$p.value
moran_test_mod_post5_gwr <- data.frame(moran_test_mod_post5_gwr)
rownames(moran_test_mod_post5_gwr) <- ""

moran_test_post_geo_gwr <- (moran.test(gwrRes$postnorm, listw = W_inv.distance))

moran_test_mod_post_geo_gwr <- moran_test_post_geo_gwr$estimate
moran_test_mod_post_geo_gwr$'standard deviate' <- moran_test_post_geo_gwr$statistic
moran_test_mod_post_geo_gwr$p.value <- moran_test_post_geo_gwr$p.value
moran_test_mod_post_geo_gwr <- data.frame(moran_test_mod_post_geo_gwr)
rownames(moran_test_mod_post_geo_gwr) <- ""

morancomb_post_gwr <- rbind(moran_test_mod_post_gwr, moran_test_mod_post3_gwr, moran_test_mod_post5_gwr, moran_test_mod_post_geo_gwr)
rownames(morancomb_post_gwr) <- c('k = 1', 'k = 3', 'k = 5', 'Inverse distance')

print(xtable::xtable(morancomb_post_gwr, label = 'morancombpostgwr', caption= "Global Moran's I test",digits=c(1,2,3,4,2,4)), file = 'Moran_test_comb_post_gwr.tex', table.placement = '!')


# Seems like the norm will, on average, make the schools directly effected perform worse than before the norm 

library(ggplot2)

data.sp$percentile<-(rank(data.sp$grpoeng)/length(data.sp$grpoeng))*100
#mid <- mean(data.sp$grpoeng)

GPA <- ggplot(as.data.frame(data.sp), aes(map_crd[,1], map_crd[,2])) +
  geom_point(aes(colour = percentile),
             size = 0.5) +
  scale_color_gradient(low="red", high="blue") +
  labs(title = 'GPA of schools', x = 'lon', y = 'lat') +
  theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_blank())

GPA
ggsave("GPA.png", width = 4.5, height = 5.2)





kart <- ggplot(gwrRes, aes(map_crd[,1], map_crd[,2])) +
geom_point(aes(colour = b.grpst),
           size = 0.5) +scale_color_gradient(low="blue", high="red") +
  labs(title = 'Local regression parameter value of group size', x = 'lon', y = 'lat') +
  theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_blank())

kart
ggsave("kart.png", width = 4.5, height = 5.2)



# Colors:
colnorm <- c('1' = 'red', '0' = 'green')

norm <- ggplot(as.data.frame(data.sp), aes(map_crd[,1], map_crd[,2])) +
  geom_point(aes(colour = factor(norm)), size = 0.5, alpha = 0.7) +
  scale_colour_manual(values = colnorm, label = c('Not directly impacted', 'Directly impacted'), name = '') +
  labs(title = '', x = 'lon', y = 'lat') +
  theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_blank())


norm
ggsave("norm.png", width = 5.2, height = 5.2)


# Output

texreg::texreg(list(mod1, mod2, mod3, mod4, mod5), file = 'table.tex',
               custom.model.names = c('OLS','Weighted OLS', 'errorsar', 'weighted errorsar', 'lagsar'),
               caption = 'Regression models using k = 1 nearest neighbors', label = 'regtable1', float.pos = 'h!')

texreg::texreg(list(mod1.3, mod2.3, mod3.3, mod4.3, mod5.3), file = 'table_k_3.tex',
               custom.model.names = c('OLS','Weighted OLS', 'errorsar', 'weighted errorsar', 'lagsar'),
               caption = 'Regression models using k = 3 nearest neighbors', label = 'regtable3', float.pos = 'h!')

texreg::texreg(list(mod1.5, mod2.5, mod3.5, mod4.5, mod5.5, mod6, mod7), file = 'table_k_5.tex',
               custom.model.names = c('OLS','Weighted OLS', 'errorsar', 'weighted errorsar', 'lagsar', 'weighted errorsar geo', 'lagsar geo'),
               caption = 'Regression models using k = 5 nearest neighbors', label = 'regtable5', float.pos = 'h!')



#### Expected change
summary(data$grpoeng[data$norm ==1]) 
summary(data$grpoeng[data$norm ==0])
# Can see that the schools directly effected by the norm already have higher GPA than the schools that won't be directly effected
data.sp$postnorm <- data.sp$grpoeng + (data.sp$grpoeng*data.sp$delta*mod4.5$coefficients[2])

moran_test_post <- (moran.test(data.sp$postnorm, listw = W_knn1_mat))

moran_test_mod_post <- moran_test_post$estimate
moran_test_mod_post$'standard deviate' <- moran_test_post$statistic
moran_test_mod_post$p.value <- moran_test_post$p.value
moran_test_mod_post <- data.frame(moran_test_mod_post)
rownames(moran_test_mod_post) <- ""

moran_test_post3 <- (moran.test(data.sp$postnorm, listw = W_knn3_mat))

moran_test_mod_post3 <- moran_test_post3$estimate
moran_test_mod_post3$'standard deviate' <- moran_test_post3$statistic
moran_test_mod_post3$p.value <- moran_test_post3$p.value
moran_test_mod_post3 <- data.frame(moran_test_mod_post3)
rownames(moran_test_mod_post3) <- ""

moran_test_post5 <- (moran.test(data.sp$postnorm, listw = W_knn5_mat))

moran_test_mod_post5 <- moran_test_post5$estimate
moran_test_mod_post5$'standard deviate' <- moran_test_post5$statistic
moran_test_mod_post5$p.value <- moran_test_post5$p.value
moran_test_mod_post5 <- data.frame(moran_test_mod_post5)
rownames(moran_test_mod_post5) <- ""

moran_test_post_geo <- (moran.test(data.sp$postnorm, listw = W_inv.distance))

moran_test_mod_post_geo <- moran_test_post_geo$estimate
moran_test_mod_post_geo$'standard deviate' <- moran_test_post_geo$statistic
moran_test_mod_post_geo$p.value <- moran_test_post_geo$p.value
moran_test_mod_post_geo <- data.frame(moran_test_mod_post_geo)
rownames(moran_test_mod_post_geo) <- ""

morancomb_post <- rbind(moran_test_mod_post, moran_test_mod_post3, moran_test_mod_post5, moran_test_mod_post_geo)
rownames(morancomb_post) <- c('k = 1', 'k = 3', 'k = 5', 'Inverse distance')

print(xtable::xtable(morancomb_post, label = 'morancombpost', caption= "Global Moran's I test",digits=c(1,2,3,4,2,4)), file = 'Moran_test_comb_post.tex', table.placement = 'H')

# Summary statistics
summarystat <- rbind(grpoeng = as.numeric(summary(data.sp$grpoeng)), grst2 = as.numeric(summary(data.sp$grpst2)),
                pop = summary(as.numeric(data.sp$pop)), logpop = as.numeric(summary(log(data.sp$pop))), edu = as.numeric(summary(data.sp$edu)),
                eleverN = as.numeric(summary(data.sp$eleverN)))
colnames(summarystat) <- names(summary(data.sp$grpst2))


print(xtable::xtable(summarystat, caption = 'Summary statistics', label = 'sumstat'), file = 'sumstat.tex', table.placement = 'H')

GPAcomp <- rbind(not = as.numeric(summary(data.sp$grpoeng[data.sp$grpst2 <= normsize])), are = as.numeric(summary(data.sp$grpoeng[data.sp$grpst2 > normsize])))
colnames(GPAcomp) <- names(summary(data.sp$grpoeng[data.sp$grpst2 <= normsize]))
rownames(GPAcomp) <- c('Not directly impacted', 'Directly impacted')
print(xtable::xtable(GPAcomp, caption = 'GPA by impact of norm', label = 'GPAcomp'), file = 'GPAcomp.tex', table.placement = 'H')

# t test of Moran's I:

t_test_var <- function(mu1, mu2, var1, var2, n1, n2){
  print((mu1-mu2)/sqrt((var1/n1)+(var2/n2)))
}


morancomb_test <- morancomb[,c(1,2)]
morancomb_test[,c(3:4)] <- morancomb_post[,c(1,2)]
morancomb_test[,5] <-t_test_var(mu1 = morancomb[,1], mu2 = morancomb_post[,1], var1 = morancomb[,3], var2 = morancomb_post[,3], n1 = 909, n2 = 909)
morancomb_test <- morancomb_test[,c(1,3,5)]
colnames(morancomb_test) <- c('MI pre', 'MI post', 't statistic')

morancomb_test

print(xtable::xtable(morancomb_test, caption = "Moran's I comparison k = 5", label = 'MIcomp'), file = 'MIcomp.tex', table.placement = 'H')

morancomb_test_gwr <- morancomb[,c(1,2)]
morancomb_test_gwr[,c(3:4)] <- morancomb_post_gwr[,c(1,2)]
morancomb_test_gwr[,5] <-t_test_var(mu1 = morancomb[,1], mu2 = morancomb_post_gwr[,1], var1 = morancomb[,3], var2 = morancomb_post_gwr[,3], n1 = 909, n2 = 909)
morancomb_test_gwr <- morancomb_test_gwr[,c(1,3,5)]
colnames(morancomb_test_gwr) <- c('MI pre', 'MI post', 't statistic')

morancomb_test_gwr

print(xtable::xtable(morancomb_test_gwr, caption = "Moran's I comparison GWR", label = 'MIcompgwr'), file = 'MIcompGWR.tex', table.placement = 'H')


groupcomp <- t.test(data.sp$grpoeng[data.sp$norm == 1], data.sp$grpoeng[data.sp$norm == 0])

groupcomptest <- cbind(groupcomp$estimate[1], groupcomp$estimate[2], groupcomp$conf.int[1], groupcomp$conf.int[2], groupcomp$statistic, groupcomp$p.value)
colnames(groupcomptest) <- c('Mean of x', 'Mean of y', '2.5% confint bound', '97.5% confint bound', 't value', 'p value')
rownames(groupcomptest) <- ' '

print(xtable::xtable(groupcomptest, caption = "T test of GPA by norm impact", label = 'tgroup'), file = 'tgroup.tex', table.placement = 'H')
