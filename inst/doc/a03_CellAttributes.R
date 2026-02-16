## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
library(ggpubr)
library(dplyr)

## -----------------------------------------------------------------------------
library(scSpatialSIM)
set.seed(333)
#create simulation object for tight clusters
sim_object_tight = CreateSimulationObject(sims = 5, cell_types = 1) %>%
  #simulate point pattern
  GenerateSpatialPattern()
#create simulation object for no clusters
sim_object_null = CreateSimulationObject(sims = 5, cell_types = 1) %>%
  #simulate point pattern
  GenerateSpatialPattern()

## -----------------------------------------------------------------------------
sim_object_tight = GenerateCellPositivity(sim_object_tight, 
                                          k = 4,
                                          sdmin = 1, sdmax = 3,
                                          density_heatmap = F,
                                          probs = c(0.0, 0.9))

sim_object_null = GenerateCellPositivity(sim_object_null, 
                                          probs = c(0.0, 0.2),
                                          no_kernel = TRUE)

## -----------------------------------------------------------------------------
library(ggplot2)
#calculate abundance for clustered samples
cluster_abundance = sapply(sim_object_tight@`Spatial Files`, function(x){
  sum(x$`Cell 1 Assignment` == "Positive")/nrow(x)
})
#calculate abundance for null/negative control samples
null_abundance = sapply(sim_object_null@`Spatial Files`, function(x){
  sum(x$`Cell 1 Assignment` == "Positive")/nrow(x)
})
#create histogram of abundances
data.frame(abundance = c(cluster_abundance, null_abundance),
           simulation = c(rep("clustered", 5), rep("null", 5))) %>%
  ggplot() + 
  geom_histogram(aes(x = abundance, fill = simulation))

## -----------------------------------------------------------------------------
#extract simulated samples to make lists
cluster_list = CreateSpatialList(sim_object_tight, single_df = FALSE)
null_list = CreateSpatialList(sim_object_null, single_df = FALSE)
#combine lists to make a single list
spatial_list = c(cluster_list, null_list)
names(spatial_list) = c(paste0("Clustered Sample ", 1:5),
                        paste0("Null Sample ", 1:5))

## -----------------------------------------------------------------------------
spat_data_distribution = GenerateDistributions(spatial_list, 
                                               positive_mean = 15,
                                               negative_mean = 5,
                                               positive_sd = 3,
                                               negative_sd = 1)

## -----------------------------------------------------------------------------
#identify neighbors
#create weight list
#calculate moran's i
library(dplyr)
library(spdep)
results = lapply(spat_data_distribution, function(dat){
  #convert data frame to an sf object compatable with other spdep functions
  sf_dat = st_as_sf(dat, coords = c("x", "y"))
  #calculate the 10 nearest neighbors
  knn = knearneigh(sf_dat, k = 10)
  #convert knn to neighbor list
  knn_nb = knn2nb(knn)
  #convert neighbor list to weight list
  knn_nb_listw = nb2listw(knn_nb,
                          style = "W")
  #calculate moran's I on the simulated protein expression "Cell 1 Var"
  res = moran.test(sf_dat$`Cell 1 Var`,
                   listw = knn_nb_listw)
  #convert results to a data frame
  data.frame(as.list(res$estimate), check.names = FALSE)
}) %>%
  #convert list results to a data frame
  bind_rows(.id = "Sample ID")
#look at the results
head(results)

## -----------------------------------------------------------------------------
results2 = results %>%
    mutate(Group = rep(c("Clustered", "Null"), each = 5))
results2 %>%
  ggplot() +
  geom_boxplot(aes(x = Group, y = `Moran I statistic`)) +
  theme_classic()

