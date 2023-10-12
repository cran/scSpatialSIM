## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## -----------------------------------------------------------------------------
library(scSpatialSIM)

## ---- creating_spatial_object-------------------------------------------------
#set seed
set.seed(333)
#create the new object
sim_object = CreateSimulationObject(sims = 5, cell_types = 1) %>%
    #produce the point pattern
    GenerateSpatialPattern()
    #make tissues
sim_object = GenerateTissue(sim_object, density_heatmap = F) %>%
    #create holes
    GenerateHoles(hole_prob = c(0.3, 0.5), density_heatmap = F) %>%
    #Create positive/negative cells
    GenerateCellPositivity(probs = c(0, 0.9))

## ----sim_object_export--------------------------------------------------------
#creating the spatial list
spatial_list = CreateSpatialList(sim_object, single_df = F)
#summarise the spatial list
summary_df = SummariseSpatial(spatial_list, markers = "Cell 1 Assignment")
head(summary_df)

## -----------------------------------------------------------------------------
library(spatialTIME)

## ----adding_spatial_df_name---------------------------------------------------
#loop over all spatial data frames and add their names
sf_names = names(spatial_list)
spatial_list = lapply(setNames(sf_names, sf_names), function(nam){
  spatial_list[[nam]] %>%
    dplyr::mutate(`Sample Tag` = nam, .before = 1)
})

## -----------------------------------------------------------------------------
summary_df$`Patient ID` = 1:5

## -----------------------------------------------------------------------------
clinical = data.frame(`Patient ID` = 1:5, check.names = F)

## -----------------------------------------------------------------------------
mif = create_mif(clinical_data = clinical,
                 sample_data = summary_df,
                 spatial_list = spatial_list,
                 patient_id = "Patient ID",
                 sample_id = "Sample Tag")
mif

## -----------------------------------------------------------------------------
mif = ripleys_k(mif = mif, mnames = "Cell 1 Assignment", r_range = seq(0, 3, 0.1), 
                permute = FALSE, edge_correction = "translation", workers = 1,
                xloc = "x", yloc = "y")

## -----------------------------------------------------------------------------
library(ggplot2)
mif$derived$univariate_Count %>%
  ggplot() +
  geom_line(aes(x = r, y = `Degree of Clustering Exact`, color = `Sample Tag`)) +
  labs(title = "Univariate Clustering - Simulation")

## -----------------------------------------------------------------------------
mif_holes = create_mif(clinical_data = clinical,
                 sample_data = summary_df,
                 spatial_list = lapply(spatial_list, function(spat){
                   spat %>%
                     dplyr::filter(`Hole Assignment` == "Keep")
                 }),
                 patient_id = "Patient ID",
                 sample_id = "Sample Tag")
mif_holes = ripleys_k(mif = mif_holes, mnames = "Cell 1 Assignment", r_range = seq(0, 3, 0.1), 
                permute = FALSE, edge_correction = "translation", workers = 1,
                xloc = "x", yloc = "y")

## -----------------------------------------------------------------------------
mif_holes$spatial %>%
  do.call(dplyr::bind_rows, .) %>%
  ggplot() +
  geom_point(aes(x = x, y = y, color = as.factor(`Cell 1 Assignment`))) +
  facet_wrap(~`Sample Tag`)

## -----------------------------------------------------------------------------
dat = do.call(dplyr::bind_rows, mif_holes$derived$univariate_Count)
dat %>%
  dplyr::mutate(`Exact-Theo` = `Exact CSR` - `Theoretical CSR`) %>%
  ggplot() +
  geom_density(aes(x = `Exact-Theo`, fill = `Sample Tag`), alpha = 0.2, adjust = 0.2)

## -----------------------------------------------------------------------------
dat %>%
    dplyr::mutate(`Exact-Theo` = `Exact CSR` - `Theoretical CSR`) %>%
    dplyr::group_by(`Sample Tag`) %>%
    dplyr::mutate(prop = ifelse(`Exact-Theo` > 0, 1/dplyr::n(), 0)) %>%
    dplyr::select(`Sample Tag`, r, `Exact-Theo`, prop) %>%
    dplyr::summarise(`Total Fraction` = sum(prop))

