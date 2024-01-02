Steps

1.

Missing from the complete project are LSOA and MSOA shapefiles which can be found here:

UK Data Service. (2022a). ‘2011 Census Geography Boundaries Lower Layour Super Output Areas and Data Zones’. Available From:
https://statistics.ukdataservice.ac.uk/dataset/2011-census-geography-boundaries-lower-layer-super-output-areas-and-data-zones [Accessed 10 December 2023].

UK Data Service. (2022b). ‘2011 Census Geography Boundaries Middle Layour Super Output Areas and Intermediate Zones’. Available From: https://statistics.ukdataservice.ac.uk/dataset/2011-census-geography-boundaries-middle-layer-super-output-areas-and-intermediate-zones [Accessed 10 December 2023].

2. Install packages
library("tmap")
library("sf")
library("spdep")
library("sp")
library("spatialreg")
library("sfdep")
library("spgwr")
library("dplyr")
library("car")

3. Run data_processing.R

4. Run regression.R
