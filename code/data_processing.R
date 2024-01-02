age_profile <- read.csv('data/aggregated_age_profile.csv')
age_counts <- age_profile %>%
  select(matches("^X[0-9][0-9]"))
boolean_condition <- as.numeric(str_extract(colnames(age_counts), "(?<=X)\\d+")) > 65
age_counts <- age_counts %>%
  select(all_of(names(age_counts)[boolean_condition]))
age_profile$older_than_65 <- rowSums(age_counts, na.rm = TRUE)
age_profile <- age_profile %>%
  select("Area.Codes", "Year", "All.Ages", "older_than_65")
age_profile$percent_over_65 <- age_profile$older_than_65/age_profile$All.Ages

age_profile <- age_profile %>%
  select("Area.Codes", "All.Ages", "Year", "percent_over_65")

imd <- read.csv('data/aggregated_imd.csv')
imd <- imd %>% 
  select("FeatureCode", "Value")

imd <- imd %>%
  rename(imd_decile = Value)

data <- merge(age_profile, imd, by.x="Area.Codes", by.y="FeatureCode")

opioid <- read.csv('data/aggregated_opioid_prescribing.csv')
opioid <- opioid %>%
  select("lsoa11", "year", "items_r")
opioid <- opioid %>%
  rename(opioid_prescribing_rate = items_r)

data <- merge(data, opioid, by.x=c("Area.Codes", "Year"), by.y=c("lsoa11", "year"))

samhi <- read.csv('data/aggregated_samhi.csv')

data <- merge(data, samhi, by.x=c("Area.Codes", "Year"), by.y=c("lsoa11", "year"))

data <- subset(data, Year == 2017)

uk_lsoa <- read_sf('data/infuse_lsoa_lyr_2011/infuse_lsoa_lyr_2011.shp', crs="OSGB36")
uk_msoa <- read_sf('data/infuse_msoa_lyr_2011/infuse_msoa_lyr_2011.shp', crs="OSGB36")

joined_soa <- st_join(uk_lsoa, uk_msoa, join = st_within)

joined_soa <- joined_soa[c("geo_code.x", "name.x", "geo_code.y", "name.y")]

data2 <- merge(data, joined_soa, by.x="Area.Codes", by.y="geo_code.x")

data2 <- data2[c("Area.Codes", "geo_code.y", "All.Ages", "name.y", "percent_over_65",
               "imd_decile", "opioid_prescribing_rate", "samhi_index")]

grouped_data <- data2 %>%
  group_by(geo_code.y, name.y) %>%
  summarize(
    mean_percent_over_65 = weighted.mean(percent_over_65, w = All.Ages),
    mean_imd_decile = weighted.mean(imd_decile, w = All.Ages),
    mean_opioid_prescribing_rate = weighted.mean(opioid_prescribing_rate, w = All.Ages),
    mean_samhi_index = weighted.mean(samhi_index, w = All.Ages),
  )

grouped_data <- grouped_data %>%
  rename(
    msoa_geo_code = geo_code.y,
    msoa_name = name.y
  )

rural <- read.csv('data/rural_urban_classification.csv')
letter_to_number <- function(letter) {
  letters <- LETTERS
  if (letter %in% letters) {
    return(match(letter, letters))
  } 
}

# Apply the mapping function to the column
rural$ru_cat <- sapply(strsplit(rural$RUC11CD, ""), function(x) {
  letter_to_number(x[1])
})

rural <- rural %>%
  select("MSOA11CD", "ru_cat")

grouped_data <- merge(grouped_data, rural, by.x=c("msoa_geo_code"), by.y=c("MSOA11CD"))

write.csv(grouped_data, 'data/merge_grouped_msoa.csv')


