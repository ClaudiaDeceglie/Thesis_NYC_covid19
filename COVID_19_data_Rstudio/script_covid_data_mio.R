#Version_Claudia : CRAN packages:
library(tidyverse)
library(convertr)
library(sp)
library(sf)
library(lubridate)
library(tidycensus)
library(ggExtra)
library(ggridges)
library(rstan)
library(drc)
library(spdep)
library(mgcv)
library(broom)
library(MASS)
library(spatialreg)
library(here)
library(pdftools)
library(matrixStats)
library(egg)
library(ggpubr)
library(scales)
library(XML)
library(rgdal)

#### SESSION CONFIGURATIONS ####
here() # current working directory
# if not already set via an environment variable, cache MTA turnstile data within the working directory
if(Sys.getenv("MTA_TURNSTILE_DATA_DIR") == ""){
  mta_dir = here("data/mta_turnstile")
  if(!dir.exists(mta_dir)) dir.create(mta_dir, recursive = TRUE)
  Sys.setenv(MTA_TURNSTILE_DATA_DIR = mta_dir)
} 
## R packages available from GitHub respositories via: 
#remotes::install_github("justlab/Just_universal", ref = "78812f519da11502706a5061e7b8bc4812e5c3b5") 
#remotes::install_github("justlab/MTA_turnstile", ref = "6c8bd7690dfa6036bf991cb4504f42631e8f6756")
library(Just.universal) 
library(MTA.turnstile) # if not already set, this will process MTA data in a temporary directory

if(!dir.exists(here("figures"))) dir.create(here("figures"))
Sys.time() # print the start time

##To generate census data, you must set an API key, which you can request here: https://api.census.gov/data/key_signup.html
census_api_key("9f161986e1eaed4ab7930d1cb16cbcf4fde050ba", install = TRUE, overwrite=TRUE) 
if(Sys.getenv("CENSUS_API_KEY")=="9f161986e1eaed4ab7930d1cb16cbcf4fde050ba.") stop("Census API Key Missing. Please see ?census_api_key()")

export.figs = FALSE #change to TRUE if you would like to save out figures 

# data will default to a subfolder "data/" within working directory
# unless 1. set by an environment variable:
data.root = Sys.getenv("COVID_DATA")
# or 2. set with an alternative path here:
if (data.root == "") data.root = here("data")
if (data.root == "data" & !dir.exists(data.root)) dir.create(here("data"))
print(paste("data being downloaded into directory", dQuote(data.root)))


##### FUNCTIONS ####

read_w_filenames <- function(flnm) {
  read_csv(flnm) %>%
    mutate(filename = flnm)
}

extract_waic <- function (stanfit){
  log_lik <- rstan::extract(stanfit, "log_lik")$log_lik
  dim(log_lik) <- if (length(dim(log_lik)) == 1) 
    c(length(log_lik), 1)
  else c(dim(log_lik)[1], prod(dim(log_lik)[2:length(dim(log_lik))]))
  S <- nrow(log_lik)
  n <- ncol(log_lik)
  lpd <- log(colMeans(exp(log_lik)))
  p_waic <- colVars(log_lik)
  elpd_waic <- lpd - p_waic
  waic <- -2 * elpd_waic
  loo_weights_raw <- 1/exp(log_lik - max(log_lik))
  loo_weights_normalized <- loo_weights_raw/matrix(colMeans(loo_weights_raw), 
                                                   nrow = S, ncol = n, byrow = TRUE)
  loo_weights_regularized <- pmin(loo_weights_normalized, sqrt(S))
  elpd_loo <- log(colMeans(exp(log_lik) * loo_weights_regularized)/colMeans(loo_weights_regularized))
  p_loo <- lpd - elpd_loo
  pointwise <- cbind(waic, lpd, p_waic, elpd_waic, p_loo, elpd_loo)
  total <- colSums(pointwise)
  se <- sqrt(n * colVars(pointwise))
  return(list(waic = total["waic"], elpd_waic = total["elpd_waic"], 
              p_waic = total["p_waic"], elpd_loo = total["elpd_loo"], 
              p_loo = total["p_loo"]))
}

quantile_split <- function (data, mix_name = mix_name, q, shift = TRUE) {
  if (shift) {
    for (i in 1:length(mix_name)) {
      dat_num = as.numeric(unlist(data[, mix_name[i]]))
      data[[mix_name[i]]] = cut(dat_num, breaks = unique(quantile(dat_num, 
                                                                  probs = seq(0, 1, by = 1/q), na.rm = TRUE)), 
                                labels = FALSE, include.lowest = TRUE) - 1
    }
  }
  else {
    for (i in 1:length(mix_name)) {
      dat_num = as.numeric(unlist(data[, mix_name[i]]))
      data[[mix_name[i]]] = cut(dat_num, breaks = unique(quantile(dat_num, 
                                                                  probs = seq(0, 1, by = 1/q), na.rm = TRUE)), 
                                labels = FALSE, include.lowest = TRUE)
    }
  }
  return(data)
}

download = function(url, to, f, ...){
  download.update.meta(data.root, url, to, f, ...)
}


##### Load Data #####

# get the Pluto dataset from #https://www1.nyc.gov/site/planning/data-maps/open-data/dwn-pluto-mappluto.page 
Pluto <- read_csv('/home/claudia/Scrivania/Tesi NYC pandemic /COVID_19_data_Rstudio/code/https:/www1.nyc.gov/assets/planning/download/zip/data-maps/open-data/nyc_pluto_20v6_csv/pluto_20v6.csv',
                     col_types = cols(spdist2 = col_character(),
                     overlay2 = col_character(),
                     zonedist4 = col_character()))

#Bldg_Footprints<-download(
#  "/home/claudia/Scrivania/Tesi NYC pandemic /COVID_19_data_Rstudio/code/https:/data.cityofnewyork.us/api/geospatial/nqwf-w8eh?method=export&format=Shapefile/",
#  "Building Footprints.zip",
#  function(p)
#    st_read(paste0("/vsizip/", p)))

#buldings
Bldg_Footprints <-st_read("/home/claudia/Scrivania/Tesi NYC pandemic /COVID_19_data_Rstudio/code/https:/data.cityofnewyork.us/api/geospatial/nqwf-w8eh?method=export&format=Shapefile/Building Footprints/geo_export_41a47b75-623c-4ea2-923b-761cca2f06a3.shp")

ZCTA_by_boro<-read_csv("/home/claudia/Scrivania/Tesi NYC pandemic /COVID_19_data_Rstudio/ZCTA_per_borough_unpivot.csv")
ZCTA_by_boro['zcta']<-lapply(ZCTA_by_boro["zcta"], as.character)


#ZCTA_by_boro.html = htmlTreeParse('/home/claudia/Scrivania/Tesi NYC pandemic /COVID_19_data_Rstudio/NYC Neighborhood ZIP Code Definitions.htm',
                         useInternal = TRUE)


ZCTA_test_series <-read_csv("/home/claudia/Scrivania/Tesi NYC pandemic /COVID_19_data_Rstudio/tests-by-zcta.csv")

ZCTAs_in_NYC <- as.character(unique(ZCTA_test_series$modzcta))

#Subway_ridership_by_UHF <- relative.subway.usage(2020L, "nhood")

#UHF definitions by zip codes
UHF_ZipCodes<-read_csv("/home/claudia/Scrivania/Tesi NYC pandemic /COVID_19_data_Rstudio/UHF NEIGHBORHOOD.csv")

#UHF shapefile 
UHF_shp<-st_read("/home/claudia/Scrivania/Tesi NYC pandemic /COVID_19_data_Rstudio/UHF_42_DOHMH_2009/UHF_42_DOHMH_2009.shp")

# NYC boroughs from NYC Open Data
NYC_basemap_shp <-st_read("/home/claudia/Scrivania/Tesi NYC pandemic /COVID_19_data_Rstudio/Borough Boundaries/geo_export_b1379fc0-f3d9-4f42-9d37-10be0cf945eb.shp", stringsAsFactors = FALSE, quiet = TRUE) %>% st_transform(., crs = 2263)

#DOHMH MODZCTA Shapefile
MODZCTA_NYC_shp <-read_sf("/home/claudia/Scrivania/Tesi NYC pandemic /COVID_19_data_Rstudio/Modified Zip Code Tabulation Areas (MODZCTA)/geo_export_efada270-fac3-4177-8416-2084afd869cf.shp")

#Food outlets 
food_retail <- read_csv("/home/claudia/Scrivania/Tesi NYC pandemic /COVID_19_data_Rstudio/Retail_Food_Stores.csv")

# Download deaths by ZCTA as of October 1st
deaths_by1Oct2020_by_zcta<-read_csv("/home/claudia/Scrivania/Tesi NYC pandemic /coronavirus-data-nychealth-october/data-by-modzcta.csv")

# download MODZCTA to ZCTA crosswalk, current version from repo
modzcta_to_zcta <-read_csv("/home/claudia/Scrivania/Tesi NYC pandemic /coronavirus-data-nychealth-october/Geography-resources/ZCTA-to-MODZCTA.csv")

##We have many sources of data, so these just help to combine the various data types
NYC_counties1 <- c("Bronx","Kings","Queens","New York","Richmond")
NYC_counties1_full <- c("Bronx County","Kings County","Queens County","New York County","Richmond County")
NYC_boro_county_match <- tibble(County = c("Bronx","Kings","Queens","New York","Richmond"), 
                                boro = c("Bronx","Brooklyn","Queens","Manhattan","Staten Island"), 
                                full_county = c("Bronx County","Kings County","Queens County","New York County","Richmond County"))

#upload the stan code alongside the disparities code 
BWQS_stan_model <- here("code", "nb_bwqs_cov.stan") 

####Census Data Collection and Cleaning####

ACS_Data <- get_acs(geography = "zcta", 
                    variables = c(medincome = "B19013_001",
                                  total_pop1 = "B01003_001",
                                  fpl_100 = "B06012_002", 
                                  fpl_100to150 = "B06012_003",
                                  median_rent = "B25031_001",
                                  total_hholds1 = "B22003_001",
                                  hholds_snap = "B22003_002",
                                  over16total_industry1 = "C24050_001",
                                  ag_industry = "C24050_002",
                                  construct_industry = "C24050_003",
                                  manufact_industry = "C24050_004",
                                  wholesaletrade_industry = "C24050_005",
                                  retail_industry = "C24050_006",
                                  transpo_and_utilities_industry = "C24050_007",
                                  information_industry = "C24050_008",
                                  finance_and_realestate_industry = "C24050_009",
                                  science_mngmt_admin_industry = "C24050_010",
                                  edu_health_socasst_industry = "C24050_011",
                                  arts_entertain_rec_accomodate_industry = "C24050_012",
                                  othsvcs_industry = "C24050_013",
                                  publicadmin_industry = "C24050_014",
                                  total_commute1 = "B08301_001",
                                  drove_commute = "B08301_002",
                                  pubtrans_bus_commute = "B08301_011",
                                  pubtrans_subway_commute = "B08301_013",
                                  pubtrans_railroad_commute = "B08301_013",
                                  pubtrans_ferry_commute = "B08301_015",
                                  taxi_commute = "B08301_016",
                                  bicycle_commute = "B08301_018",
                                  walked_commute = "B08301_019",
                                  workhome_commute = "B08301_021",
                                  unemployed = "B23025_005",
                                  under19_noinsurance = "B27010_017",
                                  age19_34_noinsurance = "B27010_033",
                                  age35_64_noinsurance = "B27010_050",
                                  age65plus_noinsurance = "B27010_066",
                                  hisplat_raceethnic = "B03002_012",
                                  nonhispLat_white_raceethnic = "B03002_003",
                                  nonhispLat_black_raceethnic = "B03002_004",
                                  nonhispLat_amerindian_raceethnic = "B03002_005",
                                  nonhispLat_asian_raceethnic = "B03002_006",
                                  age65_plus  = "B08101_008"),
                    year = 2018,
                    output = "wide",
                    survey = "acs5")

ACS_Data1 <- ACS_Data %>% #only pull out the estimates and cleaning variable names
  filter(GEOID %in% ZCTAs_in_NYC) %>%
  dplyr::select(-NAME)  %>%
  dplyr::select(GEOID, !ends_with("M")) %>%
  rename_at(vars(ends_with("E")), .funs = list(~str_sub(., end = -2)))

ACS_Data2 <- ACS_Data1 %>%
  mutate_at(vars(ends_with("_commute")), ~round((./total_commute1)*100, 2)) %>% #proportion of people relying on a given mode of transit
  mutate_at(vars(ends_with("_raceethnic")), ~round((./total_pop1)*100, 2)) %>% #proportion of ppl reporting a given race/ethncity 
  mutate(not_insured = round(((under19_noinsurance + age19_34_noinsurance + age35_64_noinsurance + age65plus_noinsurance) / total_pop1)*100, 2), #proportion uninsured
         snap_hholds = round((hholds_snap/total_hholds1)*100, 2), #proportion relying on SNAP
         fpl_150 = round(((fpl_100+fpl_100to150)/total_pop1)*100, 2), #proportion 150% or less of FPL
         unemployed = round((unemployed/over16total_industry1)*100, 2), #proportion unemployed
         not_quarantined_jobs = round(((ag_industry+(construct_industry*.25)+wholesaletrade_industry+ #an estimate of who is still leaving the house for work based on industry
                                          (edu_health_socasst_industry*.5)+transpo_and_utilities_industry)/over16total_industry1)*100, 2)) %>%
  dplyr::select(-ends_with("_noinsurance"), -fpl_100, -fpl_100to150, -ends_with("_industry"), -hholds_snap) %>%
  rename(zcta = "GEOID")

#### Estimating the mode of transportation for essential workers ####

ACS_EssentialWrkr_Commute <- get_acs(geography = "zcta", #pull down the relevant categories 
                                     variables = c(ag_car1_commute = "B08126_017",
                                                   ag_pubtrans_commute = "B08126_047",
                                                   construct_car1_commute ="B08126_018",
                                                   construct_pubtrans_commute = "B08126_048",
                                                   wholesale_car1_commute = "B08126_020",
                                                   wholesale_pubtrans_commute = "B08126_050",
                                                   transpo_car1_commute = "B08126_022",
                                                   transpo_pubtrans_commute = "B08126_052",
                                                   ed_hlthcare_car1_commute = "B08126_026",
                                                   ed_hlthcare_pubtrans_commute = "B08126_056"),
                                     year = 2018, 
                                     output = "wide",
                                     survey = "acs5")


ACS_EssentialWrkr_Commute1 <- ACS_EssentialWrkr_Commute %>% #clean data and aggregate 
  dplyr::select(-ends_with("M"), -NAME) %>%
  filter(GEOID %in% ZCTAs_in_NYC) %>%
  mutate_at(vars(starts_with("ed_hlthcare")), ~round(./2), 0) %>% #maintain same proportions as estimated nonquarintined jobs
  mutate_at(vars(starts_with("construct")), ~round(./4), 0) %>%
  mutate(essentialworker_drove = rowSums(dplyr::select(., contains("car1_commute"))), 
         essentialworker_pubtrans = rowSums(dplyr::select(., contains("pubtrans")))) %>%
  rename(zcta = GEOID) %>%
  dplyr::select(zcta, essentialworker_drove, essentialworker_pubtrans)

#### Identify the number of supermarkets/grocery stores per area ####
non_supermarket_strings <- c("DELI|TOBACCO|GAS|CANDY|7 ELEVEN|7-ELEVEN|LIQUOR|ALCOHOL|BAKERY|CHOCOLATE|DUANE READE|WALGREENS|CVS|RITE AID|RAVIOLI|WINERY|WINE|BEER|CAFE|COFFEE")

food_retail1 <- food_retail %>% #an estimate of how many supermarkets are in a ZCTA
  filter(County %in% NYC_boro_county_match$County) %>% #get rid of those that have the non_supermarket_strings words in their store & legal names
  filter(str_detect(`Establishment Type`, "J") & str_detect(`Establishment Type`, "A") & str_detect(`Establishment Type`, "C") &
           !str_detect(`Establishment Type`, "H")) %>%
  filter(!str_detect(`Entity Name`, non_supermarket_strings) & !str_detect(`DBA Name`, non_supermarket_strings)) %>%
  filter(`Square Footage`>=4500) %>%
  mutate(zcta = as.character(str_extract(Location, "[:digit:]{5}"))) %>%
  group_by(zcta) %>%
  summarise(grocers = n_distinct(`License Number`))

### Where are subway stations located? ###

#SubwayStation_shp <- as_tibble(turnstile()$stations) %>%
#  st_as_sf(., coords = c("lon", "lat"), crs = 4269) %>%
#  st_transform(., crs = 2263) %>%
#  filter(!str_detect(ca, "PTH")) #removing New Jersey PATH stations

### Calculate the residential area per ZCTA ###

Pluto_ResOnly <- Pluto %>%
  filter(landuse>="01" & landuse<="04") %>%
  mutate(base_bbl = as.character(bbl)) %>%
  dplyr::select(-bbl)

ResBBLs <- as.character(Pluto_ResOnly$base_bbl)

Res_Bldg_Footprints <- Bldg_Footprints %>%
  st_set_geometry(., NULL) %>%
  mutate(base_bbl = as.character(base_bbl)) %>%
  filter(base_bbl %in% ResBBLs &
           feat_code == "2100") %>%
  mutate(bldg_volume = shape_area * heightroof) %>%
  left_join(., Pluto_ResOnly, by = "base_bbl") %>%
  mutate(bldg_volume = if_else(is.na(bldg_volume), shape_area*numfloors*10, bldg_volume),
         res_volume = (bldg_volume/unitstotal)*unitsres, 
         zcta = as.character(zipcode)) %>%
  group_by(zcta) %>%
  summarise(total_res_volume_zcta = sum(res_volume, na.rm = TRUE)) 


#### COVID Tests  ####
MODZCTA_NYC_shp1 <- MODZCTA_NYC_shp %>%
  dplyr::select(modzcta, geometry) %>%
  rename("zcta" = "modzcta")


Oct1_tests <- ZCTA_test_series %>%
  mutate(modzcta = as.character(modzcta)) %>%
  rename(zcta = 'modzcta')
  rename(total_tests = 'Total') %>%
  dplyr::select(modzcta, Positive, total_tests, modzcta_cum_perc_pos)


ZCTA_by_boro1 <- ZCTA_by_boro %>%
  mutate(Borough = as.character(Borough))


# Supplemental Figure 1 (A - Tests)
sfig1a <- MODZCTA_NYC_shp1 %>%
  left_join(., Oct1_tests, by = "zcta") %>%
  left_join(., ACS_Data2, by = "zcta") %>%
  filter(zcta != "99999") %>%
  mutate(pos_per_100000 = (Positive/total_pop1)*100000) %>%
  ggplot() +
  geom_sf(data = NYC_basemap_shp)+
  geom_sf(aes(fill = pos_per_100000), lwd = 0.2)+
  labs(fill = "Positives per 100,000") +
  ggtitle("Cumulative Positive COVID tests by zip code (May 7, 2020)") +
  scale_fill_gradientn(colours=brewer_pal("BuPu", type = "seq")(7)) + 
  theme_bw() +
  theme_bw(base_size = 6) + 
  theme(legend.title = element_text(face = "bold", size = 6), 
        panel.background = element_rect(fill = "#dedede"), 
        legend.background = element_rect(fill = "transparent"),
        legend.position = c(0.15, 0.80),
        legend.text = element_text(size = 6),
        plot.margin = unit(c(4,0,4,0), "pt"),
        legend.key.size = unit(1.1, "lines"))
sfig1a 

#### Create data frames of all above information ####


ZCTA_ACS_COVID_shp <- MODZCTA_NYC_shp1 %>%
  st_transform(., crs = 2263) %>%
  dplyr::select(zcta, geometry) %>%
  left_join(., ACS_Data2, by = "zcta") %>%
  left_join(., Oct1_tests, by = "zcta") %>%
  left_join(., Res_Bldg_Footprints, by = "zcta") %>%
  left_join(., ACS_EssentialWrkr_Commute1, by = "zcta") %>%
  left_join(., food_retail1, by = "zcta") %>%
  mutate(pop_density = as.numeric(total_pop1/st_area(geometry)),
         avg_hhold_size = round((total_pop1/total_hholds1), 2),
         pos_per_100000 = (Positive/total_pop1)*100000,
         testing_ratio = (Total/total_pop1),
         res_vol_zctadensity = as.numeric(total_res_volume_zcta/st_area(geometry)), 
         res_vol_popdensity = as.numeric(total_pop1/total_res_volume_zcta),
         pubtrans_ferrysubway_commute = pubtrans_subway_commute + pubtrans_ferry_commute,
         grocers = replace_na(grocers, 0),
         grocers_per_1000 = (grocers/total_pop1)*1000,
         pos_per_100000 = round(pos_per_100000, 0),
         valid_var = "0",
         didnot_workhome_commute = 1/workhome_commute,
         one_over_grocers_per_1000 = if_else(is.infinite(1/grocers_per_1000), 0, 1/grocers_per_1000),
         one_over_medincome = 1/medincome) %>%
  dplyr::select(-pubtrans_subway_commute, -pubtrans_ferry_commute) %>%
  left_join(., ZCTA_by_boro1, by = "zcta") %>%
  mutate_at(vars(starts_with("essentialworker_")), ~round((./over16total_industry1)*100, 2)) %>%
  filter(zcta != "99999") #remove na


ZCTA_ACS_COVID <- ZCTA_ACS_COVID_shp %>%
  st_set_geometry(., NULL) #remove geometry

ZCTA_ACS_COVID_shp['geometry']

write.csv(ZCTA_ACS_COVID,"/home/claudia/Scrivania/Tesi NYC pandemic /COVID_19_data_Rstudio/Complete_data_per_zipcodes.csv", row.names = FALSE)

st_as_sf(ZCTA_ACS_COVID_shp)
class(ZCTA_ACS_COVID_shp)

st_write(ZCTA_ACS_COVID_shp['zcta'], "shapes", append = TRUE, driver = "ESRI Shapefile")

