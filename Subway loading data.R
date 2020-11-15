###Script for subway usage data

library(tidyverse)
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

##subway load data
Subway_ridership_by_UHF <- relative.subway.usage(2020L, "nhood")

### Where are subway stations located? ###

SubwayStation_shp <- as_tibble(turnstile()$stations) %>%
  st_as_sf(., coords = c("lon", "lat"), crs = 4269) %>%
  st_transform(., crs = 2263) %>%
  filter(!str_detect(ca, "PTH")) #removing New Jersey PATH stations

