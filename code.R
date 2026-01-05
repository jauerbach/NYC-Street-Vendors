# This code recreates tables, figures, and other results from the manuscript
# "Estimating the Number of Street Vendors in New York City:
#  Ratio Estimation with Point Process Data." 

# Manuscript: https://doi.org/10.1093/jssam/smaf051

# Contact: Jonathan Auerbach jauerba@gmu.edu

#########
# SETUP #
#########

library("tidyverse")
library("sf")
# install.packages("cmdstanr", repos = c('https://stan-dev.r-universe.dev', getOption("repos")))
library("cmdstanr")
library("posterior")

###########
# TABLE 1 #
###########

food_table <- read_csv("data/food_table.csv")
merch_table <- read_csv("data/merch_table.csv")

food_table %>% 
  transmute(Borough,
            Neighborhood, 
            `Food Respondents` = no + yes,
            `Food Population` = 5100 * `Food Respondents` / 349,
            `Food Error` = 2 * 5100 * no / 349 * 
              sqrt(1/no + 1/349 - 1/5100 + (yes * (349 - yes)) / (349 * no^2))) %>%
  left_join(
    merch_table %>% 
      transmute(Borough,
                Neighborhood,
                `Merch Respondents` = no + yes,
                `Merch Population` = 853 * `Merch Respondents` / 308,
                `Merch Error` = 2 * 853 * no /308 * 
                sqrt(1/no + 1/308 - 1/853 + (yes * (308 - yes)) / (308 * no^2)))) %>%
  mutate_if(is.numeric, list(~replace_na(., 0))) %>%
  transmute(Borough,
            Neighborhood,
            Respondents = `Food Respondents` + `Merch Respondents`,
            Population = `Food Population` + `Merch Population`,
            Error = `Food Error` + `Merch Error`) %>%
  mutate(order = ifelse(str_sub(Neighborhood, 1, 5) == "Other", 1, 
                        ifelse(Neighborhood == "Total", 2, 0))) %>%
  arrange(Borough, order, desc(Population)) %>%
  dplyr::select(-Borough, -order)


############
# FIGURE 1 #
############

map11 <- st_read("data/map11.shp")
map11 %>%
  mutate(n = ifelse(is.na(n), 0, n)) %>%
  ggplot() +
  aes(fill = log(n, base = 10)) +
  geom_sf() +
  theme(legend.position = "bottom") +
  labs(fill = "vending location of respondents\nsource: Street Vendor Project") +
  scale_fill_gradient(breaks = c(.5, 1, 1.5, 2, 2.5),
                      labels = c(1, 10, 30, 100, 300),
                      low = "white",
                      high = "#18B686",
                      na.value = "white")

map12 <- st_read("data/map12.shp")
map12 %>%
  ggplot() +
  aes(fill = log(n, base = 10)) +
  geom_sf() +
  theme(legend.position = "bottom") +
  labs(fill = "vending location of individuals summoned\nto court for first vending-related violation\nsource: OATH 2021") +
  scale_fill_gradient(#breaks = c(.5, 1, 1.5, 2, 2.5),
    #labels = c(1, 10, 30, 100, 300),
    breaks = c(.5, 1, 1.5, 2, 2.5),
    labels = c(1, 10, 30, 100, 300),
    low = "white",
    high = "#18B686",
    na.value = "white")

map21 <- st_read("data/map21.shp")
map21 %>%
  mutate(n = ifelse(is.na(n), 0, n)) %>%
  ggplot() +
  aes(fill = log(n, base = 10)) +
  geom_sf() +
  theme(legend.position = "bottom") +
  labs(fill = "residency of respondents\nsource: Street Vendor Project") +
  scale_fill_gradient(breaks = c(.5, 1, 1.5, 2, 2.5),
                      labels = c(1, 10, 30, 100, 300),
                      low = "white",
                      high = "#18B686",
                      na.value = "white")

map22 <- st_read("data/map22.shp")
map22 %>%
  ggplot() +
  aes(fill = log(n, base = 10)) +
  geom_sf() +
  theme(legend.position = "bottom") +
  labs(fill = "residency of door-to-door, news\nand street vendors, and related\nsource: ACS 2017-2021") +
  scale_fill_gradient(breaks = c(.5, 1, 1.5, 2, 2.5),
                      labels = c(1, 10, 30, 100, 300),
                      low = "white",
                      high = "#18B686",
                      na.value = "white")


####################
# FORDHAM ESTIMATE #
####################

# food vendors on "FORDHAM ROAD"
tibble(no = 9, yes = 1) %>%
  mutate(`Fordham MFV Population` = 5100 * (yes + no) / 349,
         `Food Error` = 5100 * no / 349 * 
           sqrt(1/no + (1/349 - 1/5100) + (yes * (349 - yes) / (349 * no^2))))

# merchandise vendors on "FORDHAM ROAD"
tibble(no = 7, yes = 0) %>%
mutate(`Fordham MFV Population` = 853 * (yes + no) / 303,
       `Food Error` = 853 * no / 303 * 
         sqrt(1/no + (1/308 - 1/853) + (yes * (303 - yes) / (303 * no^2))))


################
# ACS ESTIMATE #
################

# Download Census Bureau American Community Survey (ACS) 
## Public Use Microdata Sample (PUMS) file "csv_pus" from
## https://www2.census.gov/programs-surveys/acs/data/pums/2022/5-Year/ 

acs_usa <- read_csv("data/csv_pus/psam_pusa.csv") %>%
  filter(SOCP == 419091,
         POWPUMA %in% str_c(str_sub(paste0(0, map22$puma), 1, 3), "00"),
         POWSP %in% "036") %>%
  dplyr::select(PWGTP, PWGTP1:PWGTP80)

acs_usb <- read_csv("data/csv_pus/psam_pusb.csv") %>%
  filter(SOCP == 419091,
         POWPUMA %in% str_c(str_sub(paste0(0, map22$puma), 1, 3), "00"),
         POWSP %in% "036") %>%
  dplyr::select(PWGTP, PWGTP1:PWGTP80)

acs_usc <- read_csv("data/csv_pus/psam_pusc.csv") %>%
  filter(SOCP == 419091,
         POWPUMA %in% str_c(str_sub(paste0(0, map22$puma), 1, 3), "00"),
         POWSP %in% "036") %>%
  dplyr::select(PWGTP, PWGTP1:PWGTP80)

acs_usd <- read_csv("data/csv_pus/psam_pusd.csv") %>%
  filter(SOCP == 419091,
         POWPUMA %in% str_c(str_sub(paste0(0, map22$puma), 1, 3), "00"),
         POWSP %in% "036") %>%
  dplyr::select(PWGTP, PWGTP1:PWGTP80)

estimate_acs <- 
  sum(acs_usa$PWGTP) + 
  sum(acs_usb$PWGTP) + 
  sum(acs_usc$PWGTP) + 
  sum(acs_usd$PWGTP)

se_acs <- 
  acs_usa %>% 
  bind_rows(acs_usb) %>% 
  bind_rows(acs_usc) %>% 
  bind_rows(acs_usd) %>%
  dplyr::select(PWGTP, PWGTP1:PWGTP80) %>%
  mutate(across(PWGTP1:PWGTP80, function(x) x - PWGTP)) %>%
  dplyr::select(-PWGTP) %>%
  colSums() %>%
  (function(x) sqrt(4 / 80 * sum(x^2)))


#############
# WEIGHTING #
#############

weights_food_oath <- read_csv("data/weights_food_oath.csv")
weights_merch_oath <- read_csv("data/weights_merch_oath.csv")
weights_food_acs <- read_csv("data/weights_food_acs.csv")
weights_merch_acs <- read_csv("data/weights_merch_acs.csv")

weights_food_oath %>%
  #weights_food_acs %>%
  group_by(permit) %>%
  summarize(count = sum(count),
          weighted_count = sum(weighted_count),  
          weighted_count_sq = sum(weighted_count_sq)) %>%
  pivot_wider(names_from = permit, values_from = c(count, weighted_count, weighted_count_sq)) %>%
  transmute(weighted_count = weighted_count_no / weighted_count_yes * 5100 + 5100,
            se = weighted_count * sqrt(weighted_count_sq_no / weighted_count_no^2 + 
                                         weighted_count_sq_yes / weighted_count_yes^2  - 1/5100))

weights_merch_oath %>%
  #weights_merch_acs %>%
  group_by(permit) %>%
  summarize(count = n(),
          weighted_count = sum(weighted_count),
          weighted_count_sq = sum(weighted_count_sq)) %>% 
  pivot_wider(names_from = permit, values_from = c(count, weighted_count, weighted_count_sq)) %>%
  transmute(weighted_count = weighted_count_no / weighted_count_yes * 853 + 853,
            se = weighted_count * sqrt(weighted_count_sq_no / weighted_count_no^2 + weighted_count_sq_yes / weighted_count_yes^2 - 1/853))


######################
# HIERARCHICAL MODEL #
######################

food_stan <- read_csv("data/food_stan.csv")
merch_stan <- read_csv("data/merch_stan.csv")

stan_data_food <- list(
  K  = length(unique(food_stan$zip)),
  n0 = as.integer(food_stan$count[food_stan$permit == "no"]),
  n1 = as.integer(food_stan$count[food_stan$permit == "yes"]),
  N1  = as.integer(5100),
  rho = rep(10, length(unique(food_stan$zip)))
)

stan_data_merch <- list(
  K  = length(unique(merch_stan$zip)),
  n0 = as.integer(merch_stan$count[merch_stan$permit == "no"]),
  n1 = as.integer(merch_stan$count[merch_stan$permit == "yes"]),
  N1  = as.integer(853),
  rho = rep(10, length(unique(merch_stan$zip)))
)

library("cmdstanr")
library("posterior")

hierarchical_poisson <- cmdstan_model("hierarchical_poisson.stan")
#hierarchical_negbin <- cmdstan_model("hierarchical_negbin.stan")

fit_food <- hierarchical_poisson$sample(  
  #hierarchical_negbin$sample(
  data = stan_data_food,
  seed = 1,
  chains = 4,
  parallel_chains = 4,
  iter_warmup = 2000,
  iter_sampling = 2000,
  adapt_delta = .99,
  max_treedepth = 15)


fit_merch <- hierarchical_poisson$sample(  
  #hierarchical_negbin$sample(
  data = stan_data_merch,
  seed = 1,
  chains = 4,
  parallel_chains = 4,
  iter_warmup = 2000,
  iter_sampling = 2000,
  adapt_delta = .99,
  max_treedepth = 15)

lambda0_food <- as_draws_matrix(fit_food$draws("lambda0"))
lambda0_merch <- as_draws_matrix(fit_merch$draws("lambda0"))

# Posterior Quantiles
quantile(rowSums(lambda0_food) + 5100, c(.025, .5, .975)) 
quantile(rowSums(lambda0_merch) + 853 + 1000, c(.025, .5, .975))
quantile(rowSums(lambda0_food) + rowSums(lambda0_merch) + 5100 + 853 + 1000, c(.025, .5, .975))

# MAP
sum(lambda0_food[which.max(as_draws_matrix(fit_food$draws("lp__"))), ]) + 5100 + 
  sum(lambda0_merch[which.max(as_draws_matrix(fit_merch$draws("lp__"))), ]) + 853 + 1000

# Posterior Mean
mean(rowSums(lambda0_food)) + 5100 + 
  mean(rowSums(lambda0_merch)) + 853 + 1000
