library(tidyverse)
library(ConvergenceClubs)
library(sandwich)
library(lmtest)
library(broom)

gdim = read_csv("https://www.dropbox.com/scl/fi/rfhjr8vsozbm3tf9b0hl1/GDIM_2023_03.csv?rlkey=e7xkxskhd2tqrh13lp13th9xr&st=2wdvfvf0&dl=1")

# Restrict to countries with observations for 1940-1980 cohorts
keep_1940 = gdim |> group_by(country) |> count() |> filter(n == 60) |> pull(country)
gdim_1940 = gdim |> filter(country %in% keep_1940)

# Restrict to countries with observations for 1950-1980 cohorts
keep_1950 = gdim |> group_by(country) |> count() |> filter(n >= 48) |> pull(country)
gdim_1950 = gdim |> filter(country %in% keep_1950)

x_1940 = gdim_1940 |>
  filter(parent == "max", child == "all") |>
  mutate(ige = 1 - BETA) |>
  select(country, cohort, ige) |>
  pivot_wider(names_from = "cohort", values_from = "ige") |>
  arrange(country) |>
  as.data.frame()

H = computeH(x_1940[,-1], quantity = "H")
round(estimateMod(H, time_trim=1/3, HACmethod = "FQSB"), 3)

clubs = findClubs(x_1940, dataCols = 2:6, unit_names = 1, refCol = 6,
                  time_trim = 0.2, cstar = 0, HACmethod = 'FQSB')  
