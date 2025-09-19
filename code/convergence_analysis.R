library(tidyverse)
library(ConvergenceClubs)
library(sandwich)
library(lmtest)
library(broom)

gdim = read_csv("https://www.dropbox.com/scl/fi/rfhjr8vsozbm3tf9b0hl1/GDIM_2023_03.csv?rlkey=e7xkxskhd2tqrh13lp13th9xr&st=2wdvfvf0&dl=1")

# Helper to compute indicator column from measure choice
.add_indicator = function(df, measure){
  stopifnot(measure %in% c("1-beta", "1-cor", "MU050", "BHQ4"))
  if (measure == "1-beta") {
    df |> mutate(val = 1 - .data[["BETA"]])
    } else if (measure == "1-cor") {
    df |> mutate(val = 1 - .data[["COR"]])
    } else if (measure == "MU050") {
    df |> mutate(val == .data[["MU050_tiebreak"]])
    } else if (measure == "BHQ4") {
    df |> mutate(val == .data[["BHQ4_tiebreak"]])
    }
  }

# Function to build balanced panel with relevant indicators
make_panel = function(data, parent = "max", child = "all", measure = "1-beta"){
  df = gdim |> filter(.data$parent = "parent", .data$child = "child") |>
  .add_indicator("measure") |>
  select(country, cohort, val)

  wide = df |>
  pivot_wider(names_from = "cohort", values_from = "val") |>
  arrange(country) |>
  as.data.frame()

  wide
  }

# Function to run Phillips-Sul tests
run_ps = function(panel_df, dataCols, time_trime = 0.2, HACmethod = "FQSB", cstar = 0, refCol = NULL) {
  stopifnot(is.data.frame(panel_df))
  if(is.null(refCol)) {
    refCol = max(dataCols)
    }

  H = computeH(panel_df[, dataCols], quantity = "H")
  global = estimateMod(H, time_trim = time_trim, HACmethod = HACmethod)

  clubs = findClubs(panel_df, dataCols = dataCols, unit_names = 1, refCol = refCol,
                    time_trim = time_trim, cstar = cstar, HACmethod = HACmethod)

  list(H = H, global = global, clubs = clubs)
  }

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
