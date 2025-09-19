library(tidyverse)
library(ConvergenceClubs)
library(sandwich)
library(lmtest)
library(broom)

gdim = read_csv("https://www.dropbox.com/scl/fi/rfhjr8vsozbm3tf9b0hl1/GDIM_2023_03.csv?rlkey=e7xkxskhd2tqrh13lp13th9xr&st=2wdvfvf0&dl=1")

# Helper to compute indicator column from measure choice
.add_indicator = function(df, measure){
  key = tolower(gsub("[^a-z0-9]", "", measure))

  if (key %in% c("1beta","1cor","mu050","bhq4",
                 "mu050randomtiebreak","bhq4randomtiebreak")) {
    if (key == "1beta") {
      df$val = 1 - df$BETA
    } else if (key == "1cor") {
      df$val = 1 - df$COR
    } else if (key %in% c("mu050","mu050randomtiebreak")) {
      df$val = df$MU050_randomtiebreak
    } else if (key %in% c("bhq4","bhq4randomtiebreak")) {
      df$val = df$BHQ4_randomtiebreak
    }
    return(df[, c("country","cohort","val")])
  } else {
    stop("`measure` not recognized. Try one of: '1-beta','1-cor','MU050','BHQ4'.")
  }
}


# Function to build balanced panel with relevant indicators
make_panel = function(data, parent = "max", child = "all", measure = "1-beta"){
  df = gdim |> filter(.data$parent == "parent", .data$child == "child") |>
  .add_indicator("measure") |>
  select(country, cohort, val)

  wide = df |>
  pivot_wider(names_from = "cohort", values_from = "val") |>
  arrange(country) |>
  as.data.frame()

  wide
  }

# Function to run Phillips-Sul tests
run_ps = function(panel_df, dataCols, time_trim = 0.2, HACmethod = "FQSB", cstar = 0, refCol = NULL) {
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
