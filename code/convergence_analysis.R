library(tidyverse)
library(ConvergenceClubs)
library(sandwich)
library(lmtest)
library(broom)

gdim = read_csv("https://www.dropbox.com/scl/fi/rfhjr8vsozbm3tf9b0hl1/GDIM_2023_03.csv?rlkey=e7xkxskhd2tqrh13lp13th9xr&st=2wdvfvf0&dl=1")

# Helper to compute indicator column from measure choice
.add_indicator = function(df, measure){
 
  key = tolower(gsub("[^a-z0-9]", "", measure))
  
  if (key == "1beta") {
    df = df |> mutate(val = 1 - .data[["BETA"]])
  } else if (key == "1cor") {
    df = df |> mutate(val = 1 - .data[["COR"]])
  } else if (key == "mu050") {
    df = df |> mutate(val = .data[["MU050_randomtiebreak"]])
  } else if (key == "bhq4") {
    df = df |> mutate(val = .data[["BHQ4_randomtiebreak"]])
  } else if (key == "ahmp") {
    df = df |> mutate(val = .data[["CAT_ISCED0"]])
  } else if (key == "mix") {
    df = df |> mutate(val = .data[["MIX"]])
  } else {
    stop("`measure` not recognized. Use one of: '1-beta', '1-cor', 'mu050', 'bhq4', 'ahmp', 'mix'.")
  }
  
  df |> select(country, cohort, val)
}


# Function to build balanced panel with relevant indicators
make_panel = function(data, parent = "max", child = "all", measure = "1-beta"){
  
  df = data |> filter(.data$parent == !!parent, .data$child == !!child) |>
    .add_indicator(measure)
  
  wide = df |>
    pivot_wider(names_from = "cohort", values_from = "val") |>
    arrange(country) |>
    as.data.frame()
  
  wide
}

# Function to run Phillips-Sul tests
run_ps = function(data,
                  parent = "max",
                  child = "all",
                  measure = "1-beta",
                  dataCols,
                  time_trim = 0.3,
                  HACmethod = "FQSB",
                  cstar = 0,
                  refCol = NULL) {
  
  panel_df = make_panel(data, parent, child, measure)
  
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



