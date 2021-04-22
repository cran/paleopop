## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(poems)
library(paleopop)
parallel_cores <- 2
output_dir <- tempdir()

## ----siberia_raster, fig.align = "center", fig.width = 7, fig.height = 5------
library(raster)
raster::plot(paleopop::siberia_raster[[1]], main = "Siberia raster (LGM)",
             xlab = "Longitude (degrees)", ylab = "Latitude (degrees)",
             colNA = "blue")
raster::plot(paleopop::siberia_raster[[1001]], main = "Siberia raster (9000 BP)",
             xlab = "Longitude (degrees)", ylab = "Latitude (degrees)",
             colNA = "blue")

## ----paleoregion, fig.align = "center", fig.width = 7, fig.height = 5---------
region <- PaleoRegion$new(template_raster = paleopop::siberia_raster)
raster::plot(region$region_raster, main = "Bison region maximum extent (cell indices)",
             xlab = "Longitude (degrees)", ylab = "Latitude (degrees)",
             colNA = "blue")

## ----generate dispersal matrix------------------------------------------------
dispersal_gen <- DispersalGenerator$new(region = region,
                                        dispersal_max_distance = 500, # km
                                        distance_classes = seq(10, 500, 10), # divide into bins
                                        distance_scale = 1000, # sets units to km
                                        dispersal_friction = DispersalFriction$new(), # empty
                                        inputs = c("dispersal_p", "dispersal_b"),
                                        decimals = 3)
dispersal_gen$calculate_distance_data(use_longlat = TRUE) # pre-calculate matrix of distances between cells
test_dispersal <- dispersal_gen$generate(input_values =
                                                       list(dispersal_p = 0.5,
                                                            dispersal_b = 200))$dispersal_data
head(test_dispersal[[1]])

## ----hs_raster, fig.align = "center", fig.width = 7, fig.height = 5-----------
raster::plot(bison_hs_raster[[1]], main = "Bison habitat suitability (LGM)",
             xlab = "Longitude (degrees)", ylab = "Latitude (degrees)",
             colNA = "blue")

## ----carrying capacity generator, fig.align = "center", fig.width = 7, fig.height = 5----
# Convert the habitat suitability raster to matrix
bison_hs <- bison_hs_raster[region$region_indices]
# Create the generator
capacity_gen <- Generator$new(description = "Capacity generator",
                              example_hs = bison_hs, 
                              inputs = c("density_max"),
                              outputs = c("initial_abundance", "carrying_capacity"))
# indicate that initial abundance and carrying capacity are required
capacity_gen$add_generative_requirements(list(initial_abundance = "function",
                                              carrying_capacity = "function"))
# define custom generative functions
capacity_gen$add_function_template("carrying_capacity",
                                   function_def = function(params) {
                                     round(params$density_max*params$example_hs)
                                   },
                                   call_params = c("density_max", "example_hs"))
capacity_gen$add_function_template("initial_abundance",
                                   function_def = function(params) {
                                     params$carrying_capacity[,1]
                                   },
                                   call_params = c("carrying_capacity"))
# test run
test_capacity <- capacity_gen$generate(input_values = list(density_max = 2000))
raster::plot(region$raster_from_values(test_capacity$initial_abundance), main = "Initial abundance",
             xlab = "Longitude (degrees)", ylab = "Latitude (degrees)",
             colNA = "blue")

## ----paleopopmodel attributes-------------------------------------------------
model_template <- PaleoPopModel$new()
model_template$model_attributes # all possible attributes
model_template$required_attributes # required attributes to run simulations

## ----autocorrelation----------------------------------------------------------
# Distance-based environmental correlation (via a compacted Cholesky decomposition)
env_corr <- SpatialCorrelation$new(region = region, amplitude = 0.4, breadth = 500)
correlation <- env_corr$get_compact_decomposition(decimals = 2)

## ----human density------------------------------------------------------------
human_density <- array(rep(0.1), c(913, 1001)) # rows = populations and columns = timesteps

## ----model template-----------------------------------------------------------
model_template <- PaleoPopModel$new(
  simulation_function = "paleopop_simulator", # this is the default; made it explicit for the example
  region = region,
  # coordinates: you could use XY coordinates to define the region instead
  time_steps = 1001,
  years_per_step = 12, # generational length for bison
  populations = region$region_cells, # extracts number of population cells
  # initial_abundance: generated
  transition_rate = 1.0, # rate of transition between generations
  # standard_deviation: sampled
  compact_decomposition = correlation,
  # carrying_capacity: generated
  density_dependence = "logistic",
  # growth_rate_max: sampled
  harvest = TRUE,
  # harvest_max: sampled
  harvest_g = 0.4, # constant
  # harvest_z: sampled
  # harvest_max_n: sampled
  human_density = human_density,
  dispersal_target_k = 10, # minimum carrying capacity that attracts dispersers
  # dispersal_data: generated
  # abundance_threshold: sampled
  # initial_n: sampled
  occupancy_threshold = 1, # threshold: # of populations for the species to persist
  random_seed = 321,
  results_selection = c(
    # ema (expected minimum abundance), extirpation, occupancy, human_density,
    "abundance", "harvested"
    )
)

## ----latin hypercube sampling-------------------------------------------------
nsims <- 100 # adjust to run your own example if desired

lhs_generator <- LatinHypercubeSampler$new()
lhs_generator$set_uniform_parameter("standard_deviation", lower = 0, upper = sqrt(0.06))
lhs_generator$set_uniform_parameter("growth_rate_max", lower = log(1.31), upper = log(2.84))
lhs_generator$set_uniform_parameter("abundance_threshold", lower = 0, upper = 500, decimals = 0)
lhs_generator$set_uniform_parameter("harvest_max", lower = 0, upper = 0.35)
lhs_generator$set_uniform_parameter("harvest_z", lower = 1, upper = 2)
lhs_generator$set_uniform_parameter("density_max", lower = 500, upper = 3250) # alias for harvest_max_n
lhs_generator$set_uniform_parameter("dispersal_p", lower = 0.05, upper = 0.25) # for the dispersal generator
lhs_generator$set_uniform_parameter("dispersal_b", lower = 65, upper = 145) # for the dispersal generator

sample_data <- lhs_generator$generate_samples(number = nsims, random_seed = 123)
sample_data$sample <- c(1:nsims)
head(sample_data)

## ----simulate-----------------------------------------------------------------
sim_manager <- SimulationManager$new(sample_data = sample_data, 
                                     model_template = model_template,
                                     generators = list(capacity_gen,
                                                       dispersal_gen),
                                     parallel_cores = parallel_cores,
                                     results_dir = output_dir) 

sim_manager$results_filename_attributes <- c("sample", "results")

run_output <- sim_manager$run()
run_output$summary

## ----results------------------------------------------------------------------
results1 <- PaleoPopResults$new(
  results = readRDS(paste0(output_dir, "/sample_1_results.RData")),
  region = region, time_steps = 1001
  )
head(results1$extirpation)

## ----results manager----------------------------------------------------------
results_model <- PaleoPopResults$new(region = region, time_steps = 1001, trend_interval = 1:10)
metrics_manager <- ResultsManager$new(simulation_manager = sim_manager,
                                      simulation_results = results_model,
                                      generators = NULL) # don't need generators
metrics_manager$summary_metrics <- c("abundance_trend", "extinction_time")
metrics_manager$summary_functions <- list()

metrics_manager$summary_functions$extinction_time <- function(simulation_results) {
  extinction_time <- -12*(1001 - simulation_results$all$extirpation)-9001 # converts timestep to years BP
  if (is.na(extinction_time)) {
    extinction_time <- -9001 # if NA, then it survived to the end of the simulation
  }
  return(extinction_time)
}

metrics_manager$summary_functions$abundance_trend <- function(simulation_results) {
  abundance_trend <- simulation_results$all$abundance_trend
  return(abundance_trend) # this will use the trend interval specified above
}

gen_log <- metrics_manager$generate(results_dir = output_dir)

## ----histograms, fig.align = "center", fig.width = 7, fig.height = 5----------
hist(metrics_manager$summary_metric_data$abundance_trend, 
     main = "Histogram of abundance trend", xlab = "Abundance trend")
hist(metrics_manager$summary_metric_data$extinction_time, 
     main = "Histogram of extinction time", xlab = "Extinction time (years BP)")

