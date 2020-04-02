# Fishery-dependent-SDMs
Code for project looking at impact of using fishery-dependent data in SDM and forecasting distributions under future climate change

OPERATING MODELS:

SimulatedWorld_ROMS_TrophicInteraction_SamplingBias.R: this function uses temperature and chl-a (from ROMS) to build suitability for two species (predator and its prey), and generate distribution and abundance for one species adding a sampling bias to mimic data collected from a fishery (i.e. fishery-dependent observations) 

SimulatedWorld_ROMS_SamplingBias.R: this function uses temperature (from ROMS) to build suitability for a species and generate distribution and abundance that species adding a sampling bias to mimic data collected from a fishery (i.e. fishery-dependent observations) 

ESTIMATION MODELS:

ModelComparison_TrophicInteractions_SamplingBias.R: this code uses the SimulateWorld_ROMS_TrophicInteraction_SB function above to generate data, then builds an example GAM and BRT model, and makes predictions into the future 2020-2100. 

ModelComparison_SamplingBias.R: this code uses the SimulateWorld_ROMS_SB function above to generate data, then builds an example GAM and BRT model, and makes predictions into the future 2020-2100, using just temp as predictor.
