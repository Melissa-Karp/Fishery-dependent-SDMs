# Fishery-dependent-SDMs
Code for project looking at impact of using fishery-dependent data in SDM and forecasting distributions under future climate change

USE THESE FILES: 
1. DistancetoPorts.R: This code calculates the distance from every cell in the ROMS extent to 5 different fishing ports along the US West Coasts of CA, OR, and WA. 

2. Operating Model:
SimulatedWorld_ROMS_FishDep_UnequalCoverage.R: Simulates the sampling design for the different fishery location choose biases with the detection probability = 1 with no correction for env suitability of the location. That is to say that is the species is present at a location it will always be detected. This updated OM increases the strength of the fishermen preference for high suitable habitat for target species, by including a new bias siutaiton called opt_sampled. This also changes to only using dp4 and dp5 (ports in WA and OR) for the distance bias situation to simulate a situation where fishermen can only land fish in WA or OR (inc strenght/impact of the distance effect on sampling locations). 

3. Estimation Model:
ModelComparison_FishSuitability_v10_6_2020.R: this code uses the SimulatedWorld_ROMS_FishDep_NoDetProbCorrection.R function above to generate data, then builds an example GAM, makes predictions into the future 2011-2100, and plots results. 


OLD CODE: 

OPERATING MODELS:

For all operating models we are simulating a "Pelagic mobile predator" like species. For all operating models the caluclation of the species suitability functions are the same and are a function of mean spring SST, mld, and a prey species (Sp A), whose suitability is a function of zooplankton and SST. 

SimulatedWorld_ROMS_FishDep_NoDetProbCorrection.R: Simulates the sampling design for the different fishery location choose biases with the detection probability = 1 with no correction for env suitability of the location. That is to say that is the species is present at a location it will always be detected.

SimulatedWorld_ROMS_FishDep_DetProbCorrected.R: Simulates the sampling design for the different fishery location choose biases with the detection probability = 1 but corrected to be a function of habitat suitability. That is to say that is the species is present at a location the chances of it being detected (or caught be fishermen) is dep on the habitat suitability at that location. Higher habitat suit means greater abundance and therefore higher chances of species being caught. 

SimulatedWorld_ROMS_SB_FunctionFishingSuit.R: This function builds a Fishing Suitability Raster layer, by treating "Fishing" as a predator species whose 'suitable' sites/habitat is a function of (1)suitability of habitat for its prey (i.e. target species, spB), (2) Distance to port - proxy for cost, and (3) low probability of encountering bycatch species (so likes areas of low suitability for spC). This function builds 4 different Fish suitability layers - (i). Just based on habitat suit of Sp B, (ii) based on habitat suit of sp B and distance to port, (iii) hab suit sp B and low bycatch risk, and (iii) hab suit spB, distance to port, and bycatch risk. 

SimulatedWorld_ROMS_FishDepFun.R: Simulates the 5 functions below all in one function by allowing user to indicate sample size using nsamples = x. This function simulates fishing biases by building a Random Utility function. 

SimulatedWorld_ROMS_RS400.R: This function does not include any sampling bias, but samples 400 random points each year of the simulation.

SimulatedWorld_ROMS_RS100.r: This function does not include any sampling bias, but samples 100 random points each year of the simulation.

SimulatedWolrd_ROMS_PrefSamp400.R: This function includes a sampling bias where it samples 400 points each year with a 20x higher sampling intensity within the area of highest suitability of the target species (Sp B). 

SimulatedWolrd_ROMS_PrefSamp100.R: This function includes a sampling bias where it samples 100 points each year with a 20x higher sampling intensity within the area of highest suitability of the target species (Sp B). 

SimulatedWolrd_ROMS_RUM400.R: This function includes a sampling bias where it samples 400 points each year based on a Randum Utility model where "fishing location choice" is a function of "utility" which is calcuatled based on the expected catch in the previous year (based on suitability of target species from previous year of simulation) and cost, which is a function of distance from port and price of fuel. 

ESTIMATION MODELS:

ModelComparison_TrophicInteractions_SamplingBias.R: this code uses the SimulateWorld_ROMS_TrophicInteraction_SB function above to generate data, then builds an example GAM and BRT model, and makes predictions into the future 2020-2100. 

