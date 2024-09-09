- - - -
# HAL
**H**ydrodynamics-based **A**lgorithm for **L**agrangian simulations
- - - -

## Introduction

HAL is a Lagrangian model for calculating the transport of zooplankton such as copepods, the egg and larval phases of marine fish, and the early life stages of crustaceans such as krill and lobsters.

## Model Overview

HAL is written in FORTRAN and has evolved over many years to now offer a range of options for different species-specific biology and behaviours:

* Growth and development, based on temperature-dependent relationships, with default values for plaice
* Egg buoyancy, specified as either a vertical velocity or as a particle density, from which a vertical velocity is calculated based on a comparison of the particle and local water densities
* Diel vertical migration, specified as either nocturnal or diurnal (reverse) upward migration
* Directed horizontal swimming, relative to currents or bathymetric gradients
* Transport with sea ice
* Maximum and minimum depths of vertical distribution
* Seasonal diapause behaviour with separate vertical distribution limits for the surface and diapause phases

In addition, the model has been adapted for the modelling of the fate of oil spills, and includes options for buoyant oil droplets and beaching. If no behaviours are chosen, the model simulates the passive transport of any neutrally-buoyant particle.

## Installation instructions

HAL is designed to run in a Linux environment as a serial job on a dedicated processor. The following instructions assume a basic understanding of Linux commands.
* Create a home directory 'HAL'
* Copy the directory 'Source' to your home directory
* Load netCDF and HDF5 libraries and Intel compilers; the current version of HAL is tested with the following:
  * netcdf/intel/4.4.1.1
  * hdf5/intel/1.8.19
  * intel/oneapi/2021.1.0.2659
* Navigate to the Source directory and compile the code with the command _make_
* Copy the resultant executable (ibm) to the home directory

Detailed instructions on how to set up and run the model for different applications are provided in the documentation.

## Peer-reviewed articles using HAL

Young, E.F., Thorpe, S.E., Renner, A.H.H., Murphy, E.J. (2024). ‘_Environmental and behavioural drivers of Antarctic krill distribution at the South Orkney Islands: a regional perspective_’. Journal of Marine Systems, 241, https://doi.org/10.1016/j.jmarsys.2023.103920.

Trathan, P.N., Warwick-Evans, V., Young, E.F., Friedlaender, A., Kim, J.H., Kokubun, N. (2022). ‘_The ecosystem approach to management of the Antarctic krill fishery - the ‘devils are in the detail’ at small spatial and temporal scales_’. Journal of Marine Systems, 225, https://doi.org/10.1016/j.jmarsys.2021.103598.

Campanella, F., Collins, M.A., Young, E.F., Laptikhovsky, V., Whomersley, P., van der Kooij, J. (2021). ‘_First insight of Meso- and Bentho-Pelagic Fish Dynamics Around Remote Seamounts in the South Atlantic Ocean_’. Frontiers in Marine Science, 8, doi: 10.3389/fmars.2021.663278.

Silva, C.N.S., Young, E.F., Murphy, N.P., Bell, J.J., Green, B.S., Morley, S.A., Duhamel, G., Cockcroft, A.C., Strugnell, J.M. (2021). ‘_Climatic change drives dynamic source-sink relationships in marine species with high dispersal potential_’. Ecology and Evolution, doi: 10.1002/ece3.7204.

Young, E.F., Tysklind, N., Meredith, M.P., de Bruyn, M., Belchier, M., Murphy, E.J., Carvalho, G.R. (2018). ‘_Stepping stones to isolation: impacts of a changing climate on the connectivity of fragmented fish populations_’, Evolutionary Applications, doi: 10.1111/eva.12613.

Young, E.F., Belchier, M., Hauser, L., Horsburgh, G.J., Meredith, M.P., Murphy, E.J., Pascoal, S., Rock, J., Tysklind, N., Carvalho, G.R. (2015). ‘_Oceanography and life history predict contrasting genetic population structure in two Antarctic fish species_’. Evolutionary Applications, doi: 10.1111/eva.12259.

Young, E.F., Thorpe, S.E., Banglawala, N., Murphy, E.J. (2014). '_Variability in transport pathways on and around the South Georgia shelf, Southern Ocean: implications for recruitment and retention_', Journal of Geophysical Research: Oceans, 119, 1-12, doi: 10.1002/2013JC009348.

Young, E.F., Rock, J., Meredith, M.P., Belchier, M., Murphy, E.J., Carvalho, G.R. (2012). ‘_Physical and behavioural influences on larval fish retention: contrasting patterns in two Antarctic fishes_’. Marine Ecology Progress Series 465, 201-215.

Fox, C.J., McCloughrie, P., Young, E.F., Nash, R.D.M. (2006). ‘_The importance of individual behaviour for successful settlement of juvenile plaice (Pleurenectes platessa L.): a modelling and field study in the eastern Irish Sea_’. Fisheries Oceanography, 15 (4), 301-313.







