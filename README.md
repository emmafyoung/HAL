- - - -
# HAL
**H**ydrodynamics-based **A**lgorithm for **L**agrangian simulations
- - - -

## Introduction

HAL is a Lagrangian model for calculating the trajectories of the early life stages of marine biota, including the egg and larval phases of marine fish and the early life stages of crustacea such as krill and lobsters.

## Model Overview

HAL is written in FORTRAN and has evolved over many years to now offer a range of options for different species-specific biology and behaviours:

* Growth and development, based on temperature-dependent relationships, with default values for plaice
* Egg buoyancy, specified as either a vertical velocity or as a particle density, from which a vertical velocity is calculated based on a comparison of the particle and local water densities
* Diel vertical migration, either nocturnal or reverse migration
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

Young, E.F., Thorpe, S.E., Renner, A.H.H., Murphy, E.J. (2024). �Environmental and behavioural drivers of Antarctic krill distribution at the South Orkney Islands: a regional perspective�. Journal of Marine Systems, 241, https://doi.org/10.1016/j.jmarsys.2023.103920.

Trathan, P.N., Warwick-Evans, V., Young, E.F., Friedlaender, A., Kim, J.H., Kokubun, N. (2022). �The ecosystem approach to management of the Antarctic krill fishery - the �devils are in the detail� at small spatial and temporal scales�. Journal of Marine Systems, 225, https://doi.org/10.1016/j.jmarsys.2021.103598.






