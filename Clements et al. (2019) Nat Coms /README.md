# Experimental-data

## Introduction

Data from Clements, C., McCarthy, M., Blanchard, J. Early warning signals of recovery in complex systems. Nature Communications, 10:1681. These simulations rely on a stochastic version of the Mizer package (https://github.com/sizespectrum/mizer). This is currently in the pipeline to be released into the Mizer R package. It isn't currently available as source code.

## Structure

Species - Defines the species the data are for (Cod).

Year - the year of the data (from 2010 to 2100). 

Iteration - the unique simulation number (between 1 and 300 for each treatment)

Recovery.point - the inferred recovery point for each time series (see publication for further details). NA indicates no recovery, these are the control group where fishing pressures were not reduced

Rate - the time over which fishing pressures are reduced. 0 indicates no reduction in fishing pressure, 50 indicates fishing pressures were reduced from 2010 levels to zero fishing over a 50 year period.

mean.size - the mean body size of the cod population at each time step

sd.size - the standard deviation of body size of the cod population at each time step

biomass - the biomass of the cod stocks at each time step

## Treatments

Treatments consisted of a reduction in fishing pressures at various rates:

0 = control group

10 = fishing pressure reduced to zero from 2010 levels over 10 years
20 = fishing pressure reduced to zero from 2010 levels over 10 years
30 = fishing pressure reduced to zero from 2010 levels over 10 years
40 = fishing pressure reduced to zero from 2010 levels over 10 years
50 = fishing pressure reduced to zero from 2010 levels over 10 years

## Set up details

### Modelling

Here we employed a previously developed multispecies dynamic size spectrum model of the North Sea with 12 interacting species and a background resource community, which has been shown to provide realistic size spectra and population dynamics. This is a dynamic continuous time and size model based on the McKendrick von Foerster equations which is discretized for numerical approximation using finite upwind differencing methods24. The timestep was discretized at 0.1 years to ensure the integration routine was satisfactory, and size was discretized to 100 size classes. The data are yearly outputs from the model reflecting the yearly monitoring of fisheries stocks. Although the underlying model is deterministic (and is available in the R package ‘mizer’), as in Blanchard et al. we use the stochastic version that includes a log-normal error term on the recruitment for all 12 fish species. The model outputs included abundance and biomass of each species and their body size distributions through time, however in the analyses (and data presented here) we considered only the dynamics of cod. 

In line with Blanchard et al.24 a burn-in period of 300 years (between 1667 and 1967AD) was used to ensure the model reached equilibrium before carrying out time-varying fishing simulations. After this period, fishing was implemented at recorded historic levels for all of the 12-species in the model between 1967 and 2010 (Supplementary Fig. 1)24. To provide an example of an overexploited marine system, we held fishing pressures constant between 2010 indicators into a single metric of risk by summing them at each time point, as in Drake and Griffen (2010). Thus in total we tested 15 different metrics, composed of every unique combination of one to four indicators.

Subsequently, we simulated a range of post collapse scenarios (2040 onwards), where fishing pressures either remain constant (producing a persistent collapsed community) or decreased linearly across all species until fishing pressures reached 0, which allowed the system to recover from collapse. Because the rate at which pressures on a system changes can alter the prevalence of EWSCs37, we simulated five different rates of decline in the strength of fishing pressure (declining to no fishing linearly over a 10, 20, 30, 40, or 50 year period). Each of these six treatments (5 recovery and one control where fishing pressures were held at 2010 levels) were simulated 300 times over a period from 1667 to 2200, giving a total of 1800 simulations, of which 1500 recover and 300 do not. To avoid any artefacts of the data generation process, these 1800 simulations were split into three groups of equal size containing 100 simulations of each of the 6 rates of change in fishing pressure, and the stochastic simulations were initialized with the random number generators set at differing start points. We then focus our search for early warning signals solely on the simulated dynamics of cod, presenting the full community dynamics in the supplementary information.


## References

Clements, C., McCarthy, M., Blanchard, J. Early warning signals of recovery in complex systems. Nature Communications, 10:1681.