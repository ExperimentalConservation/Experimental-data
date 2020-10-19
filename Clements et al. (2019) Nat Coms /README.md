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

### to finish

.....

### Nutrients

The medium solution contained 0!1 g of crushed protozoa pellets (Carolina Biological Supply, USA) per litre of Chalkley’s medium (Thompson, Rhodes & Pettman 1988).

### Communities
S
There....

## References

Clements, C., McCarthy, M., Blanchard, J. Early warning signals of recovery in complex systems. Nature Communications, 10:1681.