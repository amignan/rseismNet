# rseismNet R Package
Earthquake Frequency-Magnitude Distribution & Network Statistics

A suite of statistical functions to describe the earthquake frequency
    magnitude distribution (FMD) statistics. This R package is based on the 2 following 
    rules: (1) the complete FMD side (m>=mc) is governed by the Gutenberg-Richter law 
    and depends on the fault network properties; (2) the incomplete FMD part (m<mc) is 
    governed by a detection function and depends on the seismic network properties 
    (with m the earthquake magnitude and mc the completeness magnitude). rseismNet also
    implements the Bayesian Magnitude of Completeness (BMC) method, which uses a seismic
    network prior to map mc.
