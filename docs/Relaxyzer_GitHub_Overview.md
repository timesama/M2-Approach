# Relaxyzer Software

Relaxyzer was initially created as a tool to speed up the analysis of NMR data obtained from Resonance Systems equipment. It has since grown into a broader and convenient software package for processing, visualizing, fitting, and comparing several types of NMR experiments.

## Key Features

### SE
Analysis of FID, Free Induction Decay, Solid Echo, and MSE data.

- Solid content calculation
- Second moment calculation for advanced analysis
- T2* calculation
- Activation energy calculation from T2* values measured at different temperatures

### DQ
Analysis of Double Quantum NMR data using the M2 approach.

- DQ distribution analysis
- T2* recalculation from the DQ distribution

### DQ(Temp)
Temperature-dependent comparison of T2* distributions.

- Overlay and comparison of distributions measured at different temperatures
- Support for thermal trend analysis

### T1T2
Analysis of T1 and T2 relaxation curves.

- Calculation of relaxation times and corresponding amplitudes
- Deconvolution into up to three components
- Support for regular T1/T2 curves
- Support for T1 curves recorded for different FID regions
- Support for T1 data measured from FFC profiles obtained with Stelar equipment

### DQMQ
Analysis of Double Quantum NMR data.

- Build-up curve construction
- Normalized double quantum signal, nDQ, calculation
- Dres distribution fitting

### Goldman-Shen
Domain size estimation using Goldman-Shen analysis.

## Documentation

Detailed user and developer documentation is available in:

`docs/Relaxyzer_User_and_Code_Documentation.tex`

This document describes the workflows, file formats, mathematical formulas, save/load behavior, and implementation details for each program tab.
