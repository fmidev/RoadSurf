# RoadSurf
Road weather model library

RoadSurf is a Fortran library for predicting road conditions. Instructions for using the library and physics documentation can be found in the included user manual

Written by: Markku Kangas, Marjo Hippi, Johanna Ruotsalainen, Martti Heikinheimo, Virve Karsisto, Mika Heiskanen and Leif Backman
1998-2023 Finnish Meteorological Institute


The library does not contain all of the features included in the original FMI road weather model. The included features are road surface temperature and storage term (water, ice, snow, deposit) calculation. Variables determined based on these values, such as friction and road condition, are excluded from the library. The library is designed to allow flexible implementation. It includes the subroutines for model initialization and calculation of the road surface temperature and storage terms, but the user can develop their own main program which implements them. The subroutines are separated to their own library to make it easier to keep the operational implementations up to date. Each user can have different kinds of input sources or even their own output variables. When the library is updated, it is easier to just update the library without the need to update the main program. One example of RoadSurf implementation is given in "examples" folder. The example utilizes C++ interface. The input and output of the model are handled by the C++, whereas the actual simulation is done with the Fortran library. The code provided has also an example for running the simulations in parallel to multiple points at the same time. 

## Licence

The library is published with MIT License

## How to contribute

Your contribution is very welcome, be it bug fixes or new features.

Small changes and bug fixes can be submitted via pull request. In larger contributions, premilinary plan is recommended (in GitHub wiki).

CLA is required in order to contribute. Please contact us for more information!

## Communication and Resources

You may contact us from following channels:

    Email: beta@fmi.fi
    Facebook: https://www.facebook.com/fmibeta/
    GitHub: issues

The physics of the library have been presented in the following publication:

Karsisto, V. E. 2024: RoadSurf 1.1: open-source road weather model library, Geosci. Model Dev., 17, 4837–4853, https://doi.org/10.5194/gmd-17-4837-2024 
