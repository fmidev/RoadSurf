# RoadSurf
Road weather model library
Contributors: Markku Kangas, Marjo Hippi, Johanna Ruotsalainen, Martti Heikinheimo, Virve Karsisto, Mika Heiskanen and Leif Backman

RoadSurf is a Fortran library for predicting road conditions. Instructions for using the library and physics documentation can be found in the included user manual

The library does not contain all of the features included in the original FMI road weather model. The included features are road surface temperature and storage term (water, ice, snow, deposit) calculation. Variables determined based on these values, such as friction and road condition, are excluded from the library. The library is designed to allow flexible implementation. It includes the subroutines for model initialization and calculation of the road surface temperature and storage terms, but the user can develop their own main program which implements them. The subroutines are separated to their own library to make it easier to keep the operational implementations up to date. Each user can have different kinds of input sources or even their own output variables. When the library is updated, it is easier to just update the library without the need to update the main program. One example of RoadSurf implementation is given in "examples" folder. The example utilizes C++ interface. The input and output of the model are handled by the C++, whereas the actual simulation is done with the Fortran library. The code provided has also an example for running the simulations in parallel to multiple points at the same time. 

## Licence

The library is published with GNU Lesser General Public License v2.1.

## How to contribute

Found a bug? Want to implement a new feature? Your contribution is very welcome!

Small changes and bug fixes can be submitted via pull request. In larger contributions, premilinary plan is recommended (in GitHub wiki).

CLA is required in order to contribute. Please contact us for more information!

## Communication and Resources

You may contact us from following channels:

    Email: beta@fmi.fi
    Facebook: https://www.facebook.com/fmibeta/
    GitHub: issues



