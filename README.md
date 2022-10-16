# Fast simulation AFP - C++
The modelling of Cherenkov based detectors is traditionally done using Geant4 toolkit. In this code, we do another method based on C++ programming language. 
As an example, we take one of the Forward Proton Detectors at the CERN LHC - ATLAS Forward Proton (AFP) Time-of-Flight, which is used to reduce the background from multiple proton-proton collisions in soft and hard diffiractive events. We describe the technical details of the fast Cherenkov model of photon generation and transportation through the optical part of the ToF detector. 
The fast simulation is revealed to be faster than the corresponding Geant4 simulation, and provides similar results concerning length and time distributions of photons. The study is meant as the first step in a construction of a building kit allowing creation of a fast simulation of an arbitrary shaped optical part of detectors.

Moreover, the process is accelerated and parallelized by using the directive OpenMP (pragma omp parallel). 
