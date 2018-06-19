# Dendro-GR
## Numerical relativity with octree based wavelet adaptive mesh refinement. 


* Dendro-GR portable and highly-scalable algorithm and framework that
targets problems in the astrophysics and numerical relativity
communities. This framework combines Dendro octree with Wavelet Adaptive 
Multiresolution (WAMR) and a physics module to solve the Einstein equations of 
general relativity in BSSN formulation.

* Wavelet adaptive multiresolution is used to provide local, 
fine-grain adaptivity, based on the fast wavelet transform of interpolating
wavelets.
* The goal of this work is to perform advanced, massively parallel 
numerical simulations of Intermediate Mass Ratio Inspirals (IMRIs)
of binary black holes with mass ratios on the order of 100:1. These studies
will be used to generate waveforms for use in LIGO data analysis and to calibrate
semi-analytical approximate methods. 

* This advanced framework is designed to easily accommodate many existing algorithms in astrophysics
for %compressible plasma dynamics and radiation hydrodynamics.

* We have designed novel algorithms to enable efficient simulations for such 
experiments and demonstrate excellent weak scalability up to 131K cores on ORNL's Titan.
