# JWLOrientedGabor
Simulation of an experiment by Marisa Carrasco and Jon Winawer measuring orientation discrimination performance for an acrhomatic, peripheral  Gabor. 

The Gabor is at 6 deg eccentricity, 6 cpd, 1.5 deg windwo, and +/- 20 deg from vertical. 

Trials are simulated as a 100 ms static exposure with scene jitter created by eye movements (sampling from a 2D isotropic Gaussian distribution, std = 2 cones).

50 trials are simulated for each of the two stimulus orientations. A classifier predicts the orientation based on the cone responses, summed across the 100 ms within each cone. Classification is also computed for several simplified RGC models. The RGC models are 2D isotropic center-surround Gaussians, convolved with cone outputs and subsampled. The RGC RFs are on-center, and the classes differ only on the size of the RF and subsampling. The subsampling is proportional to the RF size (if the center has STD of 2 cones, then we subsample every 2nd RGC in each dimension, x/y). 

The repository was initialized by copying all files in the ISETBIO repository HLColorOpponentGabor. The code was then altered to reflect differences in the experiments, including achromatic orientation discrimination (this repository) versus chromatic discrimnation (HLColorOpponentGabor), and to include eye movement traces, RGC outputs, and to *not* pool responses across cones.
