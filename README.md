# JWLOrientedGabor

Welcome to the JWLOrientedGabor repository!

This repository runs on the ISETBIO toolbox and is created by the NYU ISETBIO team (JWL stands for Jonathan Winawer Lab): Eline Kupers & Jonathan Winawer. This project is funded by the NEI grant: 'Linking brain and behavior 'around' the visual field' by Marisa Carrasco (NYU) and Jonathan Winawer (NYU)

The paper related to this project can be found on bioRxiv (https://www.biorxiv.org/content/10.1101/434514v2) and is now accepted in PLOS Computational Biology!

## Short term goals
This repository creates a computational observer model, containing multiple
stages of the front-end of the visual system for a particular visual scene 
or psychophysical experiment.

The current state of this repository construct a scene with 2 Gabor stimuli 
(oriented clockwise or counter-clockwise), that goes through the optics and
gets sampled by a small cone mosaic patch. The cone absorptions are then used 
to simulate a 2-AFC orientation discrimination task with a linear classifier.

## Long term goals
Eventually, we would like to model more stages of the visual pathway:
For example: 
* cone currents (temporal integration);
* bipolar (temporal and spatial integration);
* RGC density (temporal and spatial integration);
* cortical responses (Receptive fields, cortical magnification); 
* other behavioral tasks

We can then vary these stages systematically and compute how much of these
variations can account for variations in the computational observer
performance for the given task as a function of the visual field position of
the stimulus. 

## Dependencies
This repository is MATLAB-based (v9.1) and depends on the following toolboxes:
* ISETBIO (https://github.com/isetbio/isetbio/, works at least until version commit e620d97 - Feb 12, 2019)

MATLAB Toolboxes:
* Bioinformatics Toolbox, 4.7
* Computer Vision System Toolbox, 7.2
* Control System Toolbox, 10.1
* Curve Fitting Toolbox, 3.5.4
* DSP System Toolbox, 9.3
* Database Toolbox, 7.0
* Global Optimization Toolbox, 3.4.1
* Image Acquisition Toolbox, 5.1
* Image Processing Toolbox, 9.5
* Instrument Control Toolbox, 3.10
* Mapping Toolbox, 4.4
* Neural Network Toolbox, 9.1
* Optimization Toolbox, 7.5
* Parallel Computing Toolbox, 6.9
* Partial Differential Equation Toolbox, 2.3
* Signal Processing Toolbox, 7.3
* Statistics and Machine Learning Toolbox, 11.0
* Symbolic Math Toolbox, 7.1
* System Identification Toolbox, 9.5
* Wavelet Toolbox, 4.17

## Example scripts
```
% Run default observer model for CW/CCW Gabors at 4.5 deg eccentricity, 4 cpd, with small fixational eyemovements
runComputationalObserverModel('default')
```

