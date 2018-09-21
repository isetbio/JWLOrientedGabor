# JWLOrientedGabor

Welcome to the JWLOrientedGabor repository!

This repository runs on the ISETBIO toolbox and is created by the NYU ISETBIO team (JWL stands for Jonathan Winawer Lab): Eline Kupers & Jonathan Winawer. This project is funded by the NEI grant: 'Linking brain and behavior 'around' the visual field' by Marisa Carrasco (NYU) and Jonathan Winawer (NYU)

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
* RGC density (temporal and spatial integration;
* cortical responses (Receptive fields, cortical magnification); 
* other behavioral tasks

We can then vary these stages systematically and compute how much of these
variations can account for variations in the computational observer
performance for the given task as a function of the visual field position of
the stimulus. 

## More to follow here...
