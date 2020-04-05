function [OG, scenes, sparams] = getSceneAndStimuli()


% Temporal properties of one trial
tStep             = 0.002;                % Time step for optical image sequence (seconds)
sparams.tsamples  = (0:tStep:0.054*4);      % seconds
% sparams.timesd    = 1000.00;              % sd of temporal Gaussian window

% Scene field of view
sparams.sceneFOV  = 2;   % scene field of view in degrees (diameter)
sparams.freqCPD   = 4;   % Gabor spatial frequency (cpd)
sparams.gausSDdeg = sparams.sceneFOV/4; % Gabor SD in degrees of visual angle

% Unit converters
deg2fov = 1/sparams.sceneFOV;
fov2deg = sparams.sceneFOV;

% Gabor parameters
sparams.gabor           = harmonicP;                   % Standard Gabor
sparams.gabor.ang       = (pi/180)* 15;                % Gabor orientation (radians) - question: what is 0??
sparams.gabor.freq      = fov2deg*sparams.freqCPD;     % Spatial frequency (cycles/FOV)
sparams.gabor.contrast  = 1;                           % Presumably michelson, [0 1]
sparams.gabor.GaborFlag = sparams.gausSDdeg*deg2fov;   % Gaussian window

% Make a default Optical Image Sequence. We need some of the
% parameters from this. It will be overwritten as we change the
% contrast and/or optics.
[OG, scenes] = ogStimuli(sparams);
