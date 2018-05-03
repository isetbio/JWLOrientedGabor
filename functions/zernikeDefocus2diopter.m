function [diopters] = zernikeDefocus2diopter(zernikeDF, pupilDiameterMM)

% z4 (defocus) in microns to diopters conversion

% Inputs: zernike coefficient C2,0 (or J-number 4) and pupil diameter in
% mm 
% Output: spherical equivalent, or M in diopters

diopters = (16 * sqrt(3) * zernikeDF) / (pupilDiameterMM^2); 