function [diopters] = zernickeDefocus2diopter(zernikeDF, pupilRadiusMM)


diopters = 4*pi*sqrt(3) * zernikeDF / (pi* pupilRadiusMM^2); 