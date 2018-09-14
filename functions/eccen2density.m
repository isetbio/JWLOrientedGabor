function density = eccen2density(cMosaic, unit)

% density = eccen2density(mosaic)

% function to compute the cone density of a retinal patch located defined 
% in cMosaic struct eccentricity

% Example:
% x = 0;  y = 0;
% whichEye = 'left';
% cMosaic = coneMosaic('center', [x, y], 'whichEye', whichEye);
% density = eccen2density(cMosaic,'m')


% Check inputs
if nargin < 1 
    error('cMosaic needs to be defined')
end

if ~exist('unit','var') || isempty(cMosaic)
    unit = 'm'; % default of isetbio
end


%% Compute cones per requested unit
density = 1./(prod(cMosaic.patternSampleSize));

if strcmp(unit,'mm')
     density = density*10.^-6;
end

if strcmp(unit,'deg')
    densityMM2 = density*10.^-6;
    densityDeg2 = densityMM2/((1/.3).^2);
    density = densityDeg2;
end

return
