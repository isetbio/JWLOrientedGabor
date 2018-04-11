function spacing = getBanks1991ConeSpacing(eccentricity, varargin)

% Syntax
%   spacing = getBanks1991ConeSpacing([eccentricity], ['unit'], ['m'])

% Get cone spacing (in arc minutes), assuming a regular triangular
% arrangement based on Figure 1 of Banks et al. (1991): Peripheral spatial
% vision: limited imposed by optics, photoreceptors and receptor pooling.
% J. Opt. Soc Am.

% Inputs
% eccentricity:     integer defining the eccentricity of interest
% unit:             [=optional], string defining whether you want arc min or m  

% Outputs
% spacing:          integer spacing of cones in arc minutes [default] or
%                   meters [if requested]
%
% Example:
%  spacing = getBanks1991ConeSpacing(5, 'unit', 'm')

%% Handle input arguments
p = inputParser;
p.KeepUnmatched = true;

% Required
validEccentricities = [0, 2.5, 5, 10, 20, 40];
p.addRequired('eccentricity',@isvector);
p.addParameter('unit', 'arcmin')

% Parse
p.parse(eccentricity, varargin{:});

% Data is inferred from figure, not actual from a publicly shared dataset
coneSpacing = [0.5, 1.2, 1.8, 2.25, 3.0, 3.5];

idx = (p.Results.eccentricity == validEccentricities);
if any(idx)
     spacingArcMin = coneSpacing(idx);
else
    error('(getBanks1991ConeSpacing): No data is available for requested eccentricity');
end
    
if strcmp(p.Results.unit, 'm')
    degPerArcmin = 1/60;
    focalLength = 0.017; % meters, isetbio default
    spacing = tan(deg2rad(spacingArcMin * degPerArcmin)/2) * focalLength; % spacing in meters
else
    spacing = spacingArcMin; % spacing in arc minute
end    



return 