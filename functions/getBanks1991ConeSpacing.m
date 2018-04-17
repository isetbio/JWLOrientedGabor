function spacing = getBanks1991ConeSpacing(eccentricity, varargin)

% Syntax
%   spacing = getBanks1991ConeSpacing(eccentricity, ['eccentricityUnits'], ['deg'],['focalLength'],0.017,['spacingUnits'], ['m'])

% Get cone spacing (in arc minutes), assuming a regular triangular
% arrangement based on Figure 1 of Banks et al. (1991): Peripheral spatial
% vision: limited imposed by optics, photoreceptors and receptor pooling.
% J. Opt. Soc Am.

% Inputs
% eccentricity:         integer defining the eccentricity of interest, can 
%                       only be integer: 0, 2.5, 5, 10, 20, or 40 deg
% eccentricityUnits:    [=optional], string defining unit of input
%                       eccentricity: 'deg' [default], 'm','mm','um') will 
%                       be converted to 'deg'.
% focalLength:          length between cornea and retina in meters,
%                       necessary for calculating eccentricity based cone
%                       spacing. [isetbio default = 0.017 m)
% spacingUnits:         [=optional], string defining whether you want arc min or m  


% Outputs
% spacing:              integer spacing of cones in arc minutes [default] or
%                       meters [if requested with spacingUnits]
%
% Example:
%  spacing = getBanks1991ConeSpacing(5, 'spacingUnits', 'm')

%% Handle input arguments
p = inputParser;
p.KeepUnmatched = true;

% Required
validEccentricities = [0, 2.5, 5, 10, 20, 40];
p.addRequired('eccentricity',@isvector);
p.addParameter('eccentricityUnits','deg',@ischar);
p.addParameter('focalLength', 0.017, @isnumeric); % meters, isetbio default
p.addParameter('spacingUnits', 'arcmin', @ischar);

% Parse
p.parse(eccentricity, varargin{:});

% Get corrent eccenctricity units
switch (p.Results.eccentricityUnits)
    
    case 'm'
        eccDeg = p.Results.eccentricity*(1000/0.3);
    case 'mm'
        eccDeg = p.Results.eccentricity*(1/0.3);
    case 'um'
        eccDeg = p.Results.eccentricity*(0.001/0.3);
    case 'deg'
        eccDeg = p.Results.eccentricity;
    otherwise
        error('(getBanks1991ConeSpacing): Unknonwn units for eccentricity specified');
end

% Data is inferred from figure, not actual from a publicly shared dataset
coneSpacing = [0.5, 1.2, 1.8, 2.25, 3.0, 3.5];

idx = (eccDeg == validEccentricities);
if any(idx)
     spacingArcMin = coneSpacing(idx);
else
    error('(getBanks1991ConeSpacing): No data is available for requested eccentricity');
end
    
if strcmp(p.Results.spacingUnits, 'm')
    degPerArcmin = 1/60;
    focalLength = p.Results.focalLength;
    spacing = tan(deg2rad(spacingArcMin * degPerArcmin)/2) * focalLength; % spacing in meters
else
    spacing = spacingArcMin; % spacing in arc minute
end    



return 