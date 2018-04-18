function propCovered = getBanks1991ConeCoverage(eccentricity, varargin)

% Syntax
%   propCovered = getBanks1991ConeCoverage(eccentricity, ['eccentricityUnits'], ['deg'])

% Get proportion covered by cones for different eccentricities (ratio), by fitting
% a straight line in semi-log space to data in Figure 2 of Banks et al.
% (1991): Peripheral spatial vision: limited imposed by optics,
% photoreceptors and receptor pooling. J. Opt. Soc Am.

% Inputs
% eccentricity:         integer defining the eccentricity of interest, can 
%                       only be integer: 0, 2.5, 5, 10, 20, or 40 deg
% eccentricityUnits:    [=optional], string defining unit of input
%                       eccentricity: 'deg' [default], 'm','mm','um') will 
%                       be converted to 'deg'.


% Outputs
% propCovered:          integer of proportion covered by the cones for the
%                       requested eccentricity
%
% Example:
%  propCovered = getBanks1991ConeCoverage(4)

%% Handle input arguments
p = inputParser;
p.KeepUnmatched = true;

% Required and added
p.addRequired('eccentricity',@isvector);
p.addParameter('eccentricityUnits','deg',@ischar);

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

% set manually:
allEccentricities = 0:0.5:50;
% coverage = (exp(-1/40*allEccentricities));

banksData  = [1, 0.8, 0.5, 0.4,0.3, 0.25];
banksEccen = [0, 2, 5, 10, 20, 40];

% or with fitting function:
f = fit(banksEccen', banksData','exp2');

coverage = f.a*exp(f.b*allEccentricities) + f.c*exp(f.d*allEccentricities);
coverage(1) = 1; % reset fovea to 1.

% Plot for debugging
figure(1); clf,
plot(allEccentricities, coverage, 'o-'); set(gca, 'YScale', 'log', ...
    'XScale', 'linear', 'XLim', [0 50], 'YLim', [0.01 1]);
hold all; scatter(banksEccen, banksData, 80','k');

idx = (eccDeg == allEccentricities);
if any(idx)
     propCovered = sqrt(coverage(idx));
else
    error('(getBanks1991ConeCoverage): No data is available for requested eccentricity');
end
    

return 