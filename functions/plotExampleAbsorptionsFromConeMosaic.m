function plotExampleAbsorptionsFromConeMosaic(expName)

% Function to plot the absorptions of a stimulus. For now, one can pick
% 'conedensity' to get 3 different cone densities (fovea, 4.5 degree
% and 40 degree eccentricity along the temporal (?) side of left retina.
% Or pick 'conetypes', to plot a default mosaic versus three  mosaics with
% only one cone type (L, M or S).

%% 0. Check input parameters

p = inputParser;
p.KeepUnmatched = true;
p.addRequired('expName', @(x) (strcmp(x,'conetypes') | strcmp(x,'conedensity')));
p.parse(expName);

%% 1. Define other parameters

% Experiment params
expParams         = loadExpParams(expName, false);
whichEye          = 'left';
cparams.cmFOV     =  2; % degrees

if strcmp(expName, 'conedensity')
    eccentricities    = [0, 4.5, 40]; % deg
    nrPlots = length(eccentricities);
    figttl = 'Figure 8A - 2D absorptions for 3 cone density levels';
else
    eccentricities    = 4.5; % deg
    cparams.spatialDensity = expParams.cparams.spatialDensity;
    nrPlots = size(cparams.spatialDensity,2);
    figttl = 'Figure XX - 2D absorptions for 3 single cone type mosaics + default mosaic';
end

% Convertion deg to m
deg2m  = 1/3 * 0.001; % 3 deg per mm, .001 mm per meter

% Set up figure
figure(1); hold all;
set(gcf, 'Color', 'w', 'Position', [1000, 681, 1173, 657], 'NumberTitle', 'off', 'Name', figttl);

% Set up scene params
tStep             = 0.002;                % Time step for optical image sequence (seconds)
sparams.tsamples  = (0:tStep:(0.054));      % seconds

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

% Add defocus
sparams.oi = oiDefocus(0);

%% 2. Create Optical Image Sequence
OG = ogStimuli(sparams);

%% 3. Define cone mosaic parameters
for idx = 1:nrPlots
    
    if strcmp(expName, 'conedensity')
        % Specify retinal location where stimulus is presented
        cparams.eccentricity = eccentricities(idx);            % Visual angle of stimulus center, in deg
        spatialDensity       = [0 0.6 0.3 0.1];                % Reset to default mosaic: Blank, L, M, S cone type ratio within cone mosaic
    elseif strcmp(expName, 'conetypes')
        cparams.eccentricity = eccentricities;                 % Visual angle of stimulus center, in deg
        spatialDensity       = cparams.spatialDensity(idx,:);  % Blank, L, M, S cone type ratio within cone mosaic
    end
    
    % Specify retinal location where stimulus is presented
    cparams.polarAngle   = expParams.polarAngle;   % Polar angle (radians): 0 is right, pi/2 is superior, pi is left, 3*pi/2 inferior
      
    % Compute x,y position in m of center of retinal patch from ecc and angle
    [x, y] = pol2cart(cparams.polarAngle, cparams.eccentricity);
    x = x * deg2m;  y = y * deg2m;
    
    %% 4. Create cone mosaic
    cMosaic = coneMosaic('center', [x, y], 'whichEye', whichEye, 'spatialDensity', spatialDensity);
    
    % Set the field of view (degrees)
    cMosaic.setSizeToFOV(cparams.cmFOV);
    
    if strcmp(expName, 'conedensity')
        % Set proportion covered
        propCovered = getBanks1991ConeCoverage(cparams.eccentricity, 'verbose', false);
        cMosaic.pigment.pdWidth  = cMosaic.pigment.width*propCovered;
        cMosaic.pigment.pdHeight = cMosaic.pigment.height*propCovered;
    end
    
    % Create empty array of eye movements
    emPaths = zeros(expParams.nTrials, length(sparams.tsamples), 2);
    
    %% 5. Compute absorptions for one stimulus
    absorptions = cMosaic.compute(OG(1), 'currentFlag', false, 'emPaths', emPaths, 'seed', 1);
    
    if strcmp(expName, 'conedensity')
        ttl = sprintf('Cone density %d', eccen2density(cMosaic, 'deg'));
    elseif strcmp(expName, 'conetypes')
        ttl = sprintf('L:M:S ratio: %1.1f:%1.1f:%1.1f', spatialDensity(2), spatialDensity(3), spatialDensity(4));
    end
    
    %% 6. Visualize
    subplot(1,nrPlots, idx); hold all;
    imagesc(squeeze(absorptions(1,:,:,1))./2);
    colormap gray; axis image; axis off; colorbar;
    title(ttl ,'Fontsize',12)
    
end

return