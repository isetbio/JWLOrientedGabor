function [xUnits, colors, labels, M] = loadWeibullPlottingParams(expName)

% Define Weibull plotting parameters for a specific experiment

%       [xUnits, colors, labels, M] = loadWeibullPlottingParams(expName)

% INPUTS: 
%   expName        : String defining which type of experiment parameters to use. 
%                       Choose from 'default',
%                                   'eyemov',
%                                   'eyemovEnhanced',
%                                   'coneDensity',
%                                   'defocus',
%                                   'ConeTypes',
%                     (default= 'default')

%% Check input arguments
if isempty(expName) || ~exist('expName', 'var')
    expName = 'default'; 
end

% Get general condition parameters
expParams                = loadExpParams(expName, false);
M                        = [];

% Define plotting parameters
switch lower(expName)
    case 'default'
        colors              = [0 0 0];
        labels              = {'Fit','data'};
        xUnits              = linspace(min(expParams.contrastLevels),max(expParams.contrastLevels), 100);
        
    case 'conetypes'
        colors              = copper(size(expParams.cparams.spatialDensity,1));
        labels              = {'LMS default cone mosaic','L only cone mosaic','M only cone mosaic','S only cone mosaic'};
        xUnits              = linspace(min(expParams.contrastLevels),max(expParams.contrastLevels), 100);
        
    case 'conetypesmixed'
        nrMixtures          = size(expParams.cparams.spatialDensity,1);
        colors              = jet(nrMixtures);
        labels              = cell(nrMixtures,1);
        for ii = 1:nrMixtures
            labels{ii}      =  sprintf('LMS cone ratio = %1.1f:%1.1f:%1.1f', expParams.cparams.spatialDensity(ii,2:4));
        end
        xUnits              = linspace(min(expParams.contrastLevels),max(expParams.contrastLevels), 100);
        
    case 'eyemov'
        colors              = [0 0 0;  0.5000, 1.0000, 0.5000; 1 0 0];
        labels              = cell(1,length(colors));
        for emIdx = 1:length(colors)
            thisCondition = expParams.eyemovement(:,emIdx)';
            if all(thisCondition == [0 0])
                labels{:,emIdx} = 'No eyemovements';
            elseif  all(thisCondition == [1 0])
                labels{:,emIdx} = 'Drift';
            elseif  all(thisCondition == [0 1])
                labels{:,emIdx} = 'MS';
            elseif  all(thisCondition == [1 1])
                labels{:,emIdx} = 'Drift and MS';
            end
        end
        xUnits          = linspace(min(expParams.contrastLevels),max(expParams.contrastLevels), 100);
        
    case 'eyemovenhanced'
        colors              = copper(size(expParams.eyemovement,2));
        labels              = cell(1,length(colors));
        for emIdx = 1:length(colors)
            thisCondition = expParams.eyemovement(:,emIdx)';
            if all(thisCondition == [0 0 0])
                labels{:,emIdx} = 'No eyemovements';
            elseif  all(thisCondition == [2 0 0])
                labels{:,emIdx} = '2xTremor';
            elseif  all(thisCondition == [0 2 0])
                labels{:,emIdx} = '2xDrift';
            elseif  all(thisCondition == [2 2 0])
                labels{:,emIdx} = '2xTremor+Drift';
            elseif  all(thisCondition == [0 0 2])
                labels{:,emIdx} = '2xMicrosaccades';
            end
        end
        
        xUnits          = linspace(min(expParams.contrastLevels),max(expParams.contrastLevels), 100);
        
    case {'conedensity','conedensitynoeyemov','eccbasedcoverage'}
        
        colors              = jet(length(expParams.eccentricities));
        
        % Change x labels to density
        whichEye          = 'left';
        cparams.cmFOV     =  2; % degrees
        
        % Convertion deg to m
        deg2m  = 0.3 * 0.001; % .3 deg per mm, .001 mm per meter
        
        % Predefine density vector
        allDensity = nan(length(expParams.eccentricities),1);
        labels = cell(length(expParams.eccentricities),1);
        
        for ec = expParams.eccentricities
            % Specify retinal location where stimulus is presented
            cparams.eccentricity = ec;                     % Visual angle of stimulus center, in deg
            cparams.polarAngle   = expParams.polarAngle;   % Polar angle (radians): 0 is right, pi/2 is superior, pi is left, 3*pi/2 inferior
            
            % Compute x,y position in m of center of retinal patch from ecc and angle
            [x, y] = pol2cart(cparams.polarAngle, cparams.eccentricity);
            x = x * deg2m;  y = y * deg2m;
            cMosaic = coneMosaic('center', [x, y], 'whichEye', whichEye);
            
            % Set the field of view (degrees)
            cMosaic.setSizeToFOV(cparams.cmFOV);
            allDensity(ec==expParams.eccentricities,:) = eccen2density(cMosaic, 'deg');
            
            labels{ec==expParams.eccentricities} = sprintf('%1.3f x10^4 cells/deg2', allDensity(ec==expParams.eccentricities)/10.^4);
        end
        
        xUnits              = linspace(min(expParams.contrastLevels),max(expParams.contrastLevels), 100); % 
        M = allDensity;
        
    case 'defocus'
        
        for df = expParams.defocusLevels
            
            % compute defocus
            pupilRadiusMM = 1.5; % mm
            M(df==expParams.defocusLevels) = wvfDefocusMicronsToDiopters(df,pupilRadiusMM*2); % convert to diopters
            labels{df==expParams.defocusLevels} = sprintf('%2.2f Diopters of Defocus',M(df==expParams.defocusLevels));
        end
        
        colors              = jet(length(expParams.defocusLevels));
        xUnits              = linspace(min(expParams.contrastLevels),max(expParams.contrastLevels), 100);
        
        
    case 'conetypeseccen'
        
        colors              = jet(length(expParams.eccentricities));
        
        % Change x labels to density
        whichEye          = 'left';
        cparams.cmFOV     =  1; % degrees
        
        % Convertion deg to m
        deg2m  = 0.3 * 0.001; % .3 deg per mm, .001 mm per meter
        
        % Predefine density vector
        allDensity = nan(length(expParams.eccentricities),1);
        labels = cell(length(expParams.eccentricities),1);
        
        for ec = expParams.eccentricities
            % Specify retinal location where stimulus is presented
            cparams.eccentricity = ec;                     % Visual angle of stimulus center, in deg
            cparams.polarAngle   = expParams.polarAngle;   % Polar angle (radians): 0 is right, pi/2 is superior, pi is left, 3*pi/2 inferior
            
            % Compute x,y position in m of center of retinal patch from ecc and angle
            [x, y] = pol2cart(cparams.polarAngle, cparams.eccentricity);
            x = x * deg2m;  y = y * deg2m;
            cMosaic = coneMosaic('center', [x, y], 'whichEye', whichEye);
            
            % Set the field of view (degrees)
            cMosaic.setSizeToFOV(cparams.cmFOV);
            allDensity(ec==expParams.eccentricities,:) = eccen2density(cMosaic, 'deg');
            
            labels{ec==expParams.eccentricities} = sprintf('%1.2f x10^4 cells/deg2', allDensity(ec==expParams.eccentricities)/10.^4);
        end
        
        xUnits              = linspace(min(expParams.contrastLevels),max(expParams.contrastLevels), 100); % 
        M = allDensity;
        
end

