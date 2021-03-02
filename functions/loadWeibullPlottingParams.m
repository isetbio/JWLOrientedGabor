function [xUnits, colors, labels, xThresh, lineStyles] = loadWeibullPlottingParams(expName)
% Define Weibull plotting parameters for a specific experiment
%
%       [xUnits, colors, labels, M] = loadWeibullPlottingParams(expName)
%
% INPUTS: 
%   expName        : String defining which type of experiment parameters to use. 
%                       Choose from 'default', [=default]
%                                   'conetypes',
%                                   'conetypesmixed',
%                                   'eyemov',
%                                   'eyemovnophaseshift',
%                                   'eyemovenhanced' [NB: currently not implemented as an experiment] 
%                                   'idealobserver',
%                                   'defaultnophaseshift',
%                                   'defocus',
%                                   'conedensity',
%                                   'conetypeseccen', [NB: currently not implemented as an experiment] 
%  
%
% Example:
% [xUnits, colors, labels, xThresh, lineStyles] = loadWeibullPlottingParams('conetypes')

%% Check input arguments
if isempty(expName) || ~exist('expName', 'var')
    expName = 'default'; 
end

% Get general condition parameters
expParams    = loadExpParams(expName, false);
xThresh      = []; % x units for plotting thresholds 

% Define plotting parameters
switch lower(expName)
    case 'default'
        colors          = [0 0 0];
        labels          = {'Fit','data'};
        xUnits          = linspace(min(expParams.contrastLevels),max(expParams.contrastLevels), 200);
        lineStyles      = repmat({'-'},size(colors));

    case 'conetypes'
        colors          = [0 0 0; 1 0 0; 0 0.5 0.25; 0 0 1];
        labels          = {'LMS default cone mosaic','L only cone mosaic','M only cone mosaic','S only cone mosaic'};
        xUnits          = linspace(min(expParams.contrastLevels),max(expParams.contrastLevels), 200);
        lineStyles      = repmat({'-'},size(colors));

    case 'conetypesmixed'
        nrMixtures      = size(expParams.cparams.spatialDensity,1);
        colors          = [linspace(1,0,nrMixtures); linspace(0,.5,nrMixtures); linspace(0, .25, nrMixtures)]';
        labels          = cell(nrMixtures,1);
        for ii = 1:nrMixtures
            labels{ii}  =  sprintf('LMS cone ratio = %1.1f:%1.1f:%1.1f', expParams.cparams.spatialDensity(ii,2:4));
        end
        xUnits          = linspace(min(expParams.contrastLevels),max(expParams.contrastLevels), 200);
        xThresh         = 100:-10:0;
        lineStyles      = {'-',':', ':', ':', '-'};

    case 'eyemov'
        colors          = [0 0 0; 1 0 0; 0.5000, 1.0000, 0.5000];
        lineStyles      = {'-', '-', ':'};
        labels          = cell(1,length(colors));
        xUnits          = linspace(min(expParams.contrastLevels),max(expParams.contrastLevels), 200);

        for emIdx = 1:size(expParams.eyemovement,2)
            thisCondition = expParams.eyemovement(:,emIdx)';
            if all(thisCondition == [0 0])
                labels{:,emIdx} = 'No eye movements';
            elseif  all(thisCondition == [1 0])
                labels{:,emIdx} = 'Drift';
            elseif  all(thisCondition == [0 1])
                labels{:,emIdx} = 'MS';
            elseif  all(thisCondition == [1 1])
                labels{:,emIdx} = 'Drift and MS';
            end
        end

    case 'eyemovnophaseshift'
        colors          = [0 0 0; 1 0 0; 0.5000, 1.0000, 0.5000];
        lineStyles      = {'-', '-', ':'};
        xUnits          = linspace(min(expParams.contrastLevels),max(expParams.contrastLevels), 1000);
        labels          = {'No eye movements, no stim phase shift', ...
                            'Drift, no stim phase shift', ...
                            'Drift and MS, no stim phase shift'};
        
    case 'idealobserver'
        xUnits          = linspace(min(expParams.contrastLevels),max(expParams.contrastLevels), 1000);
        colors          = [0 0 0; 0 0 0; 0.5000, 1.0000, 0.5000; 1 0 0];
        lineStyles      = {'-', ':', '-', '-'};
        labels          = {'Ideal observer (Analytical)', ...
                            'Ideal observer (Simulation, 200 trials per stimulus class)', ...
                            'Computational observer (SVM Classifier, 200 trials per stimulus class)', ...
                            'Computational observer (SVM Classifier, 800 trials per stimulus class)'};
              
     case 'defaultnophaseshift'
        xUnits          = linspace(min(expParams.contrastLevels),max(expParams.contrastLevels), 1000);
        colors          = [0 0 0];
        lineStyles      = {'-'};
        labels          = {'No eye movements, no phase shift'};
        
    case 'eyemovenhanced'
        colors          = copper(size(expParams.eyemovement,2));
        labels          = cell(1,length(colors));
        xUnits          = linspace(min(expParams.contrastLevels),max(expParams.contrastLevels), 100);
        lineStyles      = repmat({'-'},size(colors));
        
        for emIdx = 1:length(colors)
            thisCondition = expParams.eyemovement(:,emIdx)';
            if all(thisCondition == [0 0 0])
                labels{:,emIdx} = 'No eye movements';
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

    case {'conedensity','conedensitynoeyemov','eccbasedcoverage'}
        
        colors          = jet(length(expParams.eccentricities));
        xUnits          = linspace(min(expParams.contrastLevels),max(expParams.contrastLevels), 100); % 
        
        % Get parameters to compute cone density levels
        whichEye        = 'left';
        cparams.cmFOV   =  2; % degrees
        
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
        
        xThresh         = allDensity;
        lineStyles      = repmat({'-'},size(colors));

    case 'defocus'
        
        colors          = jet(length(expParams.defocusLevels));
        xUnits          = linspace(min(expParams.contrastLevels),max(expParams.contrastLevels), 200);
        lineStyles      = repmat({'-'},1,length(expParams.defocusLevels));

        % Compute defocus levels
        for df = expParams.defocusLevels            
            pupilRadiusMM = 1.5; % mm
            xThresh(df==expParams.defocusLevels) = wvfDefocusMicronsToDiopters(df,pupilRadiusMM*2); % convert to diopters
            labels{df==expParams.defocusLevels} = sprintf('%2.2f Diopters of Defocus',xThresh(df==expParams.defocusLevels));
        end
      
    case 'conetypeseccen'
        
        colors           = jet(length(expParams.eccentricities));
        xUnits           = linspace(min(expParams.contrastLevels),max(expParams.contrastLevels), 200); % 
        lineStyles       = repmat({'-'},size(colors));

        % Get parameters to compute cone density levels
        whichEye         = 'left';
        cparams.cmFOV    =  1; % degrees
        
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
        
        xThresh           = allDensity;
        
     case 'rgcratios'
        cone2RGCRatios  = 1:5;
        colors          = parula(length(cone2RGCRatios)+1);
        xUnits          = linspace(min(expParams.contrastLevels),max(expParams.contrastLevels), 200);
        lineStyles      = repmat({'-'},1,5);

        % Compute defocus levels
        for r = cone2RGCRatios   
            ratio = 2/r; % mm
            xThresh(r) = ratio; % convert to diopters
            labels{r} = sprintf('mRGC:Cone = %2.1f : 1.0',xThresh(r));
        end

        case 'current'
            
        polarAngles = [0, pi/2, pi, 3*pi/2];
        polarAngleLabels = {'nasal','superior','temporal','inferior'};
        
        % Get parameters to compute cone density levels
        whichEye         = 'left';
        cparams.cmFOV    =  1; % degrees
        
        % Convertion deg to m
        deg2m  = 0.3 * 0.001; % .3 deg per mm, .001 mm per meter
        
        for pa = 1:length(polarAngles)
            % Specify retinal location where stimulus is presented
            cparams.eccentricity = 4.5;                     % Visual angle of stimulus center, in deg
            cparams.polarAngle   = polarAngles(pa);   % Polar angle (radians): 0 is right, pi/2 is superior, pi is left, 3*pi/2 inferior
            
            % Compute x,y position in m of center of retinal patch from ecc and angle
            [x, y] = pol2cart(cparams.polarAngle, cparams.eccentricity);
            x = x * deg2m;  y = y * deg2m;
            cMosaic = coneMosaic('center', [x, y], 'whichEye', whichEye);
            
            % Set the field of view (degrees)
            cMosaic.setSizeToFOV(cparams.cmFOV);
            allDensity(pa) = eccen2density(cMosaic, 'deg');
            
            labels{pa} = sprintf('%s (%d cells/deg2)', polarAngleLabels{pa}, round(allDensity(pa)));
        end
            
        colors          = parula(length(labels)+1);
        colors          = colors(1:length(labels),:);
        xUnits          = linspace(min(expParams.contrastLevelsPC),max(expParams.contrastLevelsPC), 200);
        lineStyles      = repmat({'-'},1,5);
        xThresh           = allDensity;
        
end

