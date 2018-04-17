function expParams = loadExpParams(expName, saveParams)

% Define parameters for a specific experiment, and possibility to save in a mat file 
% so it can be reproduced later:

%       expParams = loadExpParams(expName, saveParams)

% INPUTS: 
%   ExpNames        : String defining which type of experiment parameters to use. 
%                       Choose from 'default',
%                                   'eyemov',
%                                   'eyemovEnhanced',
%                                   'coneDensity',
%                                   'defocus',
%                                   'ConeTypes',
%                                   'spatFreq'
%                     (default= 'default')
%   saveParams      : Boolean defining whether to save or not expParams in
%                      matfile (default = true)

%% Check input arguments
if isempty(expName) || ~exist('expName', 'var')
    expName = 'default'; 
end

if isempty(saveParams) || ~exist('saveParams', 'var')
    saveParams = true; 
end

expParams          = struct;
expParams.name     = expName;
expParams.nTrials  = 25;               % Number of trials per stimulus condition
expParams.verbose  = true;             % Print out images for debugging, or not


switch lower(expName)

    case 'default'  
        expParams.contrastLevels = [0:0.005:0.04, 0.05:0.01:0.1]; % Stimulus contrast levels (Michelson)
        expParams.eyemovement     = [1 1 0]';         % Type of eye movements: tremor, drift or ms?
        expParams.eccentricities  = 4.5;              % Eccentricity (deg);
        expParams.spatFreq        = 4;                % Spatial frequency (cycles/deg);
        expParams.polarAngle      = 0;                % Polar angle (radians): 0 is right, pi/2 is superior, pi is left, 3*pi/2 inferior
        expParams.defocusLevels   = 0;                % Value of first Zernike coeff (= defocus in units of ??) 
        expParams.verbose         = true;             % Print out images for debugging, or not

     case 'eyemov'
        expParams.contrastLevels  = [0:0.005:0.04, 0.05:0.01:0.1]; % Stimulus contrast levels (Michelson)
        expParams.eyemovement     = [0 0 0; 1 0 0; 0 1 0; 1 1 0]'; % None, only tremor, only drift, both tremor and drift
        expParams.eccentricities  = 4.5;              % Eccentricity (deg);
        expParams.spatFreq        = 4;                % Spatial frequency (cycles/deg);
        expParams.polarAngle      = 0;                % Polar angle (radians): 0 is right, pi/2 is superior, pi is left, 3*pi/2 inferior
        expParams.defocusLevels   = 0;                % Value of first Zernike coeff (= defocus in units of ??) 
        expParams.verbose         = true;             % Print out images for debugging, or not
     
      case 'conedensity'
        expParams.contrastLevels  = [0:0.005:0.04, 0.05:0.01:0.1];     % Stimulus contrast levels (Michelson)
        expParams.eyemovement     = [1 1 0]';         % Type of eye movements
        expParams.eccentricities  = [0 0.5 1 2 4.5 5 10 20 40]; % Eccentricity (deg);
        expParams.spatFreq        = 4;                % Spatial frequency (cycles/deg);
        expParams.polarAngle      = 0;                % Polar angle (radians): 0 is right, pi/2 is superior, pi is left, 3*pi/2 inferior
        expParams.defocusLevels   = 0;                % Value of first Zernike coeff (= defocus in units of ??) 
        expParams.verbose         = true;
        
      case 'conedensitydoeyemov'
        expParams.contrastLevels  = [0:0.005:0.04, 0.05:0.01:.1];     % Stimulus contrast levels (Michelson)
        expParams.eyemovement     = [0 0 0]';         % Type of eye movements
        expParams.eccentricities  = [0 0.5 1 2 4.5 5 10 20 40]; % Eccentricity (deg);
        expParams.spatFreq        = 4;                % Spatial frequency (cycles/deg);
        expParams.polarAngle      = 0;                % Polar angle (radians): 0 is right, pi/2 is superior, pi is left, 3*pi/2 inferior
        expParams.defocusLevels   = 0;                % Value of first Zernike coeff (= defocus in units of ??) 
        expParams.verbose         = true;
        
    case 'eyemovenhanced'
        expParams.contrastLevels  = [0:0.005:0.04, 0.05:0.01:.1]; % Stimulus contrast levels (Michelson)
        expParams.eyemovement     = [0 0 0; 2 0 0; 0 2 0; 2 2 0]';  % Which type of eye movments, emFlag will be turned into doubling amplitude or speed
        expParams.eccentricities  = 4.5;              % Eccentricity (deg);
        expParams.spatFreq        = 4;                % Spatial frequency (cycles/deg);
        expParams.polarAngle      = 0;                % Polar angle (radians): 0 is right, pi/2 is superior, pi is left, 3*pi/2 inferior
        expParams.defocusLevels   = 0;                % Value of first Zernike coeff (= defocus in units of ??)  
        expParams.verbose         = true;             % Print out images for debugging, or not
     
    case 'defocus'
        expParams.contrastLevels  = [0:0.005:0.04, 0.05:0.01:.1]; % Stimulus contrast levels (Michelson)
        expParams.eyemovement     = [1 1 0]';         % Which type of eye movments, emFlag will be turned into doubling amplitude or speed
        expParams.eccentricities  = 4.5;              % Eccentricity (deg);
        expParams.spatFreq        = 4;                % Spatial frequency (cycles/deg);
        expParams.polarAngle      = 0;                % Polar angle (radians): 0 is right, pi/2 is superior, pi is left, 3*pi/2 inferior
        expParams.defocusLevels   = [0 0.5 1 1.5 2];  % Value of first Zernike coeff (= defocus in units of ??)  
        expParams.verbose         = true;             % Print out images for debugging, or not
        
     case 'conetypes'
        expParams.contrastLevels  = [0:0.005:0.04, 0.05:0.01:.1]; % Stimulus contrast levels (Michelson)
        expParams.eyemovement     = [1 1 0]';         % Which type of eye movments, emFlag will be turned into doubling amplitude or speed
        expParams.eccentricities  = 4.5;              % Eccentricity (deg);
        expParams.spatFreq        = 4;                % Spatial frequency (cycles/deg);
        expParams.polarAngle      = 0;                % Polar angle (radians): 0 is right, pi/2 is superior, pi is left, 3*pi/2 inferior
        expParams.defocusLevels   = 0;                % Value of first Zernike coeff (= defocus in units of ??)  
        expParams.cparams.spatialDensity = [0 0.6 0.3 0.1; 0 1 0 0; 0 0 1 0; 0 0 0 1]; % Blank, L, M, S cone probabilities; Using default, or only one cone time at a time        

    case 'eccbasedconespacing'
        expParams.contrastLevels  = [0:0.005:0.04, 0.05:0.01:.1]; % Stimulus contrast levels (Michelson)
        expParams.eyemovement     = [1 1 0]';         % Which type of eye movments, emFlag will be turned into doubling amplitude or speed
        expParams.eccentricities  = [2.5, 10, 40];%[0, 2.5, 5, 10, 20, 40];              % Eccentricity (deg);
        expParams.spatFreq        = 4;                % Spatial frequency (cycles/deg);
        expParams.polarAngle      = 0;                % Polar angle (radians): 0 is right, pi/2 is superior, pi is left, 3*pi/2 inferior
        expParams.defocusLevels   = 0;                % Value of first Zernike coeff (= defocus in units of ??)         
        
        
        
        % NOT READY YET
%     case 'spatFreq' 
%         expParams.nTrials         = 25;               % Number of trials per stimulus condition
%         expParams.contrastLevels = [0:0.01:0.1 0.2]; % Stimulus contrast levels (Michelson)
%         expParams.eyemovement     = [1 1 0]';         % Which type of eye movments, emFlag will be turned into doubling amplitude or speed
%         expParams.eccentricities  = 4.5;              % Eccentricity (deg);
%         expParams.spatFreq        = [0.25, 0.4, 0.65, 1, 1.6, 2.6, 4, 8, 10, 16, 26]; % Spatial frequency (cycles/deg);
%         expParams.gausSDdeg       = [];
%         expParams.polarAngle     = 0;                % Polar angle (radians): 0 is right, pi/2 is superior, pi is left, 3*pi/2 inferior
%         expParams.defocuslevels   = 0;                % Value of first Zernike coeff (= defocus in units of ??)  
%         expParams.verbose         = true;             % Print out images for debugging, or not


end

% Save stimulus mat-file
if saveParams
    savePth = fullfile(ogRootPath, 'Stimulus');
    if ~exist(savePth, 'dir'); mkdir(savePth); end;
    save(fullfile(ogRootPath, 'Stimulus', sprintf('ogDefault_%s.mat',expName)),'expParams')
end

return