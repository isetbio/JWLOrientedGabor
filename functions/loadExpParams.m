function expParams = loadExpParams(expName, saveParams)
% Function to define parameters for a specific experiment, and possibility
% to save in a mat file so it can be reproduced later:
%
%       expParams = loadExpParams(expName, saveParams)
%
% INPUTS: 
%   ExpNames        : String defining which type of experiment parameters to use. 
%                       Choose from:
%                       'default'         - drift and ms, 4.5 deg eccen, 4 cpd, no defocus
%                       'eyemov',         - no eye movements, drift only, drift and ms, all at 4.5 deg eccen, 4 cpd, no defocus          
%                       'eyemovEnhanced', - NOT IMPLEMENTED ANYMORE due to new code implementation of eye movements
%                       'conedensity',    - sweep through eccentricities to get different cone densities in mosaic, all with drift and ms, 4 cpd, no defocus
%                       'defocus',        - sweep through defocus levels, all with drift and ms, at 4.5 deg eccen, 4 cpd
%                       'conetypes',      - create single cone type cone mosaics and default mosaic, all with drift and ms, at 4.5 deg eccen, 4 cpd, no defocus
%                       'conetypeslm90'   - create cone mosaics with LM cone type ratio 0.9:0.1 and vv, all with drift and ms, at 4.5 deg eccen, 4 cpd, no defocus
%                       'spatFreq'        - sweep through different spatial freq of stim, all with drift and ms, at 4.5 deg eccen, 4 cpd, no defocus
%                     (default= 'default')
%   saveParams      : Boolean defining whether to save or not expParams in
%                      matfile (default = true)
% OUTPUTS:
%   expParams       : struct with parameters to run a specific simulated
%                       experiment with the given 'expName'
%
% Example:
%   expParams = loadExpParams('default', false)
%
%
%% Check input arguments
if isempty(expName) || ~exist('expName', 'var')
    expName = 'default'; 
end

if isempty(saveParams) || ~exist('saveParams', 'var')
    saveParams = true; 
end

expParams                  = struct;
expParams.name             = expName;
expParams.nTrials          = 100;                                       % Number of trials per stimulus condition
expParams.verbose          = true;                                      % Print out images for debugging, or not
expParams.contrastLevels   = [0:0.005:0.04, 0.05:0.01:0.1];             % Stimulus contrast levels (Michelson)
expParams.contrastLevelsPC = [0:0.005:0.04, 0.05:0.01:.1, 0.2:0.1:1];   % Stimulus contrast levels when computing photo current (Michelson)
expParams.cparams.spatialDensity = [0 0.6 0.3 0.1];                     % Blank, L, M, S cone probabilities; Using default
expParams.cparams.noise    = 'random';

switch lower(expName)

    case 'default'  
        expParams.eyemovement         = [1 1]';           % Type of eye movements: tremor, drift or ms?
        expParams.eccentricities      = 4.5;              % Eccentricity (deg);
        expParams.spatFreq            = 4;                % Spatial frequency (cycles/deg);
        expParams.polarAngle          = 0;                % Polar angle (radians): 0 is right, pi/2 is superior, pi is left, 3*pi/2 inferior
        expParams.defocusLevels       = 0;                % Value of first Zernike coeff (= defocus in units of ??) 
        expParams.sparams.noStimPhase = false;            % Do NOT remove difference in phase of CCW and CW Gabor stimuli
     
    case 'eyemov'
        expParams.eyemovement         = [0 0; 1 0; 1 1]'; % None, only drift, both drift and MS
        expParams.eccentricities      = 4.5;              % Eccentricity (deg);
        expParams.spatFreq            = 4;                % Spatial frequency (cycles/deg);
        expParams.polarAngle          = 0;                % Polar angle (radians): 0 is right, pi/2 is superior, pi is left, 3*pi/2 inferior
        expParams.defocusLevels       = 0;                % Value of first Zernike coeff (= defocus in units of ??) 
        expParams.sparams.noStimPhase = false;            % Do NOT remove difference in phase of CCW and CW Gabor stimuli

      case 'conedensity'
        expParams.eyemovement         = [1 1]';                    % Which type of eye movements, drift and microsaccades
        expParams.eccentricities      = [0 0.5 1 2 4.5 5 10:5:40]; % Eccentricity (deg);
        expParams.spatFreq            = 4;                         % Spatial frequency (cycles/deg);
        expParams.polarAngle          = 0;                         % Polar angle (radians): 0 is right, pi/2 is superior, pi is left, 3*pi/2 inferior
        expParams.defocusLevels       = 0;                         % Value of first Zernike coeff (= defocus in units of ??)         
        expParams.sparams.noStimPhase = false;                     % Do NOT remove difference in phase of CCW and CW Gabor stimuli
 
      case 'conedensitynoeyemov'
        expParams.eyemovement         = [0 0]';                     % Type of eye movements
        expParams.eccentricities      = [0 0.5 1 2 4.5 5 10:5:40];  % Eccentricity (deg);
        expParams.spatFreq            = 4;                          % Spatial frequency (cycles/deg);
        expParams.polarAngle          = 0;                          % Polar angle (radians): 0 is right, pi/2 is superior, pi is left, 3*pi/2 inferior
        expParams.defocusLevels       = 0;                          % Value of first Zernike coeff (= defocus in units of ??) 
        expParams.sparams.noStimPhase = false;                      % Do NOT remove difference in phase of CCW and CW Gabor stimuli

     case 'eyemovenhanced'
         error('(%s): Experiment with enhanced eye movements is currently not implemented', mfilename);
%         expParams.eyemovement     = [0 0; 2 0; 2 2]';           % Which type of eye movments, emFlag will be turned into doubling amplitude or speed
%         expParams.eccentricities  = 4.5;                        % Eccentricity (deg);
%         expParams.spatFreq        = 4;                          % Spatial frequency (cycles/deg);
%         expParams.polarAngle      = 0;                          % Polar angle (radians): 0 is right, pi/2 is superior, pi is left, 3*pi/2 inferior
%         expParams.defocusLevels   = 0;                          % Value of first Zernike coeff (= defocus in units of ??)  
%      
    case 'defocus'
        expParams.eyemovement         = [1 1]';           % Which type of eye movements: drift and microsaccades
        expParams.eccentricities      = 4.5;              % Eccentricity (deg);
        expParams.spatFreq            = 4;                % Spatial frequency (cycles/deg);
        expParams.polarAngle          = 0;                % Polar angle (radians): 0 is right, pi/2 is superior, pi is left, 3*pi/2 inferior
        expParams.defocusLevels       = [0:0.25:2];       % Value of first Zernike coeff (= defocus in units of um)  
        expParams.sparams.noStimPhase = false;            % Do NOT remove difference in phase of CCW and CW Gabor stimuli

     case 'conetypes'
        expParams.eyemovement            = [1 1]';        % Which type of eye movements: drift and microsaccades
        expParams.eccentricities         = 4.5;           % Eccentricity (deg);
        expParams.spatFreq               = 4;             % Spatial frequency (cycles/deg);
        expParams.polarAngle             = 0;             % Polar angle (radians): 0 is right, pi/2 is superior, pi is left, 3*pi/2 inferior
        expParams.defocusLevels          = 0;             % Value of first Zernike coeff (= defocus in units of ??)  
        expParams.cparams.spatialDensity = [0 0.6 0.3 0.1; 0 1 0 0; 0 0 1 0; 0 0 0 1]; % Blank, L, M, S cone probabilities; Using default, or only one cone time at a time              
        expParams.sparams.noStimPhase    = false;         % Do NOT remove difference in phase of CCW and CW Gabor stimuli

     case 'conetypesmixed'
        expParams.eyemovement            = [1 1]';              % Which type of eye movements: drift and microsaccades
        expParams.eccentricities         = 4.5;                 % Eccentricity (deg);
        expParams.spatFreq               = 4;                   % Spatial frequency (cycles/deg);
        expParams.polarAngle             = 0;                   % Polar angle (radians): 0 is right, pi/2 is superior, pi is left, 3*pi/2 inferior
        expParams.defocusLevels          = 0;                   % Value of first Zernike coeff (= defocus in units of ??)  
        expParams.sparams.noStimPhase    = false;               % Do NOT remove difference in phase of CCW and CW Gabor stimuli
        expParams.cparams.spatialDensity = [ 0 1   0   0; ...   % L:M = 1:0,  % Blank, L, M, S cone probabilities                             
                                            0 0.9 0.1 0; ...    % L:M = 0.9:0.1 
                                            0 0.8 0.2 0; ...    % L:M = 0.8:0.2
                                            0 0.7 0.3 0; ...    % L:M = 0.7:0.3
                                            0 0.6 0.4 0; ...    % L:M = 0.6:0.4
                                            0 0.5 0.5 0; ...    % L:M = 0.5:0.5
                                            0 0.4 0.6 0; ...    % L:M = 0.4:0.6
                                            0 0.3 0.7 0; ...    % L:M = 0.3:0.7
                                            0 0.2 0.8 0; ...    % L:M = 0.2:0.8
                                            0 0.1 0.9 0; ...    % L:M = 0.1:0.9
                                            0 0   1   0];       % L:M = 0:1  
    case 'idealobserver'        
        expParams.eyemovement         = [0 0]';           % No eye movements
        expParams.eccentricities      = 4.5;              % Eccentricity (deg);
        expParams.spatFreq            = 4;                % Spatial frequency (cycles/deg);
        expParams.polarAngle          = 0;                % Polar angle (radians): 0 is right, pi/2 is superior, pi is left, 3*pi/2 inferior
        expParams.defocusLevels       = 0;                % Value of first Zernike coeff (= defocus in units of ??) 
        expParams.cparams.noise       = 'none';           % no photon noise
        expParams.sparams.noStimPhase = true;             % Remove difference in phase of CCW and CW Gabor stimuli
        expParams.contrastLevels      = [0:0.0001:0.001, 0.002:0.001:0.01, 0.02:0.01:0.1]; % Stimulus contrast levels (Michelson)

    case 'defaultnophaseshift'        
        expParams.eyemovement         = [0 0]';           % No eye movements
        expParams.eccentricities  	  = 4.5;              % Eccentricity (deg);
        expParams.spatFreq            = 4;                % Spatial frequency (cycles/deg);
        expParams.polarAngle          = 0;                % Polar angle (radians): 0 is right, pi/2 is superior, pi is left, 3*pi/2 inferior
        expParams.defocusLevels       = 0;                % Value of first Zernike coeff (= defocus in units of ??) 
        expParams.cparams.noise       = 'random';         % poisson photon noise
        expParams.sparams.noStimPhase = true;             % Remove difference in phase of CCW and CW Gabor stimuli
        expParams.contrastLevels      = [0:0.0001:0.001, 0.002:0.001:0.01, 0.02:0.01:0.1]; % Stimulus contrast levels (Michelson)
        expParams.contrastLevelsPC    = [0:0.0001:0.001, 0.002:0.001:0.01, 0.02:0.01:0.1]; % Stimulus contrast levels (Michelson)
        
    case 'eyemovnophaseshift'        
        expParams.eyemovement         = [0 0; 1 0; 1 1]'; % No eye movements, drift only, drift and MS
        expParams.eccentricities      = 4.5;              % Eccentricity (deg);
        expParams.spatFreq            = 4;                % Spatial frequency (cycles/deg);
        expParams.polarAngle          = 0;                % Polar angle (radians): 0 is right, pi/2 is superior, pi is left, 3*pi/2 inferior
        expParams.defocusLevels       = 0;                % Value of first Zernike coeff (= defocus in units of ??) 
        expParams.cparams.noise       = 'random';         % no photon noise
        expParams.sparams.noStimPhase = true;             % Remove difference in phase of CCW and CW Gabor stimuli
        expParams.contrastLevels      = [0:0.0001:0.001, 0.002:0.001:0.01, 0.02:0.01:0.1];  % Stimulus contrast levels (Michelson)

    case 'conetypeseccen'   
        expParams.eyemovement            = [1 1]';                      % Which type of eye movements: drift and microsaccades
        expParams.eccentricities         = [0 0.5 1 2 4.5 5 10:5:40];   % Eccentricity (deg);
        expParams.spatFreq               = 4;                           % Spatial frequency (cycles/deg);
        expParams.polarAngle             = 0;                           % Polar angle (radians): 0 is right, pi/2 is superior, pi is left, 3*pi/2 inferior
        expParams.defocusLevels          = 0;                           % Value of first Zernike coeff (= defocus in units of ??)  
        expParams.cparams.spatialDensity = [0 1 0 0; 0 0 1 0; 0 0 0 1]; % Blank, L, M, S cone probabilities; Using default, or only one cone time at a time  
        expParams.sparams.noStimPhase    = false;                       % Do NOT remove difference in phase of CCW and CW Gabor stimuli

    case 'spatfreq' 
        expParams.eyemovement         = [1 1]';           % Which type of eye movements: drift and microsaccades
        expParams.eccentricities      = 4.5;              % Eccentricity (deg);
        expParams.spatFreq            = [0.25, 0.4, 0.65, 1, 1.6, 2.6, 4, 8, 10, 16, 26]; % Spatial frequency (cycles/deg);
        expParams.polarAngle          = 0;                % Polar angle (radians): 0 is right, pi/2 is superior, pi is left, 3*pi/2 inferior
        expParams.defocusLevels       = 0;                % Value of first Zernike coeff (= defocus in units of ??)  
        expParams.sparams.noStimPhase = false;            % Do NOT remove difference in phase of CCW and CW Gabor stimuli
end

% Save stimulus mat-file
if saveParams
    savePth = fullfile(ogRootPath, 'Stimulus');
    if ~exist(savePth, 'dir'); mkdir(savePth); end
    save(fullfile(ogRootPath, 'Stimulus', sprintf('ogDefault_%s.mat',expName)),'expParams')
end

% Get toolbox versions used
out = ver;

if saveParams
    save(fullfile(ogRootPath, 'toolboxVersions.mat'),'out')
end

return