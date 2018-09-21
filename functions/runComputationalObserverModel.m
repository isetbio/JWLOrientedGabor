function runComputationalObserverModel(expName, saveFolder)
%% runComputationalObserverModel(expName, saveFolder, seed)

% ---------------------- Description ----------------------

% Function based on the old script s_ogRGC.m
% This function computes multiple stages of the front-end of the visual system
% for a particular visual scene or psychophysical experiment.

% In this case we construct a scene with 2 Gabor stimuli (oriented clockwise
% or counter-clockwise), that goes through the optics and gets sampled by
% a small cone mosaic patch. The cone absorptions are then used to simulate
% a 2-AFC orientation discrimination task with a linear classifier.

% Long-term goals:
%   Model cone current, bipolar, RGC, cortical, and behavioral (computational
% observer) responses for the experiment on measuring orientation
% discrimination thresholds of an achromatic, peripheral Gabor. Account for
% variation in optics, cone density, and RGC density/RF size as a function
% of visual field position.

% ---------------------- Script outline ----------------------
%
%   Stages:
%       * Check inputs and define experimental manipulation
%       * Create default scene and achromatic Gabor patches (OIS, optical image sequence).
%       * Create cone mosaic
%       * Create optics 
%       * Create eyemovements
%       * Update scene and stimuli by setting contrast and spatial frequency
%       * Compute cone absorptions
%       * Classify absorptions
%
% EK/JW/ NYU ISETBIO Team, Copyright 2018


%% 0. Check inputs and define experimental parameters
p = inputParser;
p.addRequired('expName', @isstring);
p.addRequired('saveFolder', @isstring);
p.addParameter('seed', 1, @(x) (isstring(x) | isscalar(x)));
p.addParameter('currentFlag', false, islogical);
p.parse(expName, saveFolder, varargin{:});

% Specify experiment parameters
expParams = loadExpParams(expName, false);   % (false argument is for not saving params in separate matfile)

% Define deg2m converter
expParams.deg2m   = 0.3 * 0.001;     % (default in isetbio)

% Set random number generator seed for reproducibility
rng(p.Results.seed);
expParams.seed = p.Results.seed;

% Check if current is requested, then we want to add more contrast levels
if p.Results.currentFlag
    theseContrasts = expParams.contrastLevelsPC;
    expParams.currentFlag = currentFlag;
else
    theseContrasts = expParams.contrastLevels;
    expParams.currentFlag = currentFlag;
end

%% ------------------- DEFAULT SCENE and STIMULI -------------------
% Get scene radiance and default stimulus (OG, i.e. the Oriented Gabors)
if expParams.verbose; fprintf('(%s): Creating scene.\n', mfilename); end
[OG, scenes, sparams] = getSceneAndStimuli;

% We make sure that the number of time points in the eye movement sequence
% matches the number of time points in the optical image sequence
tSamples         = OG(1).length;

% Loop over conditions, generating cone absorptions for each condition

for eccen = expParams.eccentricities
    
    %% ------------------- MOSAIC -------------------
    if expParams.verbose; fprintf('(%s): Creating mosaic.\n',mfilename); end
    [cMosaic, cparams] = getConeMosaic(eccen, expParams);
    
    % Integration time can be defined independently from OIS time step.
    % Prefered to be 5 ms or lower (1 or 2 ms preferred)
    cMosaic.integrationTime = 0.002; % ms
    
    % Change cone spacing based on eccentricity
    if strcmp(expName,'eccbasedcoverage')
        propCovered = getBanks1991ConeCoverage(eccen); % proportion
        
        cMosaic.pigment.pdWidth  = cMosaic.pigment.width*propCovered; % meters
        cMosaic.pigment.pdHeight = cMosaic.pigment.height*propCovered; % meters
    end
    
    for defocus = expParams.defocusLevels
        
        %% ------------------- OPTICS -------------------
        if expParams.verbose; fprintf('(%s): Setting defocus to %s (microns)\n', mfilename, defocus); end
        sparams.oi = oiDefocus(defocus); % input is Zernicke defocus coeff
        
        
        %% ------------------- EYE MOVEMENTS -------------------
        for emIdx = 1:size(expParams.eyemovement,2)
            
            if expParams.verbose; fprintf('(%s): Defining eyemovements as %s (=drift, ms)..\n', mfilename, mat2str(expParams.eyemovement(:,emIdx))); end
            
            % Get the eyemovements
            [emPaths, cMosaic] = getEyemovements(OG, cMosaic, expParams, sparams);
            
            % Add emPaths (which are in terms of cones shifted) to cMosaic struct
            cMosaic.emPositions = emPaths;
            
            
            for c = theseContrasts
                
                for sf = expParams.spatFreq
                    
                    
                    %% ------------------- UPDATE SCENE and STIMULI  -------------------
                    if expParams.verbose; fprintf('(%s): Computing absorptions for stimulus contrast %4.3f, polar angle %d, eccen %1.2f\n', mfilename, c, expParams.polarAngle, eccen); end                   
                    fname = sprintf('OGconeOutputs_contrast%1.3f_pa%d_eye%d%d_eccen%1.2f_defocus%1.2f_noise-%s_sf%1.2f.mat',...
                        c,expParams.polarAngle,expParams.eyemovement(1,emIdx),expParams.eyemovement(2,emIdx), eccen, defocus, cMosaic.noiseFlag, sf);                    
                    if expParams.verbose;  fprintf('(%s): File will be saved as %s\n', mfilename, fname); end
                    
                    % Update the stimulus contrast & spatial frequency
                    if expParams.verbose; fprintf('(%s): Recomputing scene for current sf: %f and c: %1.2f..\n', mfilename, sf, c); end

                    sparams.gabor.contrast  = c;  % Michelson, [0 1]
                    sparams.freqCPD         = sf; % Cycles/degree        
                    [OG,scenes,tseries] = ogStimuli(sparams);
                    
                    
                    %% ------------------- COMPUTE ABSORPTIONS  -------------------
                    
                    % Compute absorptions for multiple trials
                    absorptions = zeros(expParams.nTrials,cMosaic.rows,cMosaic.cols, tSamples, length(OG));
                    current     = absorptions;
                    
                    if expParams.verbose;  fprintf('(%s): Compute absorptions.\n', fname); end

                    for s = 1:length(OG)
                        if expParams.currentFlag
                            [absorptions(:,:,:,:,s), current(:,:,:,:,s), interpFilters, meanCur] = cMosaic.compute(OG(s), 'currentFlag', expParams.currentFlag, ...
                                'emPaths', emPaths, 'seed', seed);
                        else
                            absorptions(:,:,:,:,s) = cMosaic.compute(OG(s), 'currentFlag', false, ...
                                'emPaths', emPaths, 'seed', seed);
                        end
                    end
                    
                    % Save data
                    if expParams.verbose; fprintf('(%s): Saving data..\n', mfilename); end
                    savePth = fullfile(ogRootPath, 'data', expName, subFolderName); if ~exist(savePth,'dir'); mkdir(savePth); end;
                    parsave(fullfile(savePth, fname), 'absorptions', absorptions, 'sparams', sparams, 'cparams', cparams, 'expParams', expParams, 'emPaths', emPaths);
                    if expParams.currentFlag; parsave(fullfile(savePth, ['current_' fname]), 'current', current, 'sparams', sparams, 'cparams', cparams, 'expParams', expParams, 'emPaths', emPaths); end
                end % sf
            end % contrast
        end % defocus
    end % eyemovements
end % eccentricities


%% ------------------- Classify absorptions  ------------------- 


return
