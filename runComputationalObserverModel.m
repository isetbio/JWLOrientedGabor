function runComputationalObserverModel(expName, varargin)
%% runComputationalObserverModel(expName, [saveFolder], [], [seed], 1, [currentFlag], 0)
%
% ---------------------------------------------------------
% ---------------------- Description ----------------------
% ---------------------------------------------------------
% Function based on the old script s_ogRGC.m
% This function computes multiple stages of the front-end of the visual system
% for a particular visual scene or psychophysical experiment.
%
% In this case we construct a scene with 2 Gabor stimuli (oriented clockwise
% or counter-clockwise), that goes through the optics and gets sampled by
% a small cone mosaic patch. The cone absorptions are then used to simulate
% a 2-AFC orientation discrimination task with a linear classifier.
%
% Long-term goals:
%   Model cone current, bipolar, RGC, cortical, and behavioral (computational
% observer) responses for the experiment on measuring orientation
% discrimination thresholds of an achromatic, peripheral Gabor. Account for
% variation in optics, cone density, and RGC density/RF size as a function
% of visual field position.
%
% ---------------------------------------------------------
% ---------------------- Script outline -------------------
% ---------------------------------------------------------
% Stages:
%       * Check inputs and define experimental manipulation
%       * Create default scene and achromatic Gabor patches (OIS, optical image sequence).
%       * Create cone mosaic
%       * Create optics
%       * Create eyemovements
%       * Update scene and stimuli by setting contrast and spatial frequency
%       * Compute cone absorptions
%       * Classify absorptions
%
% ---------------------------------------------------------
% ---------------------- Inputs ---------------------------
% ---------------------------------------------------------
% expName       : (string) name of experimental condition to simulate
%                 following options are available:
%                   - 'default' (4.5 deg eccen, typical optics, eye movements, and cone mosaic)
%                   - 'eyemov'  (vary eye movements: no eye movement, drift or drift and microsaccades)
%                   - 'defocus' (vary defocus levels in optics: 0-2 diopters)
%                   - 'conedensity' (vary cone density levels in mosaic: from density at fovea up to 30 deg eccentricity)
%                   - 'conetypes'   (vary homogeneity of cone types in mosaic: L-, M- or S-cone only vs default LMS)
%                   - 'conetypesmixed' (vary mixture of LM cone ratio in mosaic: from 100% L-cones to 100% M-cones)
%                   - 'idealobserver' simulate experiment without photon noise, eye movements or phase shifts in Gabor stimuli
%                 for more info and options, see loadExpParams.m
% [saveFolder]   : (string) name of folder to save simulated cone absorptions (default =[])
% [seed]         : (int) integer to reset random number generator to reproduce results (default = 1)
% [currentFlag]  : (bool) flag to compute cone current in addition to cone absorptions (default = false)
%
%
% The results of the computational observer model are written to file and
% can be then be used to reproduce several of the published figures from the
% paper the paper:
%
% Eline R. Kupers, Marisa Carrasco, Jonathan Winawer. (XXXX) Modeling
% visual performance differences 'around' the visual field: A computational
% observer approach.
%
% ---------------------------------------------------------
% ---------------------- Examples -------------------------
% ---------------------------------------------------------
%
% Example 1: Compute and classify cone absorptions for CW and CCW Gabor 
% stimuli, with a cone mosaic at 4.5 deg eccentricity, typical human optics 
%   runComputationalObserverModel('default')
% Example 2: Quantify the effect of defocus in human optics
%   runComputationalObserverModel('defocus')
% Example 3: Quantify the effect of cone density, with fixed rng seed
%   runComputationalObserverModel('conedensity', 'saveFolder', 'rng1', 'seed', 1)
% Example 4: Quantify the effect of cone types, with reset rng seed
%   runComputationalObserverModel('conetypes', 'saveFolder', 'rng1', 'seed', 'shuffle')
%
%
% EK/JW/ NYU ISETBIO Team, Copyright 2018

%% 0. Check inputs and define experimental parameters
p = inputParser;
p.KeepUnmatched = true;
p.addRequired('expName', @ischar);
p.addParameter('saveFolder', [], @ischar);
p.addParameter('seed', 1, @(x) (isstring(x) | isscalar(x)));
p.addParameter('currentFlag', false, @islogical);
p.parse(expName, varargin{:});

% Check and create folder to save absorption data
currDate = datestr(datetime,'yyyymmdd_HHMMSS');

if ~isempty(p.Results.saveFolder)
    saveFolder = fullfile(ogRootPath, 'data', expName, p.Results.saveFolder);
    saveFolderClassification = fullfile(ogRootPath, 'data', 'classification', expName, p.Results.saveFolder);
else
    saveFolder = fullfile(ogRootPath, 'data', expName, currDate);
    saveFolderClassification = fullfile(ogRootPath, 'data', 'classification', expName, currDate);
end

% Create folders if they don't exist
if ~exist('saveFolder', 'dir'); mkdir(saveFolder); end
if ~exist('saveFolderClassification', 'dir'); mkdir(saveFolderClassification); end

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
    expParams.currentFlag = p.Results.currentFlag;
else
    theseContrasts = 1.0; %expParams.contrastLevels;
    expParams.currentFlag = p.Results.currentFlag;
end

% Check if ideal observer is requested
% expParams.idealObserver = p.Results.idealObserver;

if expParams.verbose
    fH = figure(99); clf; hold all;
end

%% ------------------- DEFAULT SCENE and STIMULI -------------------
% Get scene radiance and default stimulus (OG, i.e. the Oriented Gabors)
if expParams.verbose; fprintf('(%s): Creating scene.\n', mfilename); end
[OG, scenes, sparams] = getSceneAndStimuli; %#ok<ASGLU>

% We make sure that the number of time points in the eye movement sequence
% matches the number of time points in the optical image sequence
tSamples         = OG(1).length;

% Loop over conditions, generating cone absorptions for each condition
% 1. Eccentricity (hack to get different cone densities)
% 2. Optics
% 3. Eye movements
% 4. Contrasts
% 5. Spatial frequencies
% 6. OIS (Oriented Gabors) when computing absorptions

for eccen = expParams.eccentricities  % loop over eccentricity (aka cone density) levels
    
    for lmsIdx = 1:size(expParams.cparams.spatialDensity,1)
        lmsRatio = expParams.cparams.spatialDensity(lmsIdx,:);
        
        %% ------------------- MOSAIC -------------------
        if expParams.verbose; fprintf('(%s): Creating mosaic.\n',mfilename); end
        [cMosaic, cparams] = getConeMosaic(eccen, expParams, sparams, lmsRatio);
        
        % Integration time can be defined independently from OIS time step.
        % Prefered to be 5 ms or lower (1 or 2 ms preferred)
        cMosaic.integrationTime = 0.002; % ms
        
        % Change cone spacing based on eccentricity
        if strcmp(expName,'conedensity')
            propCovered = getBanks1991ConeCoverage(eccen); % proportion
            
            cMosaic.pigment.pdWidth  = cMosaic.pigment.width*propCovered; % meters
            cMosaic.pigment.pdHeight = cMosaic.pigment.height*propCovered; % meters
        end
        
        for defocus = expParams.defocusLevels  % loop over defocus conditions
            
            %% ------------------- OPTICS -------------------
            if expParams.verbose; fprintf('(%s): Adding optics.\n', mfilename); end
            sparams.oi = oiDefocus(defocus, expParams.verbose); % input is Zernicke defocus coeff
            
            
            %% ------------------- EYE MOVEMENTS -------------------
            for emIdx = 1:size(expParams.eyemovement,2) % loop over eye movement conditions
                
                accuracy = NaN(size(theseContrasts));
                
                for c = 1 %theseContrasts % loop over contrasts
                    
                    for sf = expParams.spatFreq % loop over spatial frequencies
                        
                        
                        %% ------------------- UPDATE SCENE and STIMULI (Contrast and SF) -------------------
                        if expParams.verbose; fprintf('(%s): Computing absorptions for stimulus contrast %1.4f, polar angle %d, eccen %1.2f, LMS ratio %1.1f:%1.1f:%1.1f\n', mfilename, c, expParams.polarAngle, eccen, lmsRatio(2),lmsRatio(3),lmsRatio(4)); end
                        fname = sprintf('OGconeOutputs_contrast%1.4f_pa%d_eye%d%d_eccen%1.2f_defocus%1.2f_noise-%s_sf%1.2f_lms-%1.1f%1.1f%1.1f.mat',...
                            c,expParams.polarAngle,expParams.eyemovement(1,emIdx),expParams.eyemovement(2,emIdx), eccen, defocus, cMosaic.noiseFlag, sf, lmsRatio(2),lmsRatio(3),lmsRatio(4));
                        if expParams.verbose;  fprintf('(%s): File will be saved as %s\n', mfilename, fname); end
                        
                        % Update the stimulus contrast & spatial frequency
                        if expParams.verbose; fprintf('(%s): Recomputing scene for current sf: %1.2f and c: %1.4f..\n', mfilename, sf, c); end
                        
                        
                        sparams.gabor.contrast  = c;  % Michelson, range = [0 1]
                        sparams.freqCPD         = sf; % Cycles/degree
                        if expParams.sparams.noStimPhase; sparams.noStimPhase = expParams.sparams.noStimPhase; end
                        [OG,scenes,tseries]     = ogStimuli(sparams);
                        
                        
                        %% ------------------- COMPUTE ABSORPTIONS  -------------------
                        if expParams.verbose;  fprintf('(%s): Compute absorptions.\n', mfilename); end
                        
                        % Allocate space for absorptions
                        absorptions = zeros(expParams.nTrials,cMosaic.rows,cMosaic.cols, tSamples, length(OG));
                        current     = absorptions;
                        
                        for s = 1:length(OG) % loop over OIS' (CW Gabor x 2 phases and CCW Gabor x 2 phases)
                            
                            
                            %% ------------------- EYE MOVEMENTS -------------------
                            if expParams.verbose; fprintf('(%s): Adding eyemovements.\n', mfilename); end
                            
                            [emPaths, cMosaic] = getEyemovements(OG, cMosaic, expParams, sparams, emIdx, expParams.seed+s);
                            
                            % Add emPaths (which are in terms of cones shifted) to cMosaic struct
                            cMosaic.emPositions = emPaths;
                            
                            if expParams.currentFlag
                                [absorptions(:,:,:,:,s), current(:,:,:,:,s), interpFilters, meanCur] = cMosaic.compute(OG(s), 'currentFlag', expParams.currentFlag, ...
                                    'emPaths', emPaths, 'seed', expParams.seed+s); % add new +1 to seed for new set of eye movement trials
                            else
                                absorptions(:,:,:,:,s) = cMosaic.compute(OG(s), 'currentFlag', expParams.currentFlag, ...
                                    'emPaths', emPaths, 'seed', expParams.seed+s); % add new +1 to seed for new set of eye movement trials
                            end
                        end
                        
                        % Save cone absorption data
                        if expParams.verbose; fprintf('(%s): Saving cone absorption data.\n', mfilename); end
                        parsave(fullfile(saveFolder, fname), ...
                            'absorptions', absorptions, ...
                            'sparams', sparams, ...
                            'cparams', cparams, ...
                            'expParams', expParams, ...
                            'emPaths', emPaths);
                    
                        % If cone current was requested, also save this array
                        if expParams.currentFlag
                            fprintf('(%s): Saving cone current data..\n', mfilename);
                            parsave(fullfile(saveFolder, ['current_' fname]), ...
                                'current', current, ...
                                'interpFilters', interpFilters, ...
                                'meanCur', meanCur, ...
                                'sparams', sparams, ...
                                'cparams', cparams, ...
                                'expParams', expParams, ...
                                'emPaths', emPaths);
                        end
                        
                        %% ------------------- Classify absorptions  -------------------
                        fnameClassify = sprintf(...
                            'Classify_coneOutputs_contrast%1.4f_pa%d_eye%s_eccen%1.2f_defocus%1.2f_noise-%s_sf%1.2f_lms-%1.1f%1.1f%1.1f',...
                            c, expParams.polarAngle,sprintf('%i',expParams.eyemovement(:,emIdx)), eccen, defocus, cMosaic.noiseFlag, sf, lmsRatio(2),lmsRatio(3),lmsRatio(4));
                        
                        if expParams.verbose
                            fprintf('(%s): Classify cone absorption data..\n', mfilename);
                            fprintf('(%s): File will be saved as %s\n', mfilename, fnameClassify);
                        end
                        
                        if expParams.currentFlag
                            accuracy(c==theseContrasts) = getClassifierAccuracy(current);
                            fnameClassify = ['current_' fnameClassify]; %#ok<AGROW>
                        else
                            accuracy(c==theseContrasts) = getClassifierAccuracy(absorptions); % truncate time samples (only include stimulus on period)
                        end
                        
                        if expParams.verbose; fprintf('(%s): Classifier accuracy for stim contrast %1.4f is %3.2f..\n', mfilename, c, accuracy(c==theseContrasts)); end
                        
                    end % sf
                end % contrast
                
                % Save
                parsave(fullfile(saveFolderClassification, sprintf('%s.mat', fnameClassify)),'accuracy',accuracy, 'expParams', expParams);
                
                % Visualize if verbose
                if expParams.verbose; set(0, 'CurrentFigure', fH); plot(theseContrasts, accuracy,'o-', 'LineWidth',2); drawnow; end
                
            end % eyemovements
        end % defocus
    end % cone types
end % eccentricities

return