%% Oriented Gabor Discrimination
% This script supercedes t_orientedGaborDiscrimination.m, which no longer runs,
% due to changes in isetbio.
%
% Long-term goal:
%   Model cone, bipolar, RGC, cortical, and behavioral (computational
% observer) responses for the experiment on measuring orientation
% discrimination thresholds of an achromatic, peripheral Gabor. Account for
% variation in optics, cone density, and RGC density/RF size as a function
% of visual field position.
%
% So far:
%   Model cone, bipolar, RGC responses

% A temporary offsheet of s_ogRGC for testing purposes. Merge or delete
% when finished testing.

%% ---------------------- Script outline ----------------------
%
%   Build structures:
%
%       * Build scene as an achromatic Gabor patch with imageHarmonic.
%       * Build OIS (optical image sequence)
%       * Cone mosaic
%       * Bipolars
%       * RGCs
%       * Cortex
%       * Classifier

%% ---------------------- Questions ----------------------
%
%       * Spatial and temporal pooling for cortex? (Contact Nicholas C!)
%       * TODO: Incorporate off-axis optical data from Pablo Artal
%               Example: Figure 5 from Jaeken and Artal, 2012
%           It appears that:
%               RGC density is cell-type dependent (e.g. Offs denser than
%               Ons, which is good). The functions of position come from
%               EJ. What abou Curcio? Watson (2014)? EJ, for example, has
%               inferior and superior matched. But Watson (using Curcio)
%               does not.
%       * How is RGC rf size determined as a function of eccentricity?
%
%% ---------------------- TODO ---------------------------
%       * How to implement position-specific optics?
%       * Check position specific cone and RGC densities
%           For LE at 6º ecc and 2º FOV and polar angle of
%           -   0 deg (nasal),    the mosaic is 71x71
%           -  90 deg (superior), the mosaic is 64x64
%           - 270 deg (inferior), the mosaic is 65x65
%       * How to validate bipolar responses (including for different class
%           types?)
%
%
%% ---------------------- Experiment ----------------------
%
%       *   Loop over 4 spatial positions (left, right, lower, upper)
%       *   For each position, loop over 2 stimulus orientations and n
%               trials per orientation
%       *   Within each trial, add eye movements to sensor (or OIS?)
%       *   Build the outer segment object with linear filters and compute its
%              response with the sensor structure.
%       *   Build the bipolar layer and compute responses (which cell types?)
%       *   Build the RGC layer and compute responses (which cell types?)
%       *   Compute cortical responses from spatial / temporal responses
%       *   Compute linear discriminant of stimulus orientation using
%              outputs at various stages
%
% EK/JW/ NYU ISETBIO Team, Copyright 2017


%% Specify experiment parameters

% Load experiment parameters
expName = 'eccbasedcoverage';
expParams = loadExpParams(expName, false);

% Temporal properties of one trial
tStep             = 0.002;                % Time step for optical image sequence (seconds)
sparams.tsamples  = (0:tStep:0.054);      % seconds
sparams.timesd    = 1000.00;              % sd of temporal Gaussian window

% Scene field of view
sparams.sceneFOV  = 2;   % scene field of view in degrees (diameter)
sparams.freqCPD   = 4;   % Gabor spatial frequency (cpd)
sparams.gausSDdeg = sparams.sceneFOV/4; % Gabor SD in degrees of visual angle

% Unit converters
deg2fov = 1/sparams.sceneFOV;
fov2deg = sparams.sceneFOV;
deg2m   = 0.3 * 0.001;          % (we first choose 3 deg per mm, .001 mm per meter, but now adjusted to the default in isetbio)

% Gabor parameters
sparams.gabor           = harmonicP;                   % Standard Gabor
sparams.gabor.ang       = (pi/180)* 15;                % Gabor orientation (radians) - question: what is 0??
sparams.gabor.freq      = fov2deg*sparams.freqCPD;     % Spatial frequency (cycles/FOV)
sparams.gabor.contrast  = 1;                           % Presumably michelson, [0 1]
sparams.gabor.GaborFlag = sparams.gausSDdeg*deg2fov;   % Gaussian window

% Make a default Optical Image Sequence. We need some of the
% parameters from this. It will be overwritten as we change the
% contrast and/or optics.
OG = ogStimuli(sparams);

% We make sure that the number of time points in the eye movement sequence
% matches the number of time points in the optical image sequence
tSamples         = OG(1).length;

% Set if we want to compute cone current from cone absorptions
currentFlag = true;

%% Loop over conditions, generating cone absorptions for each condition

% Eccentricity loop
for eccen = expParams.eccentricities
    
    % ----- CONE MOSAIC -----------------------------------------
    % Make CONE MOSAIC for a given eccentricity and polar angle
    whichEye = 'left';
    
    % Specify retinal location where stimulus is presented
    cparams.eccentricity      = eccen;             % Visual angle of stimulus center, in deg
    cparams.polarAngle        = expParams.polarAngle; % Polar angle (radians): 0 is right, pi/2 is superior, pi is left, 3*pi/2 inferior
    
    % Cone mosaic field of view in degrees
    cparams.cmFOV        = sparams.sceneFOV; % degrees
    
    % Compute x,y position in m of center of retinal patch from ecc and angle
    [x, y] = pol2cart(cparams.polarAngle, cparams.eccentricity);
    x = x * deg2m;  y = y * deg2m;
    
    % Create coneMosaic for particular location and eye
    cMosaic = coneMosaic('center', [x, y], 'whichEye', whichEye);
    
    % Set the field of view (degrees)
    cMosaic.setSizeToFOV(cparams.cmFOV);
    
    % Add photon noise
    cMosaic.noiseFlag = 'random'; % 'random' 'frozen' 'none'
    
    % CURRENT: Set outer segment to be computed with linear filters
    cMosaic.os = osLinear;
    
    for defocus = expParams.defocusLevels
        
        % ---- Add optics blur or defocus if requested
        sparams.oi = oiDefocus(defocus); % input is Zernicke defocus coeff
        
        % Change cone spacing based on eccentricity
        if strcmp(expName,'eccbasedcoverage')
            propCovered = getBanks1991ConeCoverage(eccen);

            cMosaic.pigment.pdWidth  = cMosaic.pigment.width*propCovered;
            cMosaic.pigment.pdHeight = cMosaic.pigment.height*propCovered;
        end
        
        % ----- EYE MOVEMENTS -------------------------------------
        % Make EYE MOVEMENTS for a given cone mosaic
        
        % Integration time can be defined independently from OIS time step.
        % Prefered to be 5 ms or lower (1 or 2 ms preferred)
        cMosaic.integrationTime = 0.002; %OG(1).timeStep;
        
%         for emIdx = 1:size(expParams.eyemovement,2)
%             if expParams.verbose; fprintf('Defining eyemovements as %s (=tremor, drift, ms)..\n', mat2str(expParams.eyemovement(:,emIdx))); end
            % Include tremor, drift, microsaccades?
%             cparams.em        = emCreate;
%             cparams.em.emFlag = expParams.eyemovement(:,emIdx);
            
%             if cparams.em.emFlag(1) == 2
%                 cparams.em.tremor.amplitude = cparams.em.tremor.amplitude * cparams.em.emFlag(1);
%                 cparams.em.emFlag(1) = [cparams.em.emFlag(1)./cparams.em.emFlag(1)]'; % set emFlag back to 1
%                 warning('(s_ogRGC: tremor amplitude enhanced)')
%             end
%             if cparams.em.emFlag(2) == 2
%                 cparams.em.drift.speed = cparams.em.drift.speed * cparams.em.emFlag(2);
%                 cparams.em.drift.speedSD = cparams.em.drift.speedSD * cparams.em.emFlag(2);
%                 cparams.em.emFlag(2) = [cparams.em.emFlag(2)./cparams.em.emFlag(2)]'; % set emFlag back to 1
%                 warning('(s_ogRGC: drift speed enhanced)')
%             end
            
            % Generate the eye movement paths in units of cone samples. Set the
            % time sample to 2x the actual time, so that the eye does not
            % always start at 0,0. We will then clip after generating the eye
            % movements.
            maxEyeMovementsNum = OG(1).maxEyeMovementsNumGivenIntegrationTime(cMosaic.integrationTime);
            emPaths = cMosaic.emGenSequence(maxEyeMovementsNum*2, 'nTrials', expParams.nTrials); % path is in terms of cones shifted
            emPaths = emPaths(:, end-maxEyeMovementsNum+1:end,:);
            cMosaic.emPositions = squeeze(emPaths(1,:,:));
            
            if expParams.verbose
                %plot eye movements
                figure,
                subplot(211)
                plot(sparams.tsamples, emPaths(:,:,1)')
                
                subplot(212)
                plot(sparams.tsamples, emPaths(:,:,2)')
            end
            
            
            
            % Loop over contrasts and defocus
            for c = expParams.contrastLevels
                
                for sf = expParams.spatFreq
                    
                    if expParams.verbose; fprintf('Computing absorptions for stimulus contrast %4.3f, polar angle %d, eccen %1.2f\n', c, expParams.polarAngle, eccen); end
                    fname = sprintf('OGconeOutputs_contrast%1.3f_pa%d_eye%d%d_eccen%1.2f_defocus%1.2f_noise-%s_sf%1.2f.mat',...
                        c,expParams.polarAngle,expParams.eyemovement(1,1),expParams.eyemovement(2,1), eccen, defocus, cMosaic.noiseFlag, sf);
                    if expParams.verbose;  fprintf('File will be saved as %s\n', fname); end
                    
                    % Update the stimulus contrast & spatial frequency
                    sparams.gabor.contrast  = c;  % Presumably michelson, [0 1]
                    sparams.freqCPD = sf;
                    
                    % ---- MAKE SCENE AND OIS --------------------------------------------
                    if expParams.verbose; fprintf('Recomputing scene for current sf and c..\n'); end
                    [OG,scenes,tseries] = ogStimuli(sparams);
                    
                    %             OG(1).visualize; % ccw
                    %             vcNewGraphWin;
                    %               plot(OG(1).timeAxis, OG(1).modulationFunction);
                    %               xlabel('Time (s)'); ylabel('Stimulus amplitude')
                    %
                    %             OG(2).visualize; % cw
                    %             vcNewGraphWin;
                    %               plot(OG(1).timeAxis, OG(1).modulationFunction);
                    %               xlabel('Time (s)'); ylabel('Stimulus amplitude')
                    
                    
                    % ------- Compute ABSORPTIONS -----------------------------------
                    
                    % Compute absorptions for multiple trials
                    absorptions = zeros(expParams.nTrials,cMosaic.rows,cMosaic.cols, cMosaic.tSamples, length(OG));
                    current     = absorptions;
                    
                    for s = 1:length(OG)
                        if currentFlag
                            [absorptions(:,:,:,:,s), current(:,:,:,:,s), interpFilters, meanCur] = cMosaic.compute(OG(s), 'currentFlag', currentFlag, ...
                                'emPaths', emPaths);
                        else
                            absorptions(:,:,:,:,s) = cMosaic.compute(OG(s), 'currentFlag', false, ...
                                'emPaths', emPaths);
                        end
                    end
                    
                    if expParams.verbose; fprintf('Saving data..\n'); end
                    savePth = fullfile(ogRootPath, 'data', expName); if ~exist(savePth,'dir'); mkdir(savePth); end;   
                    parsave(fullfile(savePth, fname), 'absorptions', absorptions, 'sparams', sparams, 'cparams', cparams, 'expParams', expParams, 'emPaths', emPaths);
                    if currentFlag; parsave(fullfile(savePth, ['current_' fname]), 'current', current, 'sparams', sparams, 'cparams', cparams, 'expParams', expParams, 'emPaths', emPaths); end
                end % sf
            end % contrast
        end % defocus
%     end % eyemovements
end % eccentricities

return
