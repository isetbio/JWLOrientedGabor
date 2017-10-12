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
% ----- These can be changed.
nTrials         = 25;        % Number of trials per stimulus condition
contrast_levels = [0:0.01:0.1];% 0.2:0.1:1];%([1:6 10])/100; % Contrast levels
eyemovement     = [1 1 0]';  % Which type of eye movments
eccentricities  = [0 2 5 10 20 40]; %4.5;
spatFreq        = [0.25, 0.4, 0.65, 1, 1.6, 2.6, 4, 8, 10, 16, 26];
polarangles     = 0;
defocuslevels   = 0;         % units??  [0 0.5 1 1.5 2]
verbose         = true;

% Temporal properties of one trial
tStep            = 0.002;                % Time step for optical image sequence (seconds)
sparams.tsamples = (0:tStep:0.054);      % seconds
sparams.timesd   = 1000.00;              % sd of temporal Gaussian window

% Scene field of view
sparams.sceneFOV  = 2;   % scene field of view in degrees (diameter)
sparams.freqCPD   = 4;   % Gabor spatial frequency (cpd)
sparams.gausSDdeg = sparams.sceneFOV/4; % Gabor SD in degrees of visual angle

deg2fov = 1/sparams.sceneFOV;
fov2deg = sparams.sceneFOV;

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

%% Loop over conditions, generating cone absorptions for each condition

% Eccentricity loop
for eccen = eccentricities
    
    % Polar angle loop
    for pa = polarangles
        
        % ----- CONE MOSAIC -----------------------------------------
        % Make CONE MOSAIC for a given eccentricity and polar angle
        whichEye = 'left';
        
        deg2m = 1/3 * 0.001; % 3 deg per mm, .001 mm per meter
        
        % Specify retinal location where stimulus is presented
        cparams.eccentricity = eccen;             % Visual angle of stimulus center, in deg
        polarAngleDeg        = pa;
        cparams.polarAngle   = deg2rad(polarAngleDeg);   % Polar angle (radians): 0 is right, pi/2 is superior, pi is left, 3*pi/2 inferior
        
        % Cone mosaic field of view in degrees
        cparams.cmFOV     = sparams.sceneFOV; % degrees
        
        % Compute x,y position in m of center of retinal patch from ecc and angle
        [x, y] = pol2cart(cparams.polarAngle, cparams.eccentricity);
        x = x * deg2m;  y = y * deg2m;
        
        cMosaic = coneMosaic('center', [x, y], 'whichEye', whichEye);
        
        % Set the field of view (degrees)
        cMosaic.setSizeToFOV(cparams.cmFOV);
        
        % Add photon noise
        cMosaic.noiseFlag = 'random'; % 'random' 'frozen' 'none'
        
        % ----- EYE MOVEMENTS -------------------------------------
        % Make EYE MOVEMENTS for a given cone mosaic
        
        % Not sure why these have to match, but there is a bug if they don't.
        cMosaic.integrationTime = OG(1).timeStep;
        
        % Include tremor, drift, microsaccades?
        cparams.em        = emCreate;
        cparams.em.emFlag = eyemovement;
        %         cparams.em.tremor.amplitude = cparams.em.tremor.amplitude * 2;
        %         cparams.em.drift.speed = cparams.em.drift.speed * 2;
        %         cparams.em.drift.speedSD = cparams.em.drift.speedSD * 2;
        
        % Generate the eye movement paths in units of cone samples. Set the
        % time sample to 2x the actual time, so that the eye does not
        % always start at 0,0. We will then clip after generating the eye
        % movements.
        emPaths  = cMosaic.emGenSequence(tSamples*2, 'nTrials', nTrials, ...
            'em', cparams.em); % path is in terms of cones shifted
        emPaths = emPaths(:, end/2+1:end,:);
        cMosaic.emPositions = squeeze(emPaths(1,:,:));
        
        if verbose
            %plot eye movements
            figure,
            subplot(211)
            plot(sparams.tsamples, emPaths(:,:,1)')
            
            subplot(212)
            plot(sparams.tsamples, emPaths(:,:,2)')
        end
        
        for defocus = defocuslevels
            
            % ---- Add optics blur or defocus if requested
            sparams.oi = oiDefocus(defocus); % input is Zernicke defocus coeff
            
            % Loop over contrasts and defocus
            for c = contrast_levels
                
                for sf = spatFreq
                    
                    fprintf('Computing absorptions for stimulus contrast %4.2f, polar angle %d, eccen %1.2f\n', c, pa, eccen)
                    fname = sprintf('OGconeOutputs_contrast%1.2f_pa%d_eye%d%d%d_eccen%1.2f_defocus%1.2f_noise-%s_sf%1.2f.mat',...
                        c,pa,cparams.em.emFlag(1),cparams.em.emFlag(2),cparams.em.emFlag(3), eccen, defocus, cMosaic.noiseFlag,sf);
                    fprintf('File will be saved as %s\n', fname);
                    
                    % Update the stimulus contrast & spatial frequency
                    sparams.gabor.contrast  = c;  % Presumably michelson, [0 1]
                    sparams.freqCPD = sf;
                    
                    % ---- MAKE SCENE AND OIS --------------------------------------------
                    disp('recomputing scene')
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
                    absorptions = zeros(nTrials,cMosaic.rows,cMosaic.cols, cMosaic.tSamples, length(OG));
                    %current     = absorptions;
                    
                    for s = 1:length(OG)
                        absorptions(:,:,:,:,s) = cMosaic.compute(OG(s), 'currentFlag', false, ...
                            'emPaths', emPaths);
                    end
                    
                    
                    
                    save(fullfile(ogRootPath, 'data', fname), 'absorptions', 'sparams', 'cparams');
                    
                end
            end
        end
    end
end

return






%% BIPOLAR LAYER

% Create a bipolar layer
bpL = bipolarLayer(cMosaic, 'nTrials', nTrials);

% Now make each type of bipolar mosaics.
bpCellTypes = {'on diffuse','off diffuse','on midget','off midget','on SBC'};

bpMosaicParams.rectifyType = 1;

% Store the multi-trial time series here (number of cell types x 2 stimuli)
bpNTrials = cell(length(bpCellTypes), 2);

for idx = 1:length(bpCellTypes)
    
    fprintf('Computing bipolar responses for cell type %s\n', ...
        bpCellTypes{idx});
    
    % How should we define this??
    %     bpMosaicParams.spread  = 2;  % RF diameter w.r.t. input samples
    %     bpMosaicParams.stride  = 1;  % RF diameter w.r.t. input samples
    
    % Add the current cell type as a new bipolar mosaic
    bpL.mosaic{idx} = bipolarMosaic(cMosaic, bpCellTypes{idx}, bpMosaicParams);
    
    % Set and compute are not yet implemented in the bipolarLayer
    %     bpL.mosaic{bpCellTypeInd}.set('sRFcenter',1);
    %     bpL.mosaic{bpCellTypeInd}.set('sRFsurround',1); % Do we want the surround to be 1?
    
    bpNTrials{idx,1} = bpL.mosaic{idx}.compute('coneTrials',current.cw);
    bpNTrials{idx,2} = bpL.mosaic{idx}.compute('coneTrials',current.ccw);
end

% Have a look
bpL.window;

% bpFilter = bipolarFilter(bp{1}, cMosaic,'graph',true);
% vcNewGraphWin; plot(cMosaic.timeAxis,bpFilter,'o-');


%% Retinal ganglion cell model

% We need to think about RGC sampling in two senses: RF size and RGC
% density

% Create retina ganglion cell layer object based on bipolar layer
rgcL = rgcLayer(bpL);

% Choose cell types
rgcCellTypes = {'on parasol','off parasol','on midget','off midget'};

diameters = round([10 10 2 2 20]); % In cones. % Should we have to manually set this?

%diameters = round([
rgcParams.name = 'macaque inner retina 1'; % ?? Not sure about this: Do we want macaque or human retina?
rgcParams.eyeSide = whichEye;

% Do we want to use these parameters?
rgcParams.centerNoise = 0.2;
rgcParams.ellipseParams = [1 .8 0];  % Principle, minor and theta
rgcParams.axisVariance = .1;

% Loop over celltypes to create RGC mosaics in the rgc layer
for idx = 1:length(rgcCellTypes)
    rgcParams.rfDiameter = diameters(idx);
    rgcL.mosaic{idx} = rgcLNP(rgcL, bpL.mosaic{idx},rgcCellTypes{idx},rgcParams);
    
end

% Set nr of trials / repeats of the same stimulus
rgcL.set('numberTrials',nTrials);


%% Compute the inner retina response and visualize

% Number of trials refers to number of repeats of the same stimulus
disp('Computing rgc responses');
% rgcL = rgcL.compute;
nTrialsSpikes = rgcL.compute('bipolarScale',50,'bipolarContrast',0.2,'bipolarTrials',bpNTrials,'coupling',false);

%
% [rgcL, nTrialsSpikes] = rgcL.compute('bipolarTrials',bpNTrials,'coupling',false);
% % [rgcL, nTrialsSpikes] = rgcL.compute(bpL.mosaic{idx},'bipolarTrials',bpNTrials,'coupling',false);
% % [rgcL, nTrialsSpikes] = rgcL.compute(bpL.mosaic{1},'bipolarTrials',bpNTrials,'coupling',false,'bipolarScale',50,'bipolarContrast',0.2);
% % [rgcL, nTrialsSpikes] = rgcL.compute(bpL.mosaic,'bipolarTrials',bpNTrials,'bipolarScale',50,'bipolarContrast',0.2);


% Have a look
rgcL.window;


return
