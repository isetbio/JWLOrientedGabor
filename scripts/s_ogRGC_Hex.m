%% s_ogRGC_Hex

% Oriented Gabor Discrimination but now with Hexagonal mosaic
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

% A temporary offsheet of s_ogRGC_JW for testing purposes. Merge or delete
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
contrast_levels = 0.1;% [0:0.01:0.1];% 0.2:0.1:1];%([1:6 10])/100; % Contrast levels
eyemovement     = [1 1 0]';  % Which type of eye movments
eccentricities  = 4.5;
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
        %         cparams.cmFOV     = sparams.sceneFOV; % degrees
        
        % Compute x,y position in m of center of retinal patch from ecc and angle
        [x, y] = pol2cart(cparams.polarAngle, cparams.eccentricity);
        x = x * deg2m;  y = y * deg2m;
        
        
        %%
        
        % Mosaic Parameters
        quality.tolerance1 = 0.5;                           % larger than default tolerances to speed-up computation. For production work, either do not set, or set to equal or lower than 0.01
        quality.tolerance2 = 0.05;                          % larger than default tolerances to speed-up computation, For production work, either do not set, or set to equal or lower than 0.001
        
        %% Set import/export options
        saveMosaic = false;                                 % whether to save the mosaic
        loadMosaic = false;                                 % whether to load a previously saved mosaic
        saveMosaicPDF = false;                              % whether to save a PDF of the mosaic
        
        
        hexecc = sqrt(sum(x.^2));
        hexang = atan2d(y,x);
        [spacing, aperture, density, params, comment] = coneSizeReadData('eccentricity', hexecc,'eccentricityUnits','m','angle', hexang,'angleUnits','deg','whichEye', whichEye);
        
        
        cparams = struct('name', sprintf('HexMosaic_eccen%2.1f_pa%2.1f',eccen, pa), ...
            'resamplingFactor', 8, ...                        % controls the accuracy of the hex mosaic grid
            'eccBasedConeDensity', false, ...                  % whether to have an eccentricity based, spatially - varying density
            'sConeMinDistanceFactor', [],...3.0, ...                % Min distance between neighboring S-cones = f * local cone separation - used to make the S-cone lattice semi-regular
            'sConeFreeRadiusMicrons', [],...0.15*300, ...
            'customLambda', spacing*10.^6, ...                           % custom spacing?
            'centerInM', [x y], ...                           % mosaic eccentricity
            'spatialDensity', [0 6/10 3/10 1/10],...          % with a LMS density of of 0.6:0.3:0.1
            'fovDegs', sparams.sceneFOV, ...
            'latticeAdjustmentPositionalToleranceF', quality.tolerance1, ...
            'latticeAdjustmentDelaunayToleranceF', quality.tolerance2);
        
        theHexMosaic = coneMosaicHex(cparams.resamplingFactor, ...
            'name', cparams.name, ...
            'eccBasedConeDensity', cparams.eccBasedConeDensity, ...
            'sConeMinDistanceFactor', cparams.sConeMinDistanceFactor, ...
            'sConeFreeRadiusMicrons', cparams.sConeFreeRadiusMicrons, ...
            'customLambda', cparams.customLambda, ...
            'spatialDensity', cparams.spatialDensity,...
            'fovDegs', cparams.fovDegs, ...
            'latticeAdjustmentPositionalToleranceF', cparams.latticeAdjustmentPositionalToleranceF, ...
            'latticeAdjustmentDelaunayToleranceF', cparams.latticeAdjustmentDelaunayToleranceF);
        
        theHexMosaic.window;
        % Save the mosaic for later analysis
        if saveMosaic
            save(fullfile(ogRootPath, 'data', [cparams.name '.mat']), 'theHexMosaic','-v7.3')
        end
        
        %         cMosaic = coneMosaic('center', [x, y], 'whichEye', whichEye);
        
        % Set the field of view (degrees)
        %         cMosaic.setSizeToFOV(cparams.cmFOV);
        
        % Add photon noise
        theHexMosaic.noiseFlag = 'random'; % 'random' 'frozen' 'none'
        
        % ----- EYE MOVEMENTS -------------------------------------
        % Make EYE MOVEMENTS for a given cone mosaic
        
        % Not sure why these have to match, but there is a bug if they don't.
        theHexMosaic.integrationTime = OG(1).timeStep;
        
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
        emPaths  = theHexMosaic.emGenSequence(tSamples*2, 'nTrials', nTrials, ...
            'em', cparams.em); % path is in terms of cones shifted
        emPaths = emPaths(:, end/2+1:end,:);
        theHexMosaic.emPositions = squeeze(emPaths(1,:,:));
        
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
                
                if verbose; fprintf('Computing absorptions for stimulus contrast %4.2f, polar angle %d, eccen %1.2f\n', c, pa, eccen); end
                fname = sprintf('OGconeOutputs_contrast%1.2f_pa%d_eye%d%d%d_eccen%1.2f_defocus%1.2f_noise-%s',...
                    c,pa,cparams.em.emFlag(1),cparams.em.emFlag(2),cparams.em.emFlag(3), eccen, defocus, theHexMosaic.noiseFlag);
                if verbose; fprintf('File will be saved as %s\n', fname); end
                
                % Update the stimulus contrast
                sparams.gabor.contrast  = c;  % Presumably michelson, [0 1]
                
                
                % ---- MAKE SCENE AND OIS --------------------------------------------
                disp('recomputing scene')
                [OG,scenes,tseries] = ogStimuli(sparams);
                
                
                % ------- Compute ABSORPTIONS -----------------------------------
                
                % Compute absorptions for multiple trials
                %                 absorptions = zeros(nTrials,theHexMosaic.mosaicSize(1),theHexMosaic.mosaicSize(2), theHexMosaic.tSamples, length(OG));
                %current     = absorptions;
                
                nonNullConeIndices = theHexMosaic.pattern > 1;
                
                
                tic
                parfor s = 1:length(OG)
                    
                    absorptions2D = zeros(nTrials, size(nonNullConeIndices,1), size(nonNullConeIndices,2) ,tSamples);

                    absorptions = theHexMosaic.compute(OG(s), 'currentFlag', false, ...
                        'emPaths', emPaths);
                    
                    
                    for trial = 1:size(absorptions,1)
                        for sample = 1:size(absorptions,3)
                            
                            % initiate matrix
                            theseAbsorptions = zeros(size(nonNullConeIndices));
                            % put absorptions back into native 2d cone mosaic
                            theseAbsorptions(nonNullConeIndices==1) = squeeze(absorptions(trial,:,sample));
                            
                            absorptions2D(trial, :, :, sample) = theseAbsorptions;
                        end
                        
%                         [activationsHexImage(trial,:,:,:), activationsLMSHexImage(trial,:,:,:)] = theHexMosaic.computeActivationDensityMap(squeeze(absorptions2D(trial,:,:,:)));
                        %
                    end
                    
                    % Reduce file size
                    absorptions2D = single(absorption2D);
                    
                    if verbose; fprintf('Saving absorptions for OG%d \n', s); end
                    parsave(fullfile(ogRootPath, 'data', sprintf([fname '_OG%d.mat'], s)), 'absorptions2D', absorptions2D, 'sparams',sparams, 'cparams',cparams);
                    
                end
                toc
                
                %                 visualizedAperture = 'lightCollectingArea'; % choose between 'both', 'lightCollectingArea', 'geometricArea'
                %                 theHexMosaic.visualizeGrid(...
                %                     'visualizedConeAperture', visualizedAperture, ...
                %                     'apertureShape', 'disks', ...
                % %                     'panelPosition', [1 1], 'generateNewFigure', true);
                % %
                %                 nonNullConeIndices = theHexMosaic.pattern > 1;
                %                 absorptions2D = zeros(nTrials, size(nonNullConeIndices),tSamples);
                %
                %                 absorptions2D(nonNullConeIndices==1) = squeeze(absorptions(1,:,1));
                %
                %                 [activationsHexImage, activationsLMSHexImage] = theHexMosaic.computeActivationDensityMap(absorptions2D);
                % %                 figure;
                %                 imagesc(1:theHexMosaic.cols, 1:theHexMosaic.rows, squeeze(activationsLMSHexImage(:,:,1)));
                
                
                
                
            end
        end
    end
end

return