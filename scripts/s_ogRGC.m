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
%       * Repeated trials are handled awkwardly. For example, repeated
%           trials are stored within separate cell arrays for inner retina
%           mosaics (which is good), but are not stored at all within
%           coneMosais or bipolar mosaics. As a result, we need to store
%           the repeated trial data in separate variables, not attached to
%           the data strucutres, which causes complications when passing
%           one data structure (say a bp) into another (eg, innerretina).
%       * (SOLVED) Why does inner retina allow multiple mosaics whereas bp does not?
%       * (SOLVED) How to combine the bipolar responses into the RGC layer?
%       * (SOLVED) How should we deal with the negative responses from off midget
%       RGC that turn into no signal after rectification?
%       * TODO: Oncorporate off-axis optical data from Pablo Artal
%               Example: Figure 5 from Jaeken and Artal, 2012
%       * How is the denisty of cones, bipolars, and RGCs handled as a
%               function of eccentricity?
%           It appears that:
%               Cone density is based on Curcio. (fine)
%               Bipolar mosaic is exactly matched to the cone mosaic (fine)
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

% Number of trials per stimulus condition
nTrials  = 50;


%% SCENE AND OPTICAL IMAGE SEQUENCE
for eccen = 0:1:40
    for c = 1%[0.01:0.01:0.1, 0.2:0.1:1]%, 0.2:0.1:1]
        
        for pa = 0 % [0 90 180 270]
            fprintf('Computing absorptions for stimulus contrast %4.2f, polar angle %d, eccen %1.2f\n', c, pa, eccen)
            % ---- SCENE PARAMETERS ---------------------------------------------
            % Gaussian temporal window for stimulus
            tStep            = 0.002;                % Time step for optical image sequence (seconds)
            sparams.tsamples = (-0.200:tStep:0.200); % seconds
            sparams.timesd   = 0.100;                % sd of temporal Gaussian window
            
            % Scene field of view
            sparams.sceneFOV  = 2;   % scene field of view in degrees (diameter)
            sparams.freqCPD   = 6;   % Gabor spatial frequency (cpd)
            sparams.gausSDdeg = .25; % Gabor SD in degrees of visual angle
            
            deg2fov = 1/sparams.sceneFOV;
            fov2deg = sparams.sceneFOV;
            
            % Make a Gabor with default parameters, then update parameters. Will need
            %   to fit within sceneFOV
            %
            % Note that the angle convention for harmonics is different that spatial
            % position for the cone mosaic. For the harmonic, 0 is vertical, pi/4 is 4
            % deg up and right.
            sparams.gabor           = harmonicP;                   % Standard Gabor
            sparams.gabor.ang       = (pi/180)* 20;                % Gabor orientation (radians) - question: what is 0??
            sparams.gabor.freq      = fov2deg*sparams.freqCPD;     % Spatial frequency (cycles/FOV)
            sparams.gabor.contrast  = c;                           % Presumably michelson, [0 1]
            sparams.gabor.GaborFlag = sparams.gausSDdeg*deg2fov;   % Gaussian window
            
            % ---- MAKE SCENE AND OIS --------------------------------------------
            [OG,scenes,tseries, fname] = ogStimuli(sparams);
            
%             OG(1).visualize; % ccw
%             vcNewGraphWin;
%               plot(OG(1).timeAxis, OG(1).modulationFunction);
%               xlabel('Time (s)'); ylabel('Stimulus amplitude')
%             
%             OG(2).visualize; % cw
%             vcNewGraphWin;
%               plot(OG(1).timeAxis, OG(1).modulationFunction);
%               xlabel('Time (s)'); ylabel('Stimulus amplitude')
            
            %% CONE MOSAIC
            whichEye = 'left';
            
            deg2m = 1/3 * 0.001; % 3 deg per mm, .001 mm per meter
            
            % Specify retinal location where stimulus is presented
            cparams.eccentricity = eccen;             % Visual angle of stimulus center, in deg
            polarAngleDeg        = pa;
            cparams.polarAngle   = deg2rad(polarAngleDeg);   % Polar angle (radians): 0 is right, pi/2 is superior, pi is left, 3*pi/2 inferior
            
            % Cone mosaic field of view in degrees
            cparams.cmFOV     = 2; % degrees
            
            % Compute x,y position in m of center of retinal patch from ecc and angle
            [x, y] = pol2cart(cparams.polarAngle(1), cparams.eccentricity(1));
            x = x * deg2m;  y = y * deg2m;
            
            cMosaic = coneMosaic('center', [x, y], 'whichEye', whichEye);
            
            % Set the field of view (degrees)
            cMosaic.setSizeToFOV(cparams.cmFOV);
            
            % Add photon noise
            cMosaic.noiseFlag = 'random';
            
            %% EYE MOVEMENTS
            
            % NOTE
            %   Eye movements are currently specified in units of cones, not deg or m.
            %   If cone density changes with visual field position, then we would
            %   implicitly be assuming different size eye movements for different
            %   visual field positions.
            
            % We make sure that the number of time points in the eye movement sequence
            % matches the number of time points in the optical image sequence
            tSamples         = OG(1).length;
            
            % Not sure why these have to match, but there is a bug if they don't.
            cMosaic.integrationTime = OG(1).timeStep;
            
            
            % ----- EYE MOVEMENTS -----------------------------
            cparams.em        = emCreate;    % eye movements: consider adjusting to
            %   account for cone spacing and for data
            %   from different stimulus conditions
            cparams.em.emFlag = [1 1 0]';    % Include tremor, drift, microsaccades
            
            emPaths  = cMosaic.emGenSequence(tSamples, 'nTrials', nTrials, ...
                'em', cparams.em); % path is in terms of cones shifted
            
            
            %% ABSORPTIONS
            
            % Compute absorptions for multiple trials
            
            % ccw Gabor
            [absorptions.ccw, current.ccw, interpFilters.ccw, meanCur.ccw] = cMosaic.compute(OG(1), 'currentFlag', true, ...
                'emPaths', emPaths);
            
            % cw Gabor
            [absorptions.cw, current.cw, interpFilters.cw, meanCur.cw] = cMosaic.compute(OG(2), 'currentFlag', true, ...
                'emPaths', emPaths);
            
            % Have a look
            % cMosaic.window;
            
            % plot the mean absorptions and current
            %        sz = cMosaic.rows*cMosaic.cols;
            %        figure, plot(max(reshape(cMosaic.current, sz,[]))); hold on;
            %        plot(min(reshape(cMosaic.current, sz,[])));
            %        title('current')
            
            %        figure, plot(max(reshape(cMosaic.absorptions, sz,[]))); hold on;
            %        plot(min(reshape(cMosaic.absorptions, sz,[])));
            %        title('absorptions')
            
            save(fullfile(ogRootPath, 'data', sprintf('OGconeOutputs_contrast%1.2f_pa%d_eye%d%d%d_eccen%1.2f.mat',c,pa,cparams.em.emFlag(1),cparams.em.emFlag(2),cparams.em.emFlag(3), eccen)),...
                'absorptions', 'sparams', 'cparams');
            
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


%% VERSION 1 code:
%Loop over two stimulus classes and repeated trials and compute cone responses
%%
storedConeCurrents = cell(1,length(OG));

for angleInd = 1:length(OG)
    fprintf('\n'); % Comment out for publishing
    
    % We render a new scene for each stimulus orientation (but not for
    % repeated trials for the same orientation). We also pad the scene
    % field of view by factor padFactor so that eye movements do not cause
    % the sensor to move out of the scence field of view
    padFactor             = 2;
    theseParams           = params;
    theseParams.ang       = params.orientation(angleInd);
    theseParams.freq      = params.freq * padFactor;
    theseParams.GaborFlag = params.GaborFlag / padFactor;
    
    scene = sceneCreate('harmonic', theseParams);
    scene = sceneSet(scene, 'h fov', params.fov*padFactor);
    
    % Compute optical image - also computed once per stimulus class (not
    % re-computed for new trials within a stimulus class)
    oi = oiCompute(oi, scene);
    
    % Loop over trials. Each trial gets its own eye movments.
    for trial = 1:nTrials
        % fprintf('.'); drawnow();  % Comment out for publishing
        
        % Eye position in units of number of cones
        eyePos = eyePosFun(params.nSteps);
        
        % Loop through frames to build movie
        volts = zeros([sensorGet(sensor, 'size') params.nSteps]);
        for t = 1 : params.nSteps
            
            % Compute absorptions
            sensor = sensorSet(sensor, 'positions', eyePos(t, :));
            sensor = coneAbsorptions(sensor, oi);
            
            volts(:,:,t) = sensorGet(sensor, 'volts');
        end % t
        
        % Set the stimuls into the sensor object
        sensor = sensorSet(sensor, 'volts', volts);
        %% Train linear SVM and find cross-validated accuracy
        % Create the outer segment object
        %%
        os = osCreate('linear');
        
        % Compute the photocurrent for the whole time series
        os = osCompute(os, sensor);
        coneCurrentSignal = sum(osGet(os, 'conecurrentsignal'),3);
        params.coneArraySize = size(coneCurrentSignal);
        
        storedConeCurrents{angleInd}(trial,:) = coneCurrentSignal(:);
        
        % visualize cone response for an example trials
        if plotConeResponseFlag && trial == 1
            vcAddObject(scene); sceneWindow;
            
            vcNewGraphWin;
            coneCurrentSignal = osGet(os, 'conecurrentsignal');
            
            subplot(3,1,1)
            imagesc(mean(coneCurrentSignal,3));axis image
            title('Mean cone output across trial');
            for ii = 1:params.nSteps
                subplot(3,1,2)
                imagesc(coneCurrentSignal(:,:,ii)); axis image
                title('Cone output at single time points');
                subplot(3,1,3)
                imagesc(volts(:,:,ii));  axis image
                title('Cone voltage at single time points');
                pause(0.1);
            end % nSteps
        end
        
    end
    
end %angleInd
%% Save
%%
coneData = [storedConeCurrents{1}; storedConeCurrents{2}];
labels   = [ones(nTrials,1); -1*ones(nTrials,1)];

dataPth = fullfile(fileparts(which(mfilename)), 'data');

if saveConeCurrentsFlag
    fname = fullfile(dataPth, sprintf('coneResponses%s.mat', datestr(now, 'YYYY-mm-DD_HH.MM.SS')));
    save(fname,  'coneData', 'labels', 'params', 'eyePosFun');
end


return
%% Classifiy
%%
meridianName = cell(1,4);
cellnum = NaN(1,4);
rocStored = [];
rocSEM = [];
for meridian = 1:4
    
    d = dir(fullfile(dataPth, '*.mat'));
    [~,idx] = sort([d.datenum]);
    
    load(fullfile(dataPth, d(idx(end-meridian+1)).name));
    
    cellnum(meridian) = size(coneData,2);
    angles(meridian) = params.polarAngle;
    switch params.polarAngle
        case 0,     meridianName{meridian} = 'Nasal (HM)';
        case 90,    meridianName{meridian} = 'Superior (LVM)';
        case 180,   meridianName{meridian} = 'Temporal (HM)';
        case 270,   meridianName{meridian} = 'Inferior (UVM)';
    end
    
    % Fit a linear svm classifier between two orientations and calculate
    % cross-validated accuracy based on model
    
    % -------- first for cone signals ---------------------------------------
    m = fitcsvm(coneData, labels, 'KernelFunction', 'linear');
    cv = crossval(m,'kfold',5);
    rocAreaCones = 1-kfoldLoss(cv);
    
    % -------- then for RGC signals -----------------------------------------
    %  RGC hack - bandpass filter and subsample
    
    % RGC is a strucutred array. Each element corresponds to one model of RGCs.
    % The model is defined by a spatial receptive field 'rf' and a subsampling
    % rate 'ss'. The output of the RGCs are stored in 'data'. The inputs come
    % from the cone outputs.
    rgc = struct('ss', [], 'rf', [], 'data', []);
    
    % We will make several RGC classes defined by scaleFactor, which will scale
    % both the receptive field size and the subsampling rate (the bigger the
    % RF, the more coarsely we subsample).
    scaleFactor = linspace(1.4,1.8,4);
    scaleFactor = linspace(1.8,4, 3);%logspace(log10(1), log10(10),10);
    scaleFactor = [2 3];
    %scaleFactor = linspace(1.4,1.8,4);
    for ii = 1:length(scaleFactor)
        % Center surround RF. The receptive fields have an excitatory center
        % with std of 1 cone * scaleFactor. The inhibitory surround has twice
        % the sd of the center. The cell is balanced, so that the sum of center
        % minus surround = 0.
        rfCenter   = fspecial('Gaussian', 20,scaleFactor(ii));
        rfSurround = fspecial('Gaussian', 20,2*scaleFactor(ii));
        rgc(ii).rf =  rfCenter - rfSurround;
        rgc(ii).subsample = scaleFactor(ii);
    end
    
    %  Load stored cone data
    coneDataR = reshape(coneData, [], params.coneArraySize(1), params.coneArraySize(2));
    
    
    %  Loop across RGC types
    nBtsrp = 20;
    rocAreaRGC = NaN(nBtsrp, length(rgc));
    for ii = 1:length(rgc)
        ss  = rgc(ii).subsample;
        rf  = rgc(ii).rf;
        
        % Loop across trials, extracting RGC signals by convolution + subsample
        for trial = 1:size(coneDataR,1)
            coneCurrentSignal = squeeze(coneDataR(trial,:,:));
            tmp = conv2(coneCurrentSignal,rf, 'valid');
            tmp = tmp(round(ss:ss:end),round(ss:ss:end));
            rgc(ii).data(trial,:) = tmp(:);
        end
        
        % Classify RGC outputs
        for btstrp = 1:nBtsrp
            % scrambled = labels(randperm(length(labels)));
            m = fitcsvm(rgc(ii).data, labels, 'KernelFunction', 'linear');
            cv = crossval(m,'kfold',5);
            rocAreaRGC(btstrp, ii) = 1-kfoldLoss(cv);
        end
    end
    
    fprintf('ROC Area for cones: %4.2f\n', rocAreaCones)
    for ii = 1:length(rgc)
        fprintf('ROC Area for RGC class %d: %4.2f\n', ii, mean(rocAreaRGC(:,ii)))
    end
    
    % Plot classification accuracy as function of convergence ratio
    if meridian == 1, fH = vcNewGraphWin; end
    set(fH, 'Color', 'w')
    set(gca,'FontSize',20); hold on
    errorbar(scaleFactor, mean(rocAreaRGC), std(rocAreaRGC),'o-', 'LineWidth', 4, 'MarkerSize', 12)
    xlabel('Cone to Midget Ganglion Cell Ratio')
    ylabel('Classification Accuracy')
    
    rocStored(meridian,:) = mean(rocAreaRGC);
    rocSEM(meridian,:) = std(rocAreaRGC);
end

legend(meridianName)
%%
%
%%
fH = vcNewGraphWin;    set(fH, 'Color', 'w')
p = polar(0,.5, 'k');
hold on
set(gca,'FontSize',20);
set(p, 'visible', 'off');
for ii = 1:size(rocStored,2)
    p = polar(deg2rad([angles angles(1)]), rocStored([1:end 1],ii)'-.5, 'o-');
    set(p, 'LineWidth', 4, 'MarkerSize', 12)
    errorbar(0, rocStored(angles==90,ii)-.5, rocSEM(angles==90,ii), 'LineWidth', 3, 'Color', get(p, 'color'));
end
hgexport(gcf, sprintf('~/Desktop/polarPF%s.eps', whichEye))
title('Classification Accuracy')
