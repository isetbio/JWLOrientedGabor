%% Oriented Gabor Discrimination, take 2
% This script supercedes t_orientedGaborDiscrimination.m, which no longer runs, 
% due to changes in isetbio.
% 
% Model cone, bipolar, RGC, cortical, and behavioral (computational
% observer) responses for the experiment on measuring orientation
% discrimination thresholds of an achromatic, peripheral Gabor.
% 
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
%       * Spatial and temporal pooling for cortex?
%       * Repeated trials are handled awkardly. For example, repeated
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

%% ---------------------- Experiment ----------------------
%
%       *   Loop over 4 spatial positions (left, right, lower, upper)
%       *   For each position, loop over 2 stimulus orientations and n 
%               trials per orientation
%       *   Within each trial, add eye movements to sensor (or OIS?)
%       *   Build the outer segment object with linear filters and compute its
%              response with the sensor structure.
%       *   Build the bipolar layer and compute responses (on and off??)
%       *   Build the RGC layer and compute responses (on and off midgets?)
%       *   Compute cortical responses from spatial / temporal responses
%       *   Compute linear discriminant of stimulus orientation using
%              cortical responses
% 

% EK/JW/ NYU ISETBIO Team, Copyright 2017


% % Cone density; see coneDensity.m
% % load cone density data.  Units are cones / mm^2
% d = load('coneDensity.mat');
% figure,  plot(d.inferior.eccMM, d.inferior.density)
% hold on, plot(d.superior.eccMM, d.superior.density)
% hold on, plot(d.temporal.eccMM, d.temporal.density)
% hold on, plot(d.nasal.eccMM, d.nasal.density)
% set(gca, 'YScale', 'log')
% legend('Inferior', 'Superior', 'Temporal', 'Nasal')


%% Specify experiment parameters 

% Number of trials per stimulus condition
nTrials  = 10;
whichEye = 'left'; 


%% Specify stimulus parameters

% Gaussian temporal window for stimulus
tStep                      = 0.002; % Time step for making optical image sequence (seconds)
params.tsamples            = (-0.070:tStep:0.070); % seconds
params.timesd              = 0.100;                % seconds

% Scene field of view
params.sceneFOV  = 2;           % degrees

% Make a Gabor with default parameters, then update parameters. Will need
%   to fit within sceneFOV
params.gabor               = harmonicP;                   % Standard Gabor
params.gabor.ang           = (pi/180)* 20;                % Gabor orientation (radians)
params.gabor.freq          = 6*params.sceneFOV;           % Spatial frequency (cycles/deg)
params.gabor.contrast      = 1;                           % Presumably michelson, [0 1]
params.gabor.GaborFlag     = .25/params.sceneFOV;         % Gaussian window

% Specify retinal location where stimulus is presented
params.eccentricity        = 6;                           % Visual angle of stimulus center, in deg
params.polarAngle          = 90;                          % Polar angle (deg): 0 is right, 90 is superior, 180 is left, 270 inferior


%% Make the stimuli


[OG,scenes,tseries, fname] = ogStimuli(params);

% OG(1).visualize;
% vcNewGraphWin; 
%   plot(OG(1).timeAxis, OG(1).modulationFunction); 
%   xlabel('Time (s)'); ylabel('Stimulus amplitude')


%% Compute absorptions from multiple tirals

% Compute absorptions for multiple trials
tSamples         = OG(1).length;

% Cone mosaic field of view in degrees
params.cmFOV     = 2;           % degrees
params.em        = emCreate;    % eye movements: consider adjusting to 
                                %   account for cone spacing and for data
                                %   from different stimulus conditions
params.em.emFlag = [1 1 1]';    % Include tremor, drift, microsaccades

% Compute x,y position in m of center of retinal patch from ecc and angle
[x, y] = pol2cart(params.polarAngle(1), params.eccentricity(1));
x = x * .3 * 0.001; % .3 mm per deg, .001 mm per meter
y = y * .3 * 0.001; % .3 mm per deg, .001 mm per meter

% Create cone mosaic
cMosaic = coneMosaic('center', [x, y], 'whichEye', whichEye);

% Sometimes we set the mosaic size to 15 minutes (.25 deg) because that is
% the spatial pooling size found by Westheimer and McKee
cMosaic.setSizeToFOV(params.cmFOV);

% Not sure why these have to match, but there is a bug if they don't.
cMosaic.integrationTime = OG(1).timeStep;

% Add photon noise
cMosaic.noiseFlag = 'random';

% Add eye movements
emPaths  = cMosaic.emGenSequence(tSamples, 'nTrials', nTrials, ...
    'em', params.em); % path is in terms of cones shifted

% Compute absorptions 
[absorptions, current, interpFilters, meanCur] = cMosaic.compute(OG(1), 'currentFlag', true, ...
    'emPaths', emPaths);

% Have a look
% cMosaic.window;

%% Add bipolar cells

% Run two bipolar cell models
bpCellTypes = {'onmidget' 'offmidget'};

for bpCellTypeInd = 1:length(bpCellTypes)
    
    clear bpParams
    bpParams.cellType = bpCellTypes{bpCellTypeInd};
    bpParams.ecc = cMosaic.center(1);
    bpParams.rectifyType = 1;
    
    % Add the current cell type as a new bipolar mosaic
    bpMosaic{bpCellTypeInd} = bipolar(cMosaic, bpParams);
    bpMosaic{bpCellTypeInd}.set('sRFcenter',1);
    bpMosaic{bpCellTypeInd}.set('sRFsurround',1); % Do we want the surround to be 1?
    
    [~, bpNTrialsCenterTemp, bpNTrialsSurroundTemp] = bpMosaic{bpCellTypeInd}.compute(cMosaic,'coneTrials',current);
    bpNTrials{bpCellTypeInd} = bpNTrialsCenterTemp-bpNTrialsSurroundTemp;
    clear bpNTrialsCenterTemp bpNTrialsSurroundTemp
end

% Have a look
% bp{1}.window;
% bpFilter = bipolarFilter(bp{1}, cMosaic,'graph',true);
% vcNewGraphWin; plot(cMosaic.timeAxis,bpFilter,'o-');


%% Retinal ganglion cell model

% Choose a cell type
irCellTypes = {'onMidget' 'offMidget'};     %'onMidget'; %'OFF Midget';  % 'offParasol'; 'onMidget' ...
irParams.name = 'macaque inner retina 1'; % ?? Not sure about this: Do we want macaque or human retina?
irParams.eyeSide = whichEye;

% Create inner retina object
ecc = params.eccentricity(1); % Check if ecc should be in degrees, radius or m or mm?
ecc = cMosaic.center(1) * 1000; % in mm now..
% ecc = deg2rad(params.eccentricity(1));

irParams.eyeRadius = sqrt(sum(ecc.^2)); 
irParams.eyeAngle = 0;

% Compute inner retina with bipolar (bp) cell outputs
innerRetina = ir(bpMosaic, irParams);

% Do we want to use these parameters?
mosaicParams.centerNoise = 0.2;
mosaicParams.ellipseParams = [1 .8 0];  % Principle, minor and theta
mosaicParams.axisVariance = .1;
mosaicParams.model = 'lnp';


for irCellTypeInd = 1:length(irCellTypes)
    mosaicParams.type = irCellTypes{irCellTypeInd};
    innerRetina.mosaicCreate(mosaicParams);
end

% Set nr of trials / repeats of the same stimulus
innerRetina.set('numberTrials',nTrials);


%% Compute the inner retina response and visualize

% Number of trials refers to number of repeats of the same stimulus
disp('Computing rgc responses');
[innerRetina, nTrialsSpikes] = innerRetina.compute(bpMosaic,'bipolarTrials',bpNTrials,'coupling',false); 

% Have a look
innerRetina.mosaic{1}.window;
innerRetina.mosaic{2}.window;


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