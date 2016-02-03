%% t_orientedGaborDiscrimination
%
% Model cone responses for the experiment on measuring orientation
% discrimination thresholds of an acrhomatic, peripheral Gabor by Marisa
% Carrasco and Jon Winawer.
%
% Script outline:
%     1. Build scene as an achromatic Gabor patch with imageHarmonic.
%     2. Build oi and sensor with specified properties (eccentricity, etc.)
%     3. Loop over 2 stimulus orientations and n trials per orientation
%     4. Within each trial, add eye movements to sensor 
%     5. Build the outer segment object with linear filters and compute its
%         response with the sensor structure.
%     6. Sum the response over time for each cone as a cheap way to simulate
%          downstream temporal integration of cone outputs 
%     7. The result of each trial is a vector of cone responses. The result
%         for each stimulus orientation is a matrix of trial by cone responses
%         Theses responses are compared using a linear SVM, and the
%         cross-validated accuracy of the linear SVM is calculated.
%     
%     TODO: spatial integration of cone outputs prior to classification
%
% JRG/NC/BW ISETBIO Team, Copyright 2015

clear
ieInit

plotConeResponseFlag = false; % plot an example movie of cone responses during one trial
saveConeCurrentsFlag = true;  % save cone responses across all trials (after temporal integration within trials)

%% Specify parameters for contrast values and noise repititions

% Specify angles of Gabor stimulus orientation
angleArr = (pi/180)*[-20 20]; nAngles = length(angleArr);

% Specify retinal location where stimulus is presented (6 deg above fovea)
locationRadius = 6; % degrees of visual angle
locationAngle  = 0; % degrees of polar angle (0 is vertical)

nTrials        = 50;    % number of trials per stimulus condition

% Basic human optics parameters. 
oi  = oiCreate('wvf human');

%% Initialize the optics and the sensor

% Set parameters for building the scene/oi/sensor The
% stimulus consists of an achromatic Gabor patch at +/- 20 deg.

% stimulus parameters
params            = paramsGaborColorOpponent();

params.fov        = 1.5;            % Gabor is windowed with 1.5 deg aperture
params.freq       = 6*params.fov;   % ... and has spatial frequency of 6 cpd 
params.GaborFlag  = .25/params.fov; % ... and std of .25 deg
params.nSteps     = 100;            % 100 time steps (1 ms sampling)
params.ecc        = 6;              % 6 degrees eccentricity
params.contrast   = .25;            % Max contrat of 0.25


% We build a dummy scene here just so we can subsequently calculate
% the sensor size.  But this scene itself is not used.  Rather we
% build the scene below.
scene = sceneCreate('harmonic', params);
scene = sceneSet(scene, 'h fov', params.fov);


coneP = coneCreate; % The cone properties properties

sensor = sensorCreate('human', coneP, [locationRadius(1) locationAngle(1)]);
sensor = sensorSetSizeToFOV(sensor, params.fov, scene, oi);
sensor = sensorSet(sensor, 'exp time', params.expTime); % 1 ms
sensor = sensorSet(sensor, 'time interval', params.timeInterval); % 1 ms

% This computes with sensor and photon noise
sensor = sensorSet(sensor,'noise flag',2);

%% Loop over color and contrast to build the response level

storedConeCurrents = cell(1,nAngles); 
for angleInd = 1:nAngles    % +/- angleArr(angleInd) degs
    fprintf('\n'); 

    % We render a new scene for each stimulus orientation (but not for
    % repeated trials for the same orientation). We also pad the FOV by
    % factor padFactor so that eye movements do not cause the sensor to
    % move out of the scence FOV
    padFactor             = 2;
    theseParams           = params;
    theseParams.ang       = angleArr(angleInd);
    theseParams.freq      = params.freq * padFactor;
    theseParams.GaborFlag = params.GaborFlag / padFactor;
    
    scene = sceneCreate('harmonic', theseParams);
    scene = sceneSet(scene, 'h fov', params.fov*padFactor);
    
    % The mean scene luminance is set to 200 cd/m2, based on the calibration of the monitor
    scene = sceneAdjustLuminance(scene, 200);
    % vcAddAndSelectObject(scene); sceneWindow;
    
    % Compute optical image
    oi = oiCompute(oi, scene);
    % vcAddAndSelectObject(oi); oiWindow;
    
    for trial = 1:nTrials
        fprintf('.'); drawnow();
        
        eyePos = randn(params.nSteps, 2); % eye position in units of number of cones
        
        % for debugging, try a translation or square wave oscillations
        % eyePos =  (1:params.nSteps)' * [1 0];
        % eyePos =  square((1:params.nSteps)/params.nSteps*2*pi*5)' * [1 0];
        
        % Loop through frames to build movie
        for t = 1 : params.nSteps
            
            % Compute absorptions
            sensor = sensorSet(sensor, 'positions', eyePos(t, :));
            sensor = coneAbsorptions(sensor, oi);
            
            if t == 1
                volts = zeros([sensorGet(sensor, 'size') params.nSteps]);
            end
            volts(:,:,t) = sensorGet(sensor, 'volts');
        end % t
        
        % Set the stimuls into the sensor object
        sensor = sensorSet(sensor, 'volts', volts);
        
        
        %% Train linear SVM and find cross-validated accuracy
        % Create the outer segment object
        os = osCreate('linear');
        
        % Compute the photocurrent for the whole time series
        os = osCompute(os, sensor);
        coneCurrentSignal = sum(osGet(os, 'conecurrentsignal'),3);

% %         To implement RGC spatial pooling, go to isetbio rgc branch and
% %         uncomment to line 160: 
%         % Create rgc object for spatial pooling
%         clear params
%         params.name    = 'Macaque inner retina 1'; % This instance
%         params.model   = 'pool';    % Computational model
%         params.row     = sensorGet(sensor,'row');  % N row samples
%         params.col     = sensorGet(sensor,'col');  % N col samples
%         params.spacing = sensorGet(sensor,'width','um'); % Cone width
%         params.timing  = sensorGet(sensor,'time interval','sec'); % Temporal sampling
%         params.eyeSide   = 'left';   % Which eye
%         params.eyeRadius = 2;        % Radius in mm
%         params.eyeAngle  = 90;       % Polar angle in degrees
%         
%         rgc1 = rgcCreate(params);
%         for cellTypeInd = 1:5%length(obj.mosaic)
%             rgcSet(rgc1, 'mosaic', rgcMosaicPool(rgc1));
%         end
%         
%         rgc1 = rgcCompute(rgc1, os);
%         for cellTypeInd = 1:5
%             rgcSignal{cellTypeInd} = cellfun(@sum,mosaicGet(rgc1.mosaic{cellTypeInd}, 'linearResponse'),'uniformoutput',false);
%         end
        
        storedConeCurrents{angleInd}(trial,:) = coneCurrentSignal(:);
        
        % visualize cone response for an example trials
        if plotConeResponseFlag && trial == 1
            vcAddObject(scene); sceneWindow;
            
            vcNewGraphWin;
            coneCurrentSignal = osGet(os, 'conecurrentsignal');
            
            subplot(3,1,1)
            imagesc(mean(coneCurrentSignal,3));axis image
            for ii = 1:params.nSteps
                subplot(3,1,2)
                imagesc(coneCurrentSignal(:,:,ii)); title(ii), axis image
                subplot(3,1,3)
                imagesc(volts(:,:,ii)); title(ii), axis image
                pause(0.1);
            end % nSteps
        end
        
    end
   
end %angleInd

coneData = [storedConeCurrents{1}; storedConeCurrents{2}];
labels   = [ones(19,1); -1*ones(19,1)];

if saveConeCurrentsFlag
    savePth = fullfile(fileparts(which(mfilename)), 'data');
    save(fullfile(savePth, sprintf('coneResponses%s', datestr(now, 'YYYY-mm-DD_HH:MM:SS'))), ...
        'coneData', 'labels');
end
% Fit a linear svm classifier between two orientations

m1 = fitcsvm(coneData, labels, 'KernelFunction', 'linear');

% Calculate cross-validated accuracy based on model:
cv = crossval(m1,'kfold',5);
rocArea = 1-kfoldLoss(cv)';

return




%% Fit psychometric curve to thresholds as a function of contrast
[xData, yData] = prepareCurveData( maxContrast', squeeze(rocArea(angleInd,locationInd,:)));

% Set up fittype and options.
ft = fittype( '1 - 0.5*exp(-(x/a)^b)', 'independent', 'x', 'dependent', 'y' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.StartPoint = [0.3 0.95];

% Fit a curve between contrast level (x) and probability of correction
% detection.
% bootWeibullFit(stimLevels, nCorrect, nTrials, varargin)
% in computationaleyebrain/simulations/Pixel Visibility/ ...
[fitresult, gof] = fit( xData, yData, ft, opts );

% Plot fit with data.
figure( 'Name', 'untitled fit 1' );

h = plot( fitresult, xData, yData );
hold on; scatter(maxContrast,squeeze(rocArea(angleInd,locationInd,:)),40,'filled');
% set(gca,'xscale','log')
legend( h, 'data', 'fitted curve', 'Location', 'NorthWest' );
% Label axes
xlabel Contrast
ylabel p(Correct)
grid on
thresh1 = fitresult.a;
title(sprintf('Discrimination for Angle = %2.2f Location %d\nContrast \\alpha = %1.2f',(180/pi)*angleArr(angleInd),locationInd,(thresh1)));
set(gca,'fontsize',16)
axis([0 1 0.5 1]);

%%



%% Show the movie of volts
%
% % This should be a sensorMovie() call.
% %
% % Can we easily make that movie when we color the cones by type
% vcNewGraphWin;axis image; colormap(gray)
% for ii=1:params.nSteps
%     imagesc(volts(:,:,ii)); pause(.2);
% end

% % Time series at a point
% vcNewGraphWin; plot(squeeze(volts(1,1,:)))
%
% %% Movie of the cone absorptions over cone mosaic
% % from t_VernierCones by HM
%
% step = 1;
% tmp = coneImageActivity(sensor,[],step,false);
%
% % Show the movie
% vcNewGraphWin;
% tmp = tmp/max(tmp(:));
% for ii=1:size(tmp,4)
%     img = squeeze(tmp(:,:,:,ii));
%     imshow(img.^3); truesize;
%     title('Cone absorptions')
%     drawnow
% end
%
% %% Outer segment calculation
%
% % The outer segment converts cone absorptions into cone photocurrent.
% % There are 'linear','biophys' and 'identity' types of conversion.  The
% % linear is a standard convolution.  The biophys is based on Rieke's
% % biophysical work.  And identity is a copy operation.
% os = osCreate('linear');
%
% % Compute the photocurrent
% os = osCompute(os, sensor);
%
% % Plot the photocurrent for a pixel
% % Let's JG and BW mess around with various plotting things to check the
% % validity.
% osPlot(os,sensor);
%
% % Input = RGB
% % os = osCreate('identity');
% % os = osSet(os, 'rgbData', sceneRGB);
%
% %% Rieke biophysics case
%
% os = osCreate('biophys');
%
% % Compute the photocurrent
% os = osCompute(os, sensor);
%
% % Plot the photocurrent for a pixel
% % Let's JG and BW mess around with various plotting things to check the
% % validity.
% osPlot(os,sensor,'output')
%
% %% Build rgc
%
% eyeAngle = 180; % degrees
% eyeRadius = 3; % mm
% eyeSide = 'right';
% rgc1 = rgcCreate('GLM', scene, sensor, os, eyeSide, eyeRadius, eyeAngle);
%
% rgc1 = rgcCompute(rgc1, os);
%
% % rgcPlot(rgc1, 'mosaic');
% % rgcPlot(rgc1, 'linearResponse');
% rgcPlot(rgc1, 'spikeResponse');
% %% Build rgc response movie
% %  https://youtu.be/R4YQCTZi7s8
%
% % % osLinear
% % rgcMovie(rgc1, sensor);
%
% % % osIdentity
% % rgcMovie(rgc1, os);
%
%
