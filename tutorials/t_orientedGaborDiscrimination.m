%% t_orientedGaborDiscrimination
%
% Model cone responses for the experiment on measuring orientation
% discrimination thresholds of an acrhomatic, peripheral Gabor by Marisa
% Carrasco and Jon Winawer.
%
% Script outline:
%     1. Build scene as an achromatic Gabor patch with Horwitz Lab
%         display and imageHarmonic.
%     2. Build oi and sensor with specified properties  (pupil
%         width, macular density, eccentricity, etc.).
%     3. Add realistic eye movements to sensor using data from the
%         experiment.
%     4. Build the outer segment object with linear filters and compute its
%         response with the sensor structure.
%     5. Compute the response across each cone type.  For each noise
%         iteration (a) project the time series of each cone onto the
%         expected mean time series, (b) take that mean value for each cone
%         type to give the mean on that particular noise iteration. The
%         projection step eliminates signals that are orthogonal to the
%         expected linear temporal response.  The expected response is
%         derived by running the simulation with the noise flag set to 0
%         (photon noise only). This calculation produces an (L, M, S)
%         triplet for each noise iteration.  This is a specific
%         implementation of an ideal observer; see the paper for details.
%     6. This process is carried out for a range of retinal locations with
%         stimuli at +/- 20 deg of vertical and varying levels of
%         contrast/noise. The pooled response for the cone array at the two
%         orientations are comparedusing a linear SVM, and the
%         cross-validated accuracy of the linear SVM is calculated.
%      7. The accuracy is plotted as a function of the contrast, and a
%         psychometric function of the form
%
%               p = 1 - 0.5*exp(-(x/alpha)^beta)
%
%         is fit to the curve, where alpha is the threshold for detection
%         and beta is the slope of the curve.
%
% JRG/NC/BW ISETBIO Team, Copyright 2015

clear
ieInit
ieSessionSet('wait bar', 'off')
wFlag = ieSessionGet('wait bar');

%% Specify parameters for contrast values and noise repititions

% Specify angles of Gabor stimulus orientation
angleArr = (pi/180)*[-20 20]; nAngles = length(angleArr);

% Specify retinal location where stimulus is presented
locationRadius = 6; % degrees of visual angle
locationAngle  = 0; % degrees of polar angle (0 is vertical)

nTrials        = 50;    % number of trials per stimulus condition

% Basic human optics parameters.  Perhaps we need to get the macaque optics
% built here.
oi  = oiCreate('wvf human');

%% Initialize the optics and the sensor

% Set parameters for building the scene/oi/sensor The
% stimulus consists of an achromatic Gabor patch at +/- 20 deg.

% stimulus parameters
params            = paramsGaborColorOpponent();


params.fov        = 1.5;
params.freq       = 6;
params.GaborFlag  = .25; % standard deviation of the Gaussian window (per image)
params.nSteps     = 100; % 666; % 160+346+160
params.ecc        = 6;
params.contrast   = 0.25;

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
    params.ang        = angleArr(angleInd);
    
    for trial = 1:nTrials
        fprintf('.'); drawnow();
        
        eyePos = randn(params.nSteps, 10)*3; % eye position in units of number of cones
        %eyePos =  (1:params.nSteps)' * [1 1];
        % eyePos = zeros(size(eyePosI));
        
        scene = sceneCreate('harmonic', params);
        scene = sceneSet(scene, 'h fov', params.fov);
        
        % The mean scene luminance is set to 200 cd/m2, based on the calibration of the monitor
        scene = sceneAdjustLuminance(scene, 200);
        % vcAddAndSelectObject(scene); sceneWindow;
        
        % oi  = oiCreate('wvf human');
        % Compute optical image
        oi = oiCompute(oi, scene);
        % vcAddAndSelectObject(oi); oiWindow;
        
        % Loop through frames to build movie
        for t = 1 : params.nSteps
            
            % Compute absorptions
            sensor = sensorSet(sensor, 'positions', eyePos(t, :));
            sensor = coneAbsorptions(sensor, oi);
            
            if t == 1
                volts = zeros([sensorGet(sensor, 'size') params.nSteps]);
            end
            volts(:,:,t) = sensorGet(sensor, 'volts');
            % vcAddObject(sensor); sensorWindow;
        end % t
        
        % Set the stimuls into the sensor object
        sensor = sensorSet(sensor, 'volts', volts);
        
        
        % vcAddObject(scene); sceneWindow;
        % vcAddObject(sensor); sensorWindow;
        
        %% Train linear SVM and find cross-validated accuracy
        % Create the outer segment object
        os = osCreate('linear');
        
        % Compute the photocurrent for the whole time series
        os = osCompute(os, sensor);
        coneCurrentSignal = sum(osGet(os, 'conecurrentsignal'),3);
        
        storedConeCurrents{angleInd}(trial,:) = coneCurrentSignal(:);
        
    end
    % visualize
    if 0
        fH = vcNewGraphWin;
        coneCurrentSignal = osGet(os, 'conecurrentsignal');
        
        volts = sensorGet(sensor, 'volts');
        waitforbuttonpress;
        for ii = 1:params.nSteps
            subplot(2,1,1)
            imagesc(coneCurrentSignal(:,:,ii)); title(ii), axis image
            subplot(2,1,2)
            imagesc(volts(:,:,ii)); title(ii), axis image
            pause(0.1);
        end % nSteps
    end
end %angleInd

        
% Fit a linear svm classifier between two orientations

m1 = fitcsvm([storedConeCurrents{1}; storedConeCurrents{2}], [ones(nTrials,1); -1*ones(nTrials,1)], 'KernelFunction', 'linear');

% Calculate cross-validated accuracy based on model:
cv = crossval(m1,'kfold',5);
rocArea = 1-kfoldLoss(cv)';

return

title(sprintf('Pooled responses in LMS space, p(Correct) = %2.0f', 100*rocArea));
set(gca,'fontsize',14);



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
