%% t_colorOpponentGabor
%
% Model cone responses for the experiment on detection of Gabor color
% opponent stimuli in 
%
%  Chromatic detection from cone photoreceptors to V1 neurons to behavior
%  in rhesus monkeys
%   Horwitz, Hass, Angueyra, Lindbloom-Brown & Rieke,
%   J. Neuroscience, 2015 (submitted)
%
% Script outline:
%     1. Build scene as a color-opponent Gabor patch with Horwitz Lab
%         display and imageHarmonicColor (to integrate with
%         imageHarmonicColor). 
%     2. Build oi and sensor with specified properties from paper (pupil
%         width, macular density, eccentricity, etc.).
%     3. Create a dynamic scene by altering the Gabor phase and
%         regenerating the scene and recomputing the oi and sensor.
%     4. Build the outer segment object with linear filters and compute its
%         response with the dynamic sensor structure.
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
%     6. This process is first done for zero contrast with no Gabor
%         present, and repeated for a range of contrasts specified in the
%         contrastArr variable. The pooled response for the cone array at
%         each contrast level is compared to the pooled response for no
%         Gabor present using a linear SVM, and the cross-validated
%         accuracy of the linear SVM is calculated.
%      7. The accuracy is plotted as a function of the contrast, and a
%         psychometric function of the form
%
%               p = 1 - 0.5*exp(-(x/alpha)^beta)
%
%         is fit to the curve, where alpha is the threshold for detection
%         and beta is the slope of the curve.
%
% JRG/NC/BW ISETBIO Team, Copyright 2015


ieInit
clear all
ieSessionSet('wait bar', 'off')
wFlag = ieSessionGet('wait bar');

%% Specify parameters for contrast values and noise repititions

nContrast   = 4;
maxContrast = linspace(0,1,nContrast);

% contrastArr = [0,1];% 
% contrastArr = (0:0.25:1); % must start with 0
% nContrast = length(contrastArr);

noiseIterations = 100;    % more iterations will improve accuracy but take longer!
pooledData      = cell(1,nContrast);
rocArea         = zeros(1,nContrast);

% Load the display from the Horwitz Lab
display = displayCreate('CRT-Sony-HorwitzLab');

% Basic human optics parameters.  Perhaps we need to get the macaque optics
% built here.
oi  = oiCreate('wvf human');

%% Loop over color and contrast to build the response level
%
% The color index chooses the stimulus color.
% The four possible values from the original paper are
%  1 = S-cone isolating
%  2 = L-M opponent
%  3 = L-M+S
%  4 = L-M-S
for colorInd = 2 
    for contrastInd = 1:nContrast   % contrastInd = 1;
        tic
        
        % Set parameters for building the dynamic scene/oi/sensor The
        % stimulus consists of a color opponent Gabor patch whose phase
        % varies over time. The stimulus parameters specified in Fig. 6 of
        % the manuscript. 
        
        % parameters found in Fig. 6 caption
        params = paramsGaborColorOpponent();
        params.contrast = 0;                  % set max contrast of gabor
        params.color = colorInd;              % 1 = s_iso, 2 = L-M, 3 = LMS, 4 = L-M
        params.image_size = 64;               % scene is (image_size X image_size) pixels
        params.fov        = 0.6;
        params.nSteps     = 60; % 666; % 160+346+160
        params.contrast = maxContrast(contrastInd);
        
        % This is the contrast at each step in the simulation.  It ramps up
        % to the max (1 in this case), stays there, and then ramps back
        % down.
%         timeContrast = ...
%             [linspace(0,1,(1/4)*params.nSteps), ...
%              ones(1,params.nSteps/2),...
%              linspace(1,0,(1/4)*params.nSteps)];
        timeContrast = ones(params.nSteps,1);
        % plot(timeContrast)
            
        % We build a dummy scene here just so we can subsequently calculate
        % the sensor size.  But this scene itself is not used.  Rather we
        % build a series of scenes below.
        stimulusRGBdata = imageHarmonicColor(params); % sceneCreateGabor(params);
        scene = sceneFromFile(stimulusRGBdata, 'rgb', params.meanLuminance, display);
        vcAddObject(scene); sceneWindow;
        
        %% Initialize the optics and the sensor
              
        coneP = coneCreate; % The cone properties properties
        
        % see caption for Fig. 4 of Horwitz, Hass, Rieke, 2015, J. Neuro.
%         retinalPosDegAz = 5; retinalPosDegEl = -3.5;
%         retinalRadiusDegrees = sqrt(retinalPosDegAz^2+retinalPosDegEl^2);
%         retinalPolarDegrees  = abs(atand(retinalPosDegEl/retinalPosDegAz));
%         retinalPos = [retinalRadiusDegrees retinalPolarDegrees]; whichEye = 'right';
%         sensor = sensorCreate('human', [coneP], [retinalPos], [whichEye]);
        sensor = sensorCreate('human');
        sensor = sensorSetSizeToFOV(sensor, params.fov, scene, oi);
        sensor = sensorSet(sensor, 'exp time', params.expTime); % 1 ms
        sensor = sensorSet(sensor, 'time interval', params.timeInterval); % 1 ms
        
        % This computes with no sensor or photon noise, just the mean
        sensor = sensorSet(sensor,'noise flag',0);

        % According to the paper, cone collecting area is 0.6 um^2
        % sensor = sensorSet(sensor 'pixel pd width', 0.774e-6); % photo-detector width
        % sensor = sensorSet(sensor, 'pixel pd height', 0.774e-6);
       
        % Macular pigment is attached to the sensor
        % macular pigment absorbance was scaled to 0.35 at 460 nm
        % macular = sensorGet(sensor, 'human macular');
        % macular = macularSet(macular, 'density', 0.35);
        % sensor = sensorSet(sensor, 'human macular', macular);
        
        % Compute a dynamic set of cone absorptions
        % ieSessionSet('wait bar',true);
        fprintf('Computing dynamic scene/oi/sensor data    ');
        
        if wFlag, wbar = waitbar(0,'Stimulus movie'); end
        
        % Loop through frames to build movie
        for t = 1 : params.nSteps
            
            fprintf('\b\b\b%02d%%', round(100*t/params.nSteps));
            
            if wFlag, waitbar(t/params.nSteps,wbar); end
            
            % Update the phase of the Gabor
            params.ph = 2*pi*(t-1)/params.nSteps; % one period over nSteps
            params.contrast = params.contrast*timeContrast(t);
            
            % scene = sceneCreate('harmonic', params);
            % scene = sceneSet(scene, 'h fov', fov);
            
            % The imageHarmonicColor should be moved to ISETBIO.
            stimulusRGBdata = imageHarmonicColor(params); % sceneCreateGabor(params);
            
            scene = sceneFromFile(stimulusRGBdata, 'rgb', params.meanLuminance, display);
            scene = sceneSet(scene, 'h fov', params.fov);
            
            % The mean scene luminance is set to 200 cd/m2, based on the
            % calibration of the monitor
            scene = sceneAdjustLuminance(scene, 200);
            
            % oi  = oiCreate('wvf human');
            % Compute optical image
            oi = oiCompute(oi, scene);
            
            % Compute absorptions
            sensor = sensorCompute(sensor, oi);
            
            if t == 1
                volts = zeros([sensorGet(sensor, 'size') params.nSteps]);

            end
            
            volts(:,:,t) = sensorGet(sensor, 'volts');
            
            % In case we need to visualize, you can do this on each
            % iteration
            % vcAddObject(scene); sceneWindow
            % pause(.1);
        end
        fprintf('\n');
        toc
        if wFlag, delete(wbar); end
        
        % Set the stimuls into the sensor object
        sensor = sensorSet(sensor, 'volts', volts);
        % vcAddObject(sensor); sensorWindow;
        
        %% Train linear SVM and find cross-validated accuracy
        % Create the outer segment object
        os = osCreate('linear');
        
        % Compute the photocurrent for the whole time series
        os = osCompute(os, sensor);
        
        % Pool all of the noisy responses across each cone type
        % nNoiseIterations by nConeTypes
        tic
        pooledData{contrastInd} = pooledConeResponse(os, sensor, noiseIterations);
        toc
        if contrastInd > 1
            
            % Visualize pooled responses in LMS space
            vcNewGraphWin;
            scatter3(pooledData{1}(:,1),pooledData{1}(:,2),pooledData{1}(:,3))
            hold on; 
            scatter3(pooledData{contrastInd}(:,1),pooledData{contrastInd}(:,2),pooledData{contrastInd}(:,3))
            xlabel('pooled S cone response'); ylabel('pooled M cone response'); zlabel('pooled L cone response');
            
            % Fit a linear svm classifier between pooled responses at contrast = 0
            % and contrast = maxContrast(contrastInd):
            
            % If you have the statistics toolbox
            if exist('fitcsvm','file');
                m1 = fitcsvm([pooledData{1}; pooledData{contrastInd}], [ones(noiseIterations,1); -1*ones(noiseIterations,1)], 'KernelFunction', 'linear');
                % Calculate cross-validated accuracy based on model:
                cv = crossval(m1,'kfold',5);
                rocArea(contrastInd) = 1-kfoldLoss(cv)'
            else
                % Use the ISETBIO libsvm 
                nFoldCrossVal = 5;
                % NEEDS TO BE CHECKED.   ASK HJ.
                pd1 = pooledData{1}; pd2 = pooledData{contrastInd};
                [acc,w] = svmClassifyAcc([pd1; pd2], ...
                    [ones(noiseIterations,1); -1*ones(noiseIterations,1)], ...
                    nFoldCrossVal,'linear');
                % perfcurve for roc plot
                rocArea(contrastInd) = acc(1)
            end
            
            title(sprintf('Pooled responses in LMS space, p(Correct) = %2.0f', 100*rocArea(contrastInd)));
            set(gca,'fontsize',14);
        end
        
        % clear sensor scene oi display params os
    end %colorInd
    
    %% Fit psychometric curve to thresholds as a function of contrast
    [xData, yData] = prepareCurveData( maxContrast, rocArea);
    
    % Set up fittype and options.
    ft = fittype( '1 - 0.5*exp(-(x/a)^b)', 'independent', 'x', 'dependent', 'y' );
    opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
    opts.Display = 'Off';
    opts.StartPoint = [0.323369521886293 0.976303691832645];
    
    % Fit a curve between contrast level (x) and probability of correction
    % detection.
    % bootWeibullFit(stimLevels, nCorrect, nTrials, varargin)
    % in computationaleyebrain/simulations/Pixel Visibility/ ...
    [fitresult, gof] = fit( xData, yData, ft, opts );
    
    % Plot fit with data.
    figure( 'Name', 'untitled fit 1' );
    h = plot( fitresult, xData, yData );
    hold on; scatter(maxContrast,rocArea,40,'filled');
    % set(gca,'xscale','log')
    legend( h, 'data', 'fitted curve', 'Location', 'NorthWest' );
    % Label axes
    xlabel Contrast
    ylabel p(Correct)
    grid on
    thresh1 = fitresult.a;
    title(sprintf('Detection, \\alpha = %1.2f',(thresh1)));
    set(gca,'fontsize',16)
    axis([0 1 0.5 1]);
    

end%contrastInd

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
