%% s_exampleScript

% Testing script for isetbio (a stand alone simple experiment)
% by EK (May 29, 2018)

%% 1. Set up scene parameters

nTrials = 5;

% Scene field of view
sparams.fov      = 2;   % scene field of view in degrees (diameter)
sparams.distance = 0.57;    % Meters

% Gabor parameters
P(1:2)        = harmonicP;  % Standard Gabor
P(1).contrast = 0;

% Temporal properties of one trial
dur      = 54;   % ms
% deltaT   = 2;    % ms
gaussSd  = 15;   % ms

% modulation amplitude when blending gabor and blank
minA     = 0;
maxA     = 1; 
tseries  = ieScale(fspecial('gaussian',[1,dur],gaussSd),minA,maxA);

% pad tseries with zeros
tseries = [tseries zeros(1,dur*3)];

% Make a default Optical Image Sequence.
[ois, scene] = oisCreate('harmonic','blend', tseries, 'testParameters', P, 'sceneParameters', sparams);

% We make sure that the number of time points in the eye movement sequence
% matches the number of time points in the optical image sequence
tSamples         = ois.length;

% Check ois time axis, i.e. number of frames
oiTimeAxis       = ois.timeAxis;
nFrames          = length(oiTimeAxis);

%% 2. Make cone mosaic for a given eccentricity and polar angle

% Pick an eye
whichEye = 'left';

% Convert degrees to meters
deg2m = 0.3 * 0.001; % (adjusted to the default in isetbio)

% Specify retinal location where stimulus is presented
cparams.eccentricity = 10;             % Visual angle of stimulus center, in deg
polarAngleDeg        = 0;              % Polar angle of stimulus, in radians, so right from fixation, horizontal meridian
cparams.polarAngle   = deg2rad(polarAngleDeg);   % Polar angle (radians): 0 is right, pi/2 is superior, pi is left, 3*pi/2 inferior

% Cone mosaic field of view in degrees
cparams.cmFOV        = sparams.fov; % degrees

% Compute x,y position in m of center of retinal patch from ecc and angle
[x, y] = pol2cart(cparams.polarAngle, cparams.eccentricity);
x = x * deg2m;  y = y * deg2m;

% Create cone mosaic
cMosaic = coneMosaic('center', [x, y], 'whichEye', whichEye);

% Set the field of view (degrees)
cMosaic.setSizeToFOV(cparams.cmFOV);

% Add photon noise
cMosaic.noiseFlag = 'random';

% Set outer segment to be computed with linear filters
cMosaic.os = osLinear;

cMosaic.integrationTime = ois.timeStep;
eyeMovementsNum = ois.maxEyeMovementsNumGivenIntegrationTime(cMosaic.integrationTime);

% Change proportion covered to scale the size of the cones accordingly to
% the simulated eccentricity
% propCovered = getBanks1991ConeCoverage(cparams.eccentricity);
% cMosaic.pigment.pdWidth  = cMosaic.pigment.width*propCovered;
% cMosaic.pigment.pdHeight = cMosaic.pigment.height*propCovered;


%% 3. Add eyemovements

emPaths  = cMosaic.emGenSequence(eyeMovementsNum, 'nTrials', nTrials); % path is in terms of cones shifted


%% 4. Compute absorptions and cone current

[absorptions, current, interpFilters, meanCur] = cMosaic.computeForOISequence(ois, 'currentFlag', true, 'emPaths', emPaths);


%% 5. Visualize

% Check outputs
% cMosaic.window;

% Define response axis
responseTimeAxis = cMosaic.timeAxis;

% Indices to 3 cone classes
lcones = cMosaic.pattern == 2;
mcones = cMosaic.pattern == 3;
scones = cMosaic.pattern == 4;

absorptions1 = permute(absorptions, [2 3 4 1]);
current1 = permute(current, [2 3 4 1]);
sz = size(absorptions1);
absorptions1 = reshape(absorptions1, sz(1)*sz(2), sz(3), sz(4));
current1 = reshape(current1, sz(1)*sz(2), sz(3), sz(4));

% Compute absorptions and current per cone class as cones x tPoints x trials
l_absorptions = absorptions1(lcones,:,:);
m_absorptions = absorptions1(mcones,:,:);
s_absorptions = absorptions1(scones,:,:);

l_curr = current1(lcones,:,:);
m_curr = current1(mcones,:,:);
s_curr = current1(scones,:,:);



% Summarize each cone class as the mean time series across cones and across
% trials
l_absorptions_mn = mean(mean(l_absorptions,3));
m_absorptions_mn = mean(mean(m_absorptions,3));
s_absorptions_mn = mean(mean(s_absorptions,3));

l_current_mn = mean(mean(l_curr,3));
m_current_mn = mean(mean(m_curr,3));
s_current_mn = mean(mean(s_curr,3));

% plot
figure; 
subplot(4,1,1)
plot(oiTimeAxis,ois.modulationFunction,'k')
hold on; grid on
ylabel('magnitude')
xlabel('Time (s)')
title('Stimulus modulation function')

subplot(4,1,2)
plot(responseTimeAxis,l_absorptions_mn,'r', 'LineWidth', 1.5)
hold on; grid on
ylabel('Photon absorption count')
xlabel('Time (s)')
title('L cones')

subplot(4,1,3)
plot(responseTimeAxis,m_absorptions_mn,'g', 'LineWidth', 1.5)
hold on; grid on
ylabel('Photon absorption count')
xlabel('Time (s)')
title('M cones')

subplot(4,1,4)
plot(responseTimeAxis,s_absorptions_mn,'b','LineWidth', 1.5)
hold on; grid on
ylabel('Photon absorption count')
xlabel('Time (s)')
title('s cones')


%% plot Currents
figure; 
subplot(4,1,1)
plot(oiTimeAxis,ois.modulationFunction,'k')
grid on
ylabel('magnitude')
xlabel('Time (s)')
title('Stimulus modulation function')

subplot(4,1,2); hold on; grid on
plot(responseTimeAxis,l_current_mn,'r', 'LineWidth', 1.5)
ylabel('pAmps')
xlabel('Time (s)')
title('L cones')

subplot(4,1,3); hold on; grid on
plot(responseTimeAxis,m_current_mn,'g', 'LineWidth', 1.5)
ylabel('pAmps')
xlabel('Time (s)')
title('M cones')

subplot(4,1,4); hold on; grid on
plot(responseTimeAxis,s_current_mn,'b','LineWidth', 1.5)
ylabel('pAmps')
xlabel('Time (s)')
title('s cones')

