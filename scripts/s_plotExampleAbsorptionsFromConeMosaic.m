%% s_plotExampleAbsorptionsFromConeMosaic

% Script to plot the absorptions of a stimulus for 3 cone densities

%% 0. Define parameters

% Experiment params
expName           = 'eccbasedcoverage';
expParams         = loadExpParams(expName, false);

% Plotting params
eccentricities    = [0, 4.5, 40]; % deg
whichEye          = 'left';
cparams.cmFOV     =  2; % degrees
climsA            = [0 220];  % color bar / ylims for absorptions (photon count)

% Convertion deg to m
deg2m  = 1/3 * 0.001; % 3 deg per mm, .001 mm per meter


figure(1); hold all;
set(gcf, 'Color', 'w', 'Position', [1000, 681, 1173, 657], 'NumberTitle', 'off', 'Name', 'Figure 8A - 2D absorptions for 3 cone density levels');


for eccen = eccentricities 
    
        tStep             = 0.002;                % Time step for optical image sequence (seconds)
        sparams.tsamples  = (0:tStep:(0.054));      % seconds
        
        % Scene field of view
        sparams.sceneFOV  = 2;   % scene field of view in degrees (diameter)
        sparams.freqCPD   = 4;   % Gabor spatial frequency (cpd)
        sparams.gausSDdeg = sparams.sceneFOV/4; % Gabor SD in degrees of visual angle

        % Unit converters
        deg2fov = 1/sparams.sceneFOV;
        fov2deg = sparams.sceneFOV;

        % Gabor parameters
        sparams.gabor           = harmonicP;                   % Standard Gabor
        sparams.gabor.ang       = (pi/180)* 15;                % Gabor orientation (radians) - question: what is 0??
        sparams.gabor.freq      = fov2deg*sparams.freqCPD;     % Spatial frequency (cycles/FOV)
        sparams.gabor.contrast  = 1;                           % Presumably michelson, [0 1]
        sparams.gabor.GaborFlag = sparams.gausSDdeg*deg2fov;   % Gaussian window

        % Add defocus
        sparams.oi = oiDefocus(0);

        % Create Optical Image Sequence
        OG = ogStimuli(sparams);
    
        % Specify retinal location where stimulus is presented
        cparams.eccentricity = eccen;             % Visual angle of stimulus center, in deg
        cparams.polarAngle   = expParams.polarAngle;   % Polar angle (radians): 0 is right, pi/2 is superior, pi is left, 3*pi/2 inferior
        
        % Compute x,y position in m of center of retinal patch from ecc and angle
        [x, y] = pol2cart(cparams.polarAngle, cparams.eccentricity);
        x = x * deg2m;  y = y * deg2m;
        
        cMosaic = coneMosaic('center', [x, y], 'whichEye', whichEye);
        
        % Set the field of view (degrees)
        cMosaic.setSizeToFOV(cparams.cmFOV);
        
        % Set proportion covered
        propCovered = getBanks1991ConeCoverage(eccen, 'verbose', false);
        cMosaic.pigment.pdWidth  = cMosaic.pigment.width*propCovered;
        cMosaic.pigment.pdHeight = cMosaic.pigment.height*propCovered;
       

        % Create empty array of eye movements
        emPaths = zeros(expParams.nTrials, length(sparams.tsamples), 2);

        % Compute absorptions for one stimulus
        absorptions = cMosaic.compute(OG(1), 'currentFlag', false, 'emPaths', emPaths, 'seed', 1);
        

        subplot(1,3, find(eccen==eccentricities)); hold all;
        imagesc(squeeze(absorptions(1,:,:,1))./2);
        colormap gray; axis image; axis off; colorbar;
        title(sprintf('Cone density %d', eccen2density(cMosaic, 'deg')) ,'Fontsize',12)

end