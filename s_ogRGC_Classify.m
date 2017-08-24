%% s_ogRGC_Classify

% Script with first attempt to classify oriented gabors simualted at 4
% different polar angles

%% ------- 0. Define parameters -------

% Define paths
projectPth = '~/matlab/git/toolboxes/JWLOrientedGabor';
dataPth = fullfile(projectPth, 'data');
saveData = true;

% Define Gaussian temporal window for stimulus
tStep                      = 0.002; % Time step for making optical image sequence (seconds)
sparams.tsamples            = (-0.200:tStep:0.200); % seconds
sparams.timesd              = 0.100;                % sd of temporal Gaussian window

% Define scene field of view, spatial frequency and Gabor SD
sparams.sceneFOV  = 2;           % degrees (Diameter??)
sparams.freqCPD   = 6;           % Gabor spatial frequency (cpd)
sparams.gausSDdeg = .25;         % Gabor SD in degrees of visual angle

% Unit converters
deg2fov = 1/sparams.sceneFOV;
fov2deg = sparams.sceneFOV;
deg2m = .3 * 0.001; % .3 mm per deg, .001 mm per meter

% Polar angles to simulate
polarAngles = [0 90 180 270]; % 0 is nasal (HM), 90 is superior (LVM), 180 is temporal (HM), 270 is inferior (UVM)

% Number of trials per stimulus condition
nTrials  = 50;

% Which eye?
whichEye = 'left';

% Predefine cells to store data
storedAbsorptions = cell(2, length(polarAngles));
storedCurrents    = cell(2, length(polarAngles));


%% ------- 1. CREATE SCENE AND OPTICAL IMAGE SEQUENCE -------

% Define Gabor parameters
sparams.gabor           = harmonicP;                   % Standard Gabor
sparams.gabor.ang       = (pi/180)* 20;                % Gabor orientation (radians) - question: what is 0??
sparams.gabor.freq      = fov2deg*sparams.freqCPD;     % Spatial frequency (cycles/FOV)
sparams.gabor.contrast  = 1;                           % Presumably michelson, [0 1]
sparams.gabor.GaborFlag = sparams.gausSDdeg*deg2fov;   % Gaussian window

% ---- MAKE SCENE AND OIS --------------------------------------------
[OG,scenes,tseries, fname] = ogStimuli(sparams);

%% Set cone mosaic params

% Specify retinal location where stimulus is presented (might be a parameter we want to loop
% over later)
cparams.eccentricity = 6;             % Visual angle of stimulus center, in deg

% Cone mosaic field of view in degrees
cparams.cmFOV     = 2; % degrees

%% Loop over both OG's
for orientationInd = 1:2
    
    % Get an OG
    thisOG = OG(orientationInd);
    
    % Loop over all polar angles
    for polarAngleInd = 1:length(polarAngles)
        

        %% CONE MOSAIC

        % Specify polar angle
        cparams.polarAngle   = deg2rad(polarAngles(polarAngleInd));   % Polar angle (radians): 0 is right, pi/2 is superior, pi is left, 3*pi/2 inferior
        
        % Compute x,y position in m of center of retinal patch from ecc and angle
        [x, y] = pol2cart(cparams.polarAngle(1), cparams.eccentricity(1));
        x = x * deg2m;  y = y * deg2m;
        
        cMosaic = coneMosaic('center', [x, y], 'whichEye', whichEye);
        
        % Set the field of view (degrees)
        cMosaic.setSizeToFOV(cparams.cmFOV);
        
        % Add photon noise
        cMosaic.noiseFlag = 'random';
        
        %% DEFINE EYE MOVEMENTS FOR EACH TRIAL
        
        % Compute absorptions for multiple trials
        tSamples         = thisOG.length; % number of time points in optical image sequence

        % Not sure why these have to match, but there is a bug if they don't.
        cMosaic.integrationTime = thisOG.timeStep;

        % ----- EYE MOVEMENTS -----------------------------
        cparams.em        = emCreate;    % eye movements: consider adjusting to
        %   account for cone spacing and for data
        %   from different stimulus conditions
        cparams.em.emFlag = [1 1 1]';    % Include tremor, drift, microsaccades

        % Add eye movements (generate sequence for all trials at once)
        emPaths  = cMosaic.emGenSequence(tSamples, 'nTrials', nTrials, ...
            'em', cparams.em); % path is in terms of cones shifted
        
        %% ABSORPTIONS: Loop through frames to build movie
        
        allAbsorptions = [];
        allCurrents = [];
        
        for trial = 1:nTrials
            
            % Compute
            [absorptions, currents] = cMosaic.compute(thisOG, 'currentFlag', true, 'emPaths', emPaths(trial,:,:)); 
            
            % Concatenate trials
            allAbsorptions = cat(1,allAbsorptions,absorptions);
            allCurrents = cat(1,allCurrents,currents);

        end   
            
        % Store the different OG's and polar angles
        storedAbsorptions{orientationInd, polarAngleInd} = allAbsorptions;
        storedCurrents{orientationInd, polarAngleInd}= allCurrents;
                    
    end
end


%% Save

if saveData
    fname = fullfile(dataPth, sprintf('coneAbsorptions%s.mat', datestr(now, 'YYYY-mm-DD_HH.MM.SS')));
    save(fname,  'storedAbsorptions', 'cparams','sparams');
    
    fname = fullfile(dataPth, sprintf('coneCurrents%s.mat', datestr(now, 'YYYY-mm-DD_HH.MM.SS')));
    save(fname,  'storedCurrents', 'cparams','sparams');
end


%% Classify

labels = [ones(nTrials,1);-1*ones(nTrials,1)];

d = dir(fullfile(dataPth, '*.mat'));
[~,idx] = sort([d.datenum]);
    
load(fullfile(dataPth, d(idx).name));


for polarAngleInd = 1:length(polarAngles)

    thisData = [storedAbsorptions{1,polarAngleInd};storedAbsorptions{2,polarAngleInd}];

    m = fitcsvm(thisData, labels, 'KernelFunction', 'linear');
    cv = crossval(m,'kfold',5);
    rocAreaCones = 1-kfoldLoss(cv);
    
    fprintf('ROC Area for cones: %4.2f\n', rocAreaCones)

    
end
