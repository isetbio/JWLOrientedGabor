% Compute mean over eccbasedcoverage percent correct

%% 0. Set general experiment parameters
expName                  = 'eccbasedcoverage';
expParams                = loadExpParams(expName, false);
[xUnits, colors, labels, M] = loadWeibullPlottingParams(expName);

% Use cone current (flag = true) or cone absorptions (flag = false)
currentFlag = false;
polarAngles = expParams.polarAngle;
FFTflag     = true;
saveFig     = true;

% Where to find data and save figures
dataPth     = fullfile(ogRootPath,'data','classification',expName);
figurePth   = fullfile(ogRootPath,'figs', expName, 'average');

% Number of total trials in computational observer model (50 clockwise, 50 counterclockwise)
nTotal      = 100;%expParams.nTrials*4;

% Get nr of conditions
nrEyemovTypes    = size(expParams.eyemovement,2);
nrEccen          = length(expParams.eccentricities);
nrSpatFreq       = length(expParams.spatFreq);
nrDefocusLevels  = length(expParams.defocusLevels);

for eccen = 1:nrEccen
    
    fName   = sprintf('Classify_coneOutputs_contrast%1.3f_pa%d_eye%s0_eccen%1.2f_defocus%1.2f_noise-random_sf%1.2f.mat', ...
        max(expParams.contrastLevels),polarAngles,sprintf('%i',expParams.eyemovement(:,em)'),expParams.eccentricities(eccen),expParams.defocusLevels(df),expParams.spatFreq);
    if currentFlag; fName = ['current_' fName]; end;
    
    d = dir(fullfile(dataPth, '100trials_*'));
    
    P =[];
    for ii = 1:size(d,1)
        try
        accuracy = load(fullfile(d(ii).folder, d(ii).name, fName));
        accuracy.P = squeeze(accuracy.P);
        if size(accuracy.P,1)<size(accuracy.P,2)
            accuracy.P = accuracy.P';
        end
        
        P = [P accuracy.P];
        end
    end
    
    P = mean(P,2);
    fName   = sprintf('Classify_coneOutputs_contrast%1.3f_pa%d_eye%s0_eccen%1.2f_defocus%1.2f_noise-random_sf%1.2f_AVERAGE2.mat', ...
        max(expParams.contrastLevels),polarAngles,sprintf('%i',expParams.eyemovement(:,em)'),expParams.eccentricities(eccen),expParams.defocusLevels(df),expParams.spatFreq);
    save(fullfile(dataPth, 'average',fName),'P');
    
end