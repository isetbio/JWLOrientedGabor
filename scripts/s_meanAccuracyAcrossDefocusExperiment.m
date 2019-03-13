% Compute mean over defocus percent correct as a function of contrast

%% 0. Set general experiment parameters
expName                  = 'defocus';
expParams                = loadExpParams(expName, false);
[xUnits, colors, labels, M] = loadWeibullPlottingParams(expName);

% Use cone current (flag = true) or cone absorptions (flag = false)
polarAngles = expParams.polarAngle;
FFTflag     = true;

% Where to find data and save figures
dataPth     = fullfile(ogRootPath,'data','PF_data_alias','classification',expName, 'HPC');
figurePth   = fullfile(ogRootPath,'figs', expName, 'averageHPC');

% Number of total trials in computational observer model (50 clockwise, 50 counterclockwise)
nTotal      = 100;%expParams.nTrials*4;

% Get nr of conditions
nrDefocusLevels    = length(expParams.defocusLevels);
lmsRatio = expParams.cparams.spatialDensity;

for df = 1:nrDefocusLevels
    
    fName   = sprintf('Classify_coneOutputs_contrast%1.3f_pa%d_eye%s_eccen%1.2f_defocus%1.2f_noise-random_sf%1.2f_lms-%1.1f%1.1f%1.1f.mat', ...
        max(expParams.contrastLevels),polarAngles,sprintf('%i',expParams.eyemovement'),expParams.eccentricities,expParams.defocusLevels(df),expParams.spatFreq,lmsRatio(2),lmsRatio(3),lmsRatio(4));    
    
    d = dir(fullfile(dataPth, 'run*'));
    
    P =[];
    for ii = 1:size(d,1)
        
        accuracy = load(fullfile(d(ii).folder, d(ii).name, fName));
        accuracy.P = squeeze(accuracy.accuracy);
        if size(accuracy.P,1)<size(accuracy.P,2)
            accuracy.P = accuracy.P';
        end
        
        P = [P accuracy.P];
        
    end
    
     figure; hold all; for ii = 1:size(P,2); plot(expParams.contrastLevels,P(:,ii)); end
    P_SE = std(P,[],2)./sqrt(size(P,2));
    
    P = mean(P,2);
    fName   = sprintf('Classify_coneOutputs_contrast%1.3f_pa%d_eye%s_eccen%1.2f_defocus%1.2f_noise-random_sf%1.2f_lms-%1.1f%1.1f%1.1f_AVERAGE.mat', ...
        max(expParams.contrastLevels),polarAngles,sprintf('%i',expParams.eyemovement'),expParams.eccentricities,expParams.defocusLevels(df),expParams.spatFreq,lmsRatio(2),lmsRatio(3),lmsRatio(4));
    
    fNameSE   = sprintf('Classify_coneOutputs_contrast%1.3f_pa%d_eye%s_eccen%1.2f_defocus%1.2f_noise-random_sf%1.2f_lms-%1.1f%1.1f%1.1f_SE.mat', ...
        max(expParams.contrastLevels),polarAngles,sprintf('%i',expParams.eyemovement'),expParams.eccentricities,expParams.defocusLevels(df),expParams.spatFreq,lmsRatio(2),lmsRatio(3),lmsRatio(4));
    
    
    if ~exist(fullfile(dataPth, 'average'),'dir'), mkdir(fullfile(dataPth, 'average')); end;
    save(fullfile(dataPth, 'average',fName),'P');
    save(fullfile(dataPth, 'average',fNameSE),'P_SE');
    
end