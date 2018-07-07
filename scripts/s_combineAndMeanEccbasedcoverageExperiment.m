%% s_combineAndMeanEccbasedcoverageExperiment

%% 0. Set general experiment parameters
expName                  = 'eccbasedcoverage';
expParams                = loadExpParams(expName, false);
[xUnits, colors, labels, M] = loadWeibullPlottingParams(expName);

% Use cone current (flag = true) or cone absorptions (flag = false)
polarAngles = expParams.polarAngle;
FFTflag     = true;

% Where to find data and save figures
dataPth     = fullfile(ogRootPath,'data','PF_data_alias','classification',expName);
figurePth   = fullfile(ogRootPath,'figs', expName, 'average');

% Number of total trials in computational observer model (50 clockwise, 50 counterclockwise)
nTotal      = 100;

% Get nr of conditions
nrEccen    = length(expParams.eccentricities);


% First combine fovea    
d = dir(fullfile(dataPth, 'paddedStim*'));

for ii = 3;%2:size(d,1)
    
    P =[];
    for c = expParams.contrastLevels

        fName   = sprintf('Classify_coneOutputs_contrast%1.3f_pa%d_eye%s_eccen%1.2f_defocus%1.2f_noise-random_sf%1.2f.mat', ...
            c,polarAngles,sprintf('%i',expParams.eyemovement'),expParams.eccentricities(1),expParams.defocusLevels,expParams.spatFreq);


        accuracy = load(fullfile(d(ii).folder, d(ii).name, fName));
        P = [P accuracy.P];

    end

    fName = 'Classify_coneOutputs_contrast0.100_pa0_eye11_eccen0.00_defocus0.00_noise-random_sf4.00_v2.mat';

    save(fullfile(d(ii).folder, d(ii).name, fName), 'P')
end


% Then combine all trials
for ec = 1:nrEccen
        
        fName   = sprintf('Classify_coneOutputs_contrast%1.3f_pa%d_eye%s_eccen%1.2f_defocus%1.2f_noise-random_sf%1.2f.mat', ...
            max(expParams.contrastLevels),polarAngles,sprintf('%i',expParams.eyemovement'),expParams.eccentricities(ec),expParams.defocusLevels,expParams.spatFreq);
        
        d = dir(fullfile(dataPth, 'paddedStim*'));
        
        P =[];
        for ii = 2:size(d,1)
            
            accuracy = load(fullfile(d(ii).folder, d(ii).name, fName));
            accuracy.P = squeeze(accuracy.P);
            if size(accuracy.P,1)<size(accuracy.P,2)
                accuracy.P = accuracy.P';
            end
            
            P = [P accuracy.P];
            
        end
        
        figure; hold all; for ii = 1:size(P,2); plot(expParams.contrastLevels,P(:,ii)); end
        P_SE = std(P,[],2)./sqrt(size(P,2));
        
        P = mean(P,2);
        fName   = sprintf('Classify_coneOutputs_contrast%1.3f_pa%d_eye%s_eccen%1.2f_defocus%1.2f_noise-random_sf%1.2f_AVERAGE.mat', ...
            max(expParams.contrastLevels),polarAngles,sprintf('%i',expParams.eyemovement'),expParams.eccentricities(ec),expParams.defocusLevels,expParams.spatFreq);
        
        fNameSE   = sprintf('Classify_coneOutputs_contrast%1.3f_pa%d_eye%s_eccen%1.2f_defocus%1.2f_noise-random_sf%1.2f_SE.mat', ...
            max(expParams.contrastLevels),polarAngles,sprintf('%i',expParams.eyemovement'),expParams.eccentricities(ec),expParams.defocusLevels,expParams.spatFreq);
        
        
        if ~exist(fullfile(dataPth, 'average'),'dir'), mkdir(fullfile(dataPth, 'average')); end;
        save(fullfile(dataPth, 'average',fName),'P');
        save(fullfile(dataPth, 'average',fNameSE),'P_SE');
    
end