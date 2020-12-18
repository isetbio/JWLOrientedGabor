%% s_combineAndMeanConeDensityExperiment

%% 0. Set general experiment parameters
expName                  = 'conedensity';
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
nrEccen     = length(expParams.eccentricities);

% Get cone ratio
lmsRatio    = expParams.cparams.spatialDensity;


% First combine fovea
for ec = 1:nrEccen
    
    fName   = sprintf('Classify_coneOutputs_contrast%1.3f_pa%d_eye%s_eccen%1.2f_defocus%1.2f_noise-random_sf%1.2f_lms-%1.1f%1.1f%1.1f.mat', ...
        max(expParams.contrastLevels),polarAngles,sprintf('%i',expParams.eyemovement'),expParams.eccentricities(ec),expParams.defocusLevels,expParams.spatFreq, lmsRatio(2),lmsRatio(3),lmsRatio(4));
    
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
    
    P_AVG = mean(P,2);
    fName   = sprintf('Classify_coneOutputs_contrast%1.3f_pa%d_eye%s_eccen%1.2f_defocus%1.2f_noise-random_sf%1.2f_lms-%1.1f%1.1f%1.1f_AVERAGE.mat', ...
        max(expParams.contrastLevels),polarAngles,sprintf('%i',expParams.eyemovement'),expParams.eccentricities(ec),expParams.defocusLevels,expParams.spatFreq, lmsRatio(2),lmsRatio(3),lmsRatio(4));
    
    fNameSE   = sprintf('Classify_coneOutputs_contrast%1.3f_pa%d_eye%s_eccen%1.2f_defocus%1.2f_noise-random_sf%1.2f_lms-%1.1f%1.1f%1.1f_SE.mat', ...
        max(expParams.contrastLevels),polarAngles,sprintf('%i',expParams.eyemovement'),expParams.eccentricities(ec),expParams.defocusLevels,expParams.spatFreq, lmsRatio(2),lmsRatio(3),lmsRatio(4));
    
    
    %         if ~exist(fullfile(dataPth, 'average'),'dir'), mkdir(fullfile(dataPth, 'average')); end;
    %         save(fullfile(dataPth, 'average',fName),'P_AVG');
    %         save(fullfile(dataPth, 'average',fNameSE),'P_SE');
    
    %% Bootstrap runs with replacement,
    nboot = 1000;
    bootData = bootstrp(nboot, @mean, P');
    
    % fit Weibull to each mean and get threshold
    [xUnits, ~, ~, ~, ~] = loadWeibullPlottingParams(expName);
    
    % Prepare fit variables
    fit = [];
    
    % Set inital slope, threshold for first stage fitting
    fit.init   = [2, 0.01]; % slope, threshold at ~80%
    fit.thresh = 0.75;
    
    
    for ii = 1:nboot
        %% 4. Fit Weibull
        % Make a Weibull function first with contrast levels and then search for
        % the best fit with the classifier data
        ctrvar = fminsearch(@(x) ogFitWeibull(x, expParams.contrastLevels, bootData(ii,:)', nTotal), fit.init);
        
        % Then fit a Weibull function again, but now with the best fit parameters
        % from the previous step.
        ctrpred = ogWeibull(ctrvar, xUnits);
        
        %% 5. Find contrast threshold
        ctrthresh(ec,ii) = ctrvar(2);
    end
end

varThresh = std(ctrthresh,[],2);
fNameSEThresh = sprintf('varThresh_coneResponse_absorptionrate_%d_conedensity.mat', ec);
baseFolder = '/Volumes/server/Projects/PerformanceFields_RetinaV1Model';
if ~exist(fullfile(baseFolder,'data',expName,'thresholds'), 'dir'); mkdir(fullfile(baseFolder,'data',expName,'thresholds')); end
save(fullfile(baseFolder,'data',expName,'thresholds',fNameSEThresh),'varThresh', 'ctrthresh');



