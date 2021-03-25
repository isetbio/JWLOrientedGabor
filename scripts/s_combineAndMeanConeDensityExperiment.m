%% s_combineAndMeanConeDensityExperiment

%% 0. Set general experiment parameters
expName           = 'conedensity';
stimTemplateFlag  = true;
expParams         = loadExpParams(expName, false);
[xUnits, colors, labels, M] = loadWeibullPlottingParams(expName);

polarAngles = expParams.polarAngle;

% Where to find data
if stimTemplateFlag 
    subFolder = 'stimTemplate'; 
    templateName = '_svmEnergy'; % choose from '_svmEnergy' or '_svmLinear'
    baseFolder = '/Volumes/server/Projects/PerformanceFields_RetinaV1Model/';
    dataPth     = fullfile(baseFolder,'data',expName,'classification','absorptions', subFolder);
else 
    templateName = [];
    subFolder = '';
    dataPth = fullfile(ogRootPath,'data','classification',expName);
end

% Where to save data and figures
figurePth   = fullfile(ogRootPath,'figs', expName, ['average' templateName]);
dataSavePth = fullfile(dataPth, ['average' templateName]);


% Number of total trials in computational observer model (50 clockwise, 50 counterclockwise)
nTotal      = 100;

% Get nr of conditions
nrEccen     = length(expParams.eccentricities);

% Get cone ratio
lmsRatio    = expParams.cparams.spatialDensity;


% First combine fovea
for ec = 1:nrEccen
    
    fName   = sprintf('Classify_coneOutputs_contrast%1.4f_pa%d_eye%s_eccen%1.2f_defocus%1.2f_noise-random_sf%1.2f_lms-%1.1f%1.1f%1.1f.mat', ...
        max(expParams.contrastLevels),polarAngles,sprintf('%i',expParams.eyemovement'),expParams.eccentricities(ec),expParams.defocusLevels,expParams.spatFreq, lmsRatio(2),lmsRatio(3),lmsRatio(4));
    
    d = dir(fullfile(dataPth, 'run*'));
    
    P =[];
    for ii = 1:size(d,1)
        
        tmp = load(fullfile(d(ii).folder, d(ii).name, fName));
        
        if stimTemplateFlag
            fn = fieldnames(tmp);
            P_svm = tmp.(fn{strcmpi(fn,['P' templateName])});
        else
            P_svm = squeeze(tmp.accuracy);
        end
        
        if size(P_svm,1)<size(P_svm,2)
            P_svm = P_svm';
        end
        
        P = [P P_svm];
        
    end
    
    figure; hold all; for ii = 1:size(P,2); plot(expParams.contrastLevels,P(:,ii)); end
    P_SE = std(P,[],2)./sqrt(size(P,2));
    
    P_AVG = mean(P,2);
    fName   = sprintf('Classify_coneOutputs_contrast%1.3f_pa%d_eye%s_eccen%1.2f_defocus%1.2f_noise-random_sf%1.2f_lms-%1.1f%1.1f%1.1f_AVERAGE.mat', ...
        max(expParams.contrastLevels),polarAngles,sprintf('%i',expParams.eyemovement'),expParams.eccentricities(ec),expParams.defocusLevels,expParams.spatFreq, lmsRatio(2),lmsRatio(3),lmsRatio(4));
    
    fNameSE   = sprintf('Classify_coneOutputs_contrast%1.3f_pa%d_eye%s_eccen%1.2f_defocus%1.2f_noise-random_sf%1.2f_lms-%1.1f%1.1f%1.1f_SE.mat', ...
        max(expParams.contrastLevels),polarAngles,sprintf('%i',expParams.eyemovement'),expParams.eccentricities(ec),expParams.defocusLevels,expParams.spatFreq, lmsRatio(2),lmsRatio(3),lmsRatio(4));
    
    
    if ~exist(dataSavePth,'dir'), mkdir(dataSavePth); end;
    save(fullfile(dataSavePth, fName),'P_AVG');
    save(fullfile(dataSavePth, fNameSE),'P_SE');
    
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
if ~exist(fullfile(baseFolder,'data',expName,'thresholds',subFolder, ['stimTemplate' templateName]), 'dir'); 
    mkdir(fullfile(baseFolder,'data',expName,'thresholds',subFolder, ['stimTemplate' templateName])); end
save(fullfile(baseFolder,'data',expName,'thresholds',subFolder, ['stimTemplate' templateName], fNameSEThresh),'varThresh', 'ctrthresh');



