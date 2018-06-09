%% s_ogRGC_Classify_HPC_eccbasedcoverage_fovea_trainHighContrast

% Script with first attempt to classify oriented gabors simulated at
% 7 different contrast levels, 4 polar angles.

%% Classify

% Load experiment parameters
expName = 'eccbasedcoverage';
subFolderName_toSave = '100trials_trainHighContrast';
subFolderName_toLoad = 'paddedStim';
expParams = loadExpParams(expName, false);

% Compute accuracy for cone current as well
currentFlag    = false;

% Compute accuracy on fft component of cone absorptions
fftFlag        = true;

% Predefine matrix for predictions
nrContrasts      = length(expParams.contrastLevels);
nrEyemovTypes    = size(expParams.eyemovement,2);
nrEccen          = length(expParams.eccentricities);
nrSpatFreq       = length(expParams.spatFreq);
nrDefocusLevels  = length(expParams.defocusLevels);


savePth = fullfile(ogRootPath, 'data', 'classification', expName, subFolderName_toSave);
if ~exist('savePth', 'dir'); mkdir(savePth); end;

% Init figure
figure; clf; set(gcf,'Color','w'); hold all;
set(gca, 'XScale','log', 'XLim', [.005 max(expParams.contrastLevels)], 'XTick', [1:7, 10:10:100]/100, ...
    'YLim', [40 100], 'TickDir','out','TickLength',[.015 .015]);
ylabel('Classifier Accuracy')
xlabel('Contrast level (Michelson)')

df = 1;
em = 1;
sf = expParams.spatFreq(1);
eccen = 1;

% Train on high contrast
[data, nTrials] = loadAndPermuteData(expParams, nrContrasts(end), em, eccen, df, sf, currentFlag, subFolderName_toLoad);

% Compute fourier transform the cone array outputs
if fftFlag; data  = abs(fft2(data)); end

% reshape to all trials x [rows x colums x time] for classification
data = permute(data, [3 1 2 4]);
data = reshape(data, nTrials*2, []);

% permute the trial order within each of the two classes
idx = [randperm(nTrials) randperm(nTrials)+nTrials];

data = data(idx, :);

label = [ones(nTrials, 1); -ones(nTrials, 1)];

% Fit the SVM model.
mdl = fitcsvm(data, label, 'Standardize', true, 'KernelFunction', 'linear');

parfor c = 1:nrContrasts
    
    [data, nTrials] = loadAndPermuteData(expParams, c, em, eccen, df, sf, currentFlag, subFolderName_toLoad);
    
    % Compute fourier transform the cone array outputs
    if fftFlag; data  = abs(fft2(data)); end
    
    % reshape to all trials x [rows x colums x time] for classification
    data = permute(data, [3 1 2 4]);
    data = reshape(data, nTrials*2, []);
    
    % permute the trial order within each of the two classes
    idx = [randperm(nTrials) randperm(nTrials)+nTrials];
    
    data = data(idx, :);
    
    label = [ones(nTrials, 1); -ones(nTrials, 1)];
    
    [predictedLabel,score] = predict(mdl,data);
    
    % Fit the SVM model.
    %                     cvmdl = fitcsvm(data, label, 'Standardize', true, 'KernelFunction', 'linear', 'kFold', 10);
    
    %             cvmdl = crossval(mdl);
    
    % predict the data not in the training set.
    %                     classLoss = kfoldLoss(cvmdl);
    
    % Different type of linear classifier (faster, but less
    % accurate)
    %             mdl = fitclinear(data', label,  'KFold', 10, 'ObservationsIn', 'columns');
    %             classLoss = kfoldLoss(mdl);
    
    %                     P(c) = (1-classLoss) * 100;
    
    P = (sum(label==predictedLabel)/length(label))*100;
    
    disp(P);
    
    % Save classifier accuracy
    fname = sprintf(...
        'Classify_coneOutputs_contrast%1.3f_pa%d_eye%s_eccen%1.2f_defocus%1.2f_noise-random_sf%1.2f',...
        c, expParams.polarAngle,sprintf('%i',expParams.eyemovement(:,em)), expParams.eccentricities(eccen), expParams.defocusLevels(df), sf);
    if currentFlag; fname = ['current_' fname]; end
    parsave(fullfile(savePth, sprintf('%s.mat', fname)),'P',P)
    
end

return











