%% s_ogRGC_Classify_HPC_eyemov

% Script with first attempt to classify oriented gabors simulated at
% 7 different contrast levels, 4 polar angles.

%% Classify

% Load experiment parameters
expName = 'eyemov';
subFolderName_toSave = 'paddedStim10';
subFolderName_toLoad = 'paddedStim10';
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

parfor em = 1:nrEyemovTypes
    
    P = nan(nrContrasts,1);    
    for c = expParams.contrastLevels
        
        % Load dataset
        fname = sprintf(...
            'OGconeOutputs_contrast%1.3f_pa%d_eye%s_eccen%1.2f_defocus%1.2f_noise-random_sf%1.2f.mat',...
            c,expParams.polarAngle,sprintf('%i',expParams.eyemovement(:,em)), expParams.eccentricities(1), expParams.defocusLevels(1), expParams.spatFreq(1));
        
        if currentFlag
            fname = ['current_' fname];
        end
        
        pth = fullfile(ogRootPath, 'data', expName, subFolderName_toLoad, fname);
        if ~exist(pth, 'file'), error('The file %s is not found', fname); end
        
        tmp = load(pth);
        
        if currentFlag
            data = getfield(tmp,'current');
        else
            data = getfield(tmp,'absorptions');
        end
        
        
        fprintf('Loading and classifying %s\n', fname);
        % Get the trials and samples (should be the data for all data sets though
        nStimuli = size(data,5);
        nTrials  = size(data,1) * nStimuli/2;
        tSamples = size(data,4);
        nrows    = size(data,2);
        ncols    = size(data,3);
        % absorptions is trials x rows x cols x time points x stimuli
        
        %   permute to trials x stimuli x rows x cols x time points
        data = permute(data, [1 5 2:4]);
        
        %   reshape to (trials x stimuli) x rows x cols x time points
        data = reshape(data, [], nrows, ncols, tSamples);
        
        % permute to rows x cols x (trials x stimuli) x time points
        data  = permute(data, [2 3 1 4]);
        
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
        cvmdl = fitcsvm(data, label, 'Standardize', true, 'KernelFunction', 'linear', 'kFold', 10);
        
        % predict the data not in the training set.
        classLoss = kfoldLoss(cvmdl);

        P(c==expParams.contrastLevels) = (1-classLoss) * 100;

    end
    
    disp(P);
    
    % Save classifier accuracy
    fname = sprintf(...
        'Classify_coneOutputs_contrast%1.3f_pa%d_eye%s_eccen%1.2f_defocus%1.2f_noise-random_sf%1.2f',...
        c, expParams.polarAngle,sprintf('%i',expParams.eyemovement(:,em)), expParams.eccentricities(1), expParams.defocusLevels(1), expParams.spatFreq(1));
    if currentFlag; fname = ['current_' fname]; end
    parsave(fullfile(savePth, sprintf('%s.mat', fname)),'P',P)
    
end


return









