function accuracy = svmTemplateSimulation(expParams, expNameTemplate, expNameData, subFolderName)


% Preallocate space
percentCorrectSVMDiffTemplate  = NaN(size(expParams.contrastLevels));

for c = expParams.contrastLevels
    
    % Get file name ideal template
    fnameTemplate = sprintf('OGconeOutputs_contrast%1.4f_pa0_eye00_eccen4.50_defocus0.00_noise-none_sf4.00_lms-0.60.30.1.mat',c);
    
    % Load template
    template = load(fullfile(ogRootPath, 'data', expNameTemplate, 'idealtemplate', fnameTemplate));
    template = template.absorptions;
    
    % Get the trials and samples (should be the data for all data sets though
    nStimuli = size(template,5);
    nTrials  = size(template,1) * nStimuli/2;
    tSamples = size(template,4);
    nrows    = size(template,2);
    ncols    = size(template,3);
    
    %   permute to trials x stimuli x rows x cols x time points
    template = permute(template, [1 5 2:4]);
    
    %   reshape to (trials x stimuli) x rows x cols x time points
    template = reshape(template, [], nrows, ncols, tSamples);
    
    % Label clockwise and counterclockwise trials
    label = [ones(nTrials, 1); -ones(nTrials, 1)]; % First set is CW, second set is CCW
    
    alphaMean = template(label==1,:,:,:);
    templateCW = alphaMean(1,:,:,1);
    betaMean = template(label==-1,:,:,:);
    templateCCW = betaMean(1,:,:,1);
    
    % Get data
    fnameData = sprintf('OGconeOutputs_contrast%1.4f_pa0_eye00_eccen4.50_defocus0.00_noise-random_sf4.00_lms-0.60.30.1.mat', c);
    data = load(fullfile(ogRootPath, 'data', expNameData, subFolderName, fnameData));
    data = data.absorptions;

    %   permute to trials x stimuli x rows x cols x time points
    data = permute(data, [1 5 2:4]);

    %   reshape to (trials x stimuli) x rows x cols x time points
    data = reshape(data, [], nrows, ncols, tSamples);

    % permute to rows x cols x (trials x stimuli) x time points
    data  = permute(data, [2 3 1 4]);
    
     % take dotproduct with difference in noiseless CW/CCW template
    diffTemplate = templateCCW - templateCW;
    
    data2 = [];
    for ii = 1:nTrials*2
        for jj = 1:tSamples
            data2(:,:,ii,jj) = data(:,:,ii, jj)*squeeze(diffTemplate);
        end
    end

    % reshape to all trials x [rows x colums x time] for classification
    data2 = permute(data2, [3 1 2 4]);
    
    % permute the trial order within each of the two classes
    idx = [randperm(nTrials) randperm(nTrials)+nTrials];
    data2 = data2(idx, :);

    label = [ones(nTrials, 1); -ones(nTrials, 1)];

    % Fit the SVM model.
    cvmdl = fitcsvm(data2, label, 'Standardize', false, 'KernelFunction', 'linear', 'kFold', 10);

    % predict the data not in the training set.
    classLoss = kfoldLoss(cvmdl);

    % Get percent accuracy
    P = (1-classLoss);
    percentCorrectSVMDiffTemplate(c==expParams.contrastLevels) = P;
        
    % Report percent correct
    fprintf('For contrast %1.4f: percent correct %2.2f%%\n',c, percentCorrectSVMDiffTemplate(c==expParams.contrastLevels)*100);
    
end

% Save simulation results
saveFolderClassification = fullfile(ogRootPath, 'data', 'classification', expNameData, subFolderName);
fnameClassify = sprintf('svmtemplate_Classify_coneOutputs_contrast%1.3f_pa0_eye00_eccen4.50_defocus0.00_noise-random_sf4.00_lms-0.60.30.1', max(expParams.contrastLevels));
accuracy = percentCorrectSVMDiffTemplate.*100;
parsave(fullfile(saveFolderClassification, sprintf('%s.mat', fnameClassify)),'accuracy',accuracy);

