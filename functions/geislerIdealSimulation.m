function accuracy = geislerIdealSimulation(expParams, expNameTemplate, expNameData, subFolderName)

% preallocate space
percentCorrectSimulation = NaN(size(expParams.contrastLevels));

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
    
    alphaMean  = template(label==1,:,:,:);
    templateCW  = sum(alphaMean(1,:,:,:),4);
    templateCW  = templateCW(:);
    betaMean    = template(label==-1,:,:,:);
    templateCCW = sum(betaMean(1,:,:,:),4);
    templateCCW = templateCCW(:);
    
    % If one can't distinguish between templates, then don't calculate the
    % rest
    if templateCW==templateCCW
        percentCorrectSimulation(c==expParams.contrastLevels) = 0.5;
        
        % Report percent correct
        fprintf('For contrast %1.4f: percent correct %2.2f%%\n',c, percentCorrectSimulation(c==expParams.contrastLevels)*100);
    else
        
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
        
        % reshape to all trials x [rows x colums x time] for classification
        data = permute(data, [3 1 2 4]);
        
        CWData  = sum(data(label==1,:,:,:),4);
        CWData  = CWData(1:200,:);
        CCWData = sum(data(label==-1,:,:,:),4);
        CCWData = CCWData(1:200,:);
        
        % Compute loglikelihood per trial
        for ii = 1:nTrials
            ZGivenCW(ii) = sum(CWData(ii,:)'  .* log(templateCCW./templateCW));
            ZGivenCCW(ii)= sum(CCWData(ii,:)' .* log(templateCCW./templateCW));
        end
        
        criterion = sum(templateCCW-templateCW);
        
        numberCorrect = sum((ZGivenCW < criterion)) + sum((ZGivenCCW > criterion));
        percentCorrectSimulation(c==expParams.contrastLevels) = (numberCorrect) / (nTrials*2);
        
        % Report percent correct
        fprintf('For contrast %1.4f: percent correct %2.2f%%\n',c, percentCorrectSimulation(c==expParams.contrastLevels)*100);
    end
end

% Save simulation results
saveFolderClassification = fullfile(ogRootPath, 'data', 'classification', expNameData, subFolderName);

if ~exist(saveFolderClassification, 'dir'), mkdir(saveFolderClassification); end

fnameClassify = sprintf('ideal_Classify_coneOutputs_contrast%1.4f_pa0_eye00_eccen4.50_defocus0.00_noise-random_sf4.00_lms-0.60.30.1', max(expParams.contrastLevels));
accuracy = percentCorrectSimulation.*100;
parsave(fullfile(saveFolderClassification, sprintf('%s.mat', fnameClassify)),'accuracy',accuracy);
