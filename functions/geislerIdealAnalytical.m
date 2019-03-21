function accuracy = geislerIdealAnalytical(expParams)

% preallocate space
dprimeAnalytic = NaN(size(expParams.contrastLevels));
percentCorrectAnalytic = NaN(size(expParams.contrastLevels));


for c = expParams.contrastLevels
    % Get file name ideal template
    fnameTemplate = sprintf('OGconeOutputs_contrast%1.3f_pa0_eye00_eccen4.50_defocus0.00_noise-none_sf4.00_lms-0.60.30.1.mat',c);
    
    % Load data
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

    % Get all trials with the same label and only take one trial and one time point
    % (since data are all the same across time, without any photon noise or eyemovement or
    % phase shifts)
    alphaMean = template(label==1,:,:,:);
    templateCW = alphaMean(1,:,:,1);
    templateCW = templateCW(:);
    betaMean = template(label==-1,:,:,:);
    templateCCW = betaMean(1,:,:,1);
    templateCCW = templateCCW(:);
    
    if templateCW==templateCCW
        dprimeAnalytic(c==expParams.contrastLevels) = 0;
        percentCorrectAnalytic(c==expParams.contrastLevels) = 0.5;
        
        fprintf('Contrast %1.3f \t d-prime: %2.3f,  percent correct: %2.3f\n', c, dprimeAnalytic(c==expParams.contrastLevels),percentCorrectAnalytic(c==expParams.contrastLevels))
    else
        numerator = sum( (templateCCW-templateCW).*log(templateCCW./templateCW) );
        denominator = sqrt( 0.5* sum( (templateCW+templateCCW) .* (log(templateCCW./templateCW)).^2 ));
        
        dprime = numerator/denominator;
        
        percentCorrectAnalytic(c==expParams.contrastLevels) = normcdf(dprime/2);
        dprimeAnalytic(c==expParams.contrastLevels) = dprime;
        
        fprintf('Contrast %1.3f \t d-prime: %2.3f,  percent correct: %2.3f\n', c, dprimeAnalytic(c==expParams.contrastLevels),percentCorrectAnalytic(c==expParams.contrastLevels))
        
    end
end

saveFolderClassification = fullfile(ogRootPath, 'data', 'classification', expNameTemplate, 'idealtemplate');
fnameClassify = sprintf('ideal_Classify_coneOutputs_contrast%1.3f_pa0_eye00_eccen4.50_defocus0.00_noise-none_sf4.00_lms-0.60.30.1', max(expParams.contrastLevels));
accuracy = percentCorrectAnalytic;
parsave(fullfile(saveFolderClassification, sprintf('%s.mat', fnameClassify)),'accuracy',accuracy);
