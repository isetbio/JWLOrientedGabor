function accuracy = geislerIdealAnalytical(varargin)
% Function to compute closed form solution of ideal observer model as in Geisler (1984).
%
% accuracy = geislerIdealAnalytical(expParams)
%
% INPUTS:
%   expParams       : struct with experiment parameters (from loadExpParams)

%% Derive variables from parsed inputs
p = inputParser;
p.KeepUnmatched = true;  
p.addParameter('expParams', @isstruct);
p.parse(varargin{:});

expParams        = p.Results.expParams; 

% preallocate space
dprimeAnalytic = NaN(size(expParams.contrastLevels));
percentCorrectAnalytic = NaN(size(expParams.contrastLevels));

% define data alias
dataPathAlias = 'PF_data_alias';
dataPathFull  = fullfile(ogRootPath, 'data', dataPathAlias, 'coneabsorptions',expParams.name, 'idealtemplate');
dataPathFull  = fullfile(ogRootPath, 'data', expParams.name, 'onlyL');

saveFolderClassification = fullfile(ogRootPath, 'data', 'classification', expParams.name, 'onlyL');

%% Get d' and percent correct for every contrast level
for c = expParams.contrastLevels
    % Get file name ideal template
%     fnameTemplate = sprintf('OGconeOutputs_contrast%1.4f_pa0_eye00_eccen4.50_defocus0.00_noise-none_sf4.00_lms-0.60.30.1.mat',c);
    fnameTemplate = sprintf('OGconeOutputs_contrast%1.4f_pa0_eye00_eccen4.50_defocus0.00_noise-none_sf4.00_lms-1.00.00.0.mat',c);

    % Load data
    template = load(fullfile(dataPathFull, fnameTemplate));
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

    % Since all trials are without photon noise, phase shifts or eye
    % movements, they have the same cone absorptions. However, we need to
    % sum across all time points to have a fair comparison to the SVM
    % results.
    alphaMean  = template(label==1,:,:,1:28);
    templateCW  = sum(alphaMean(1,:,:,:),4);
    templateCW  = templateCW(:);
    betaMean    = template(label==-1,:,:,1:28);
    templateCCW = sum(betaMean(1,:,:,:),4);
    templateCCW = templateCCW(:);
    
    if templateCW==templateCCW
        dprimeAnalytic(c==expParams.contrastLevels) = 0;
        percentCorrectAnalytic(c==expParams.contrastLevels) = 0.5;
        
        fprintf('Contrast %1.4f \t d-prime: %2.3f,  percent correct: %2.3f\n', c, dprimeAnalytic(c==expParams.contrastLevels),percentCorrectAnalytic(c==expParams.contrastLevels))
    else
        numerator = sum( (templateCCW-templateCW).*log(templateCCW./templateCW) );
        denominator = sqrt( 0.5* sum( (templateCW+templateCCW) .* (log(templateCCW./templateCW)).^2 ));
        
        dprime = numerator/denominator;
        
        percentCorrectAnalytic(c==expParams.contrastLevels) = normcdf(dprime/2);
        dprimeAnalytic(c==expParams.contrastLevels) = dprime;
        
        fprintf('Contrast %1.4f \t d-prime: %2.3f,  percent correct: %2.3f\n', c, dprimeAnalytic(c==expParams.contrastLevels),percentCorrectAnalytic(c==expParams.contrastLevels))
        
    end
end

%% Save data
if ~exist(saveFolderClassification, 'dir'), mkdir(saveFolderClassification); end
fnameClassify = sprintf('ideal_Classify_coneOutputs_contrast%1.4f_pa0_eye00_eccen4.50_defocus0.00_noise-none_sf4.00_lms-1.00.00.0', max(expParams.contrastLevels));
% fnameClassify = sprintf('ideal_Classify_coneOutputs_contrast%1.4f_pa0_eye00_eccen4.50_defocus0.00_noise-none_sf4.00_lms-0.60.30.1', max(expParams.contrastLevels));
accuracy = percentCorrectAnalytic.*100;
parsave(fullfile(saveFolderClassification, sprintf('%s.mat', fnameClassify)),'accuracy',accuracy, 'expParams', expParams);
