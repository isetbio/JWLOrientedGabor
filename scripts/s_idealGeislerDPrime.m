% s_idealGeislerDPrime

% Load parameters
expNameTemplate = 'idealobserver';
expNameData     = 'defaultnophaseshift';
expParams       = loadExpParams(expNameTemplate, false);

% save figures?
saveFigs = false;
figure(1); clf; hold all;

% preallocate space
dprimeAnalytic = NaN(size(expParams.contrastLevels));
dprimeSimulation = NaN(size(expParams.contrastLevels));
percentCorrectAnalytic = NaN(size(expParams.contrastLevels));
percentCorrectSimulation = NaN(size(expParams.contrastLevels));
percentCorrectSVMDiffTemplate  = NaN(size(expParams.contrastLevels));

 %% 1. Analytical solution
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
accuracy = percentCorrectAnalytic.*100;
parsave(fullfile(saveFolderClassification, sprintf('%s.mat', fnameClassify)),'accuracy',accuracy);

logzero = 4e-4;
figure(1); plot(expParams.contrastLevels(2:end), percentCorrectAnalytic(2:end), 'ko-'); hold on;
plot(logzero, percentCorrectAnalytic(1), 'ko')
set(gca, 'XScale', 'log', 'YLim', [.4, 1])

%% 2. Simulation of analytic solution

% Load data parameters
expParams = loadExpParams(expNameData, false);
subFolderName = 'run1';

for c = expParams.contrastLevels
    
    % Get file name ideal template
    fnameTemplate = sprintf('OGconeOutputs_contrast%1.3f_pa0_eye00_eccen4.50_defocus0.00_noise-none_sf4.00_lms-0.60.30.1.mat',c);
    
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
    templateCW = templateCW(:);
    betaMean = template(label==-1,:,:,:);
    templateCCW = betaMean(1,:,:,1);
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
        
        CWData = data(label==1,:,:,1);
        CWData = CWData(1:200,:);
        CCWData = data(label==-1,:,:,1);
        CCWData = CCWData(1:200,:);
        
        % Compute loglikelihood per trial
        for ii = 1:nTrials
            ZGivenCW(ii) = sum(CWData(ii,:)'  .* log(templateCCW./templateCW));
            ZGivenCCW(ii)= sum(CCWData(ii,:)' .* log(templateCCW./templateCW));
        end
        
        criterion = sum(templateCCW-templateCW);
        
        numberCorrect = sum((ZGivenCW < criterion)) + sum((ZGivenCCW > criterion));
        percentCorrectSimulation(c==expParams.contrastLevels) = (numberCorrect) / (nTrials*2) ;
        
        % Report percent correct
        fprintf('For contrast %1.4f: percent correct %2.2f%%\n',c, percentCorrectSimulation(c==expParams.contrastLevels)*100);
    end
end

% Plot simulation results
figure(1); plot(expParams.contrastLevels(2:end), percentCorrectSimulation(2:end), 'ro-');
plot(logzero, percentCorrectSimulation(1), 'ro')
set(gca, 'XScale', 'log', 'YLim', [.4, 1])

% Save simulation results
saveFolderClassification = fullfile(ogRootPath, 'data', 'classification', expNameData, subFolderName);
fnameClassify = sprintf('ideal_Classify_coneOutputs_contrast%1.3f_pa0_eye00_eccen4.50_defocus0.00_noise-random_sf4.00_lms-0.60.30.1', max(expParams.contrastLevels));
accuracy = percentCorrectSimulation.*100;
parsave(fullfile(saveFolderClassification, sprintf('%s.mat', fnameClassify)),'accuracy',accuracy);

%% 3. SVM with data sample match to template

% Load data parameters
expParams = loadExpParams(expNameData, false);
subFolderName = 'run1';

for c = expParams.contrastLevels
    
    % Get file name ideal template
    fnameTemplate = sprintf('OGconeOutputs_contrast%1.3f_pa0_eye00_eccen4.50_defocus0.00_noise-none_sf4.00_lms-0.60.30.1.mat',c);
    
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
    cvmdl = fitcsvm(data2, label, 'Standardize', true, 'KernelFunction', 'linear', 'kFold', 10);

    % predict the data not in the training set.
    classLoss = kfoldLoss(cvmdl);

    % Get percent accuracy
    P = (1-classLoss) * 100;
    percentCorrectSVMDiffTemplate(c==expParams.contrastLevels) = P;
        
    % Report percent correct
    fprintf('For contrast %1.4f: percent correct %2.2f%%\n',c, percentCorrectSVMDiffTemplate(c==expParams.contrastLevels));
    
end

% Plot simulation results
figure(1); plot(expParams.contrastLevels(2:end), percentCorrectSVMDiffTemplate(2:end), 'ro-');
plot(logzero, percentCorrectSVMDiffTemplate(1), 'ro')
set(gca, 'XScale', 'log', 'YLim', [.4, 1])

% Save simulation results
saveFolderClassification = fullfile(ogRootPath, 'data', 'classification', expNameData, subFolderName);
fnameClassify = sprintf('svmtemplate_Classify_coneOutputs_contrast%1.3f_pa0_eye00_eccen4.50_defocus0.00_noise-random_sf4.00_lms-0.60.30.1', max(expParams.contrastLevels));
accuracy = percentCorrectSVMDiffTemplate;
parsave(fullfile(saveFolderClassification, sprintf('%s.mat', fnameClassify)),'accuracy',accuracy);






% OLD CODE 
%     % Compute mean and varianec across trials
%     meanZGivenCCW = mean(ZGivenCCW);
%     meanZGivenCW = mean(ZGivenCW);
%     varZGivenCCW = var(ZGivenCCW);
%     varZGivenCW = var(ZGivenCW);
%
%     % Compute d-prime
%     dprime = (meanZGivenCCW-meanZGivenCW) / sqrt(0.5*(varZGivenCCW+varZGivenCW));
%
%
%     % Plot loglikelihood distributions of the two stimulus classes
%     figure(2); clf; set(gcf, 'Color', 'w', 'Position', [508, 820, 1052, 518]);
%     subplot(1,2,1); hold all;
%     [nA,xA] = hist(ZGivenCW,40);
%     bar(xA,nA, 'EdgeColor', 'k', 'FaceColor', 'w')
%     predA = length(ZGivenCW)*normpdf(xA,meanZGivenCW,sqrt(varZGivenCW))*(xA(2)-xA(1));
%     plot(xA,predA,'g','LineWidth',4);
%     set(gca, 'TickDir', 'out', 'FontSize', 14);
%
%     [nB,xB] = hist(ZGivenCCW,40);
%     bar(xB,nB, 'EdgeColor', 'k', 'FaceColor', 'k');
%     predB = length(ZGivenCCW)*normpdf(xB,meanZGivenCCW,sqrt(varZGivenCCW))*(xB(2)-xB(1));
%     plot(xB,predB,'r','LineWidth',4);
%     title(sprintf('Probability of absorptions for Alpha vs Beta stim, contrast: %1.4f',c))
%     legend({'Alpha samples', 'Alpha fit', 'Beta samples', 'Beta fit'}, 'FontSize', 14, 'Location', 'best'); legend boxoff;
%
%     % Get ROC curve from Hit and False Alarm rates and sample along many criteria
%     criteria = linspace(min([ZGivenCW,ZGivenCCW]),max([ZGivenCW,ZGivenCCW]),1000);
%     for cc = 1:length(criteria)
%         HitRate(cc) = length(find(ZGivenCCW > criteria(cc))) / length(ZGivenCCW);
%         FARate(cc) = length(find(ZGivenCW > criteria(cc))) / length(ZGivenCCW);
%     end
%
%     % Get the percent correct
%     percentCorrectSim2(c==expParams.contrastLevels) = -trapz([1 FARate 0],[1 HitRate 0]);
%     subplot(1,2,2);
%     plot([1 FARate 0],[1 HitRate 0],'r','LineWidth',4);
%     xlabel('False Alarm Rate');
%     ylabel('Hit Rate');
%     title(sprintf('ROC Curve, d-prime: %2.2f, percent corrent: %2.2f',dprime,percentCorrectSim2(c==expParams.contrastLevels)*100));
%     axis('square'); axis([0 1 0 1]); set(gca, 'TickDir', 'out', 'FontSize', 14);
%     drawnow;
%     if saveFigs
%         savefig(fullfile(figurePth,sprintf('Geisler_Ideal_dprime_%s_c%1.4f.fig',expNameTemplate,c)))
%         hgexport(gcf,fullfile(figurePth,sprintf('Geisler_Ideal_dprime_%s_c%1.4f.eps',expNameTemplate,c)))
%     end