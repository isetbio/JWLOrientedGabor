function percentCorrect = getIdealObserverAccuracy(data, expName, subFolder, baseFolder, cone2RGCRatio)
% Function to train and test linear SVM classifier on cone data with cross-validation:
%       P = getClassifierAccuracy(data)
%
% INPUTS: 
%   data        : 5 dimensional array (with trials x rows x cols x time
%                                       samples x stimuli) 
% OUTPUTS:
%   P           : classifier accuracy of computational observer model in
%                   percent correct for given absorption dataset

%% 1. Transform and reshape absorption data:

% Get dimensions of data
fprintf('(%s): Loading and classifying\n', mfilename);

% Load parameters
% expName = 'idealobserver';
% subFolder = 'idealtemplate';
expParams = loadExpParams(expName, false);


%% Analytical solution
for c = expParams.contrastLevels
%     fnameTemplate = sprintf('OGconeOutputs_contrast%1.4f_pa0_eye00_eccen4.50_defocus0.00_noise-none_sf4.00_lms-1.00.00.0.mat',c);
%     template = load(fullfile(baseFolder, 'data', expName, subFolder, fnameTemplate));
%     template = template.absorptions;   

%     fnameTemplate = sprintf('rgcResponse_Cones2RGC%d_absorptionrate',cone2RGCRatio);

    % Load data
%     template = load(fullfile(baseFolder, 'data', expName, 'rgc', fnameTemplate));
%     template = template.rgcResponse;
    
    template = data{c==expParams.contrastLevels};
    
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
    
    %% 3. Analytical solution
    
    % Get all trials with the same label and only take one trial and one time point
    % (since data are all the same across time, without any photon noise or eyemovement or
    % phase shifts)
    alphaMean = template(label==1,:,:,1:28);
    templateCW = sum(alphaMean(1,:,:,:),4); % sum across all time points to have a fair comparison to the SVM results.
    templateCW = templateCW(:);
    betaMean = template(label==-1,:,:,1:28);
    templateCCW = sum(betaMean(1,:,:,:),4); % sum across all time points to have a fair comparison to the SVM results.
    templateCCW = templateCCW(:);
    
    if templateCW==templateCCW
        dprimeAnalytical(c==expParams.contrastLevels) = 0;
        percentCorrect(c==expParams.contrastLevels) = 0.5;
        
        fprintf('Contrast %1.4f \t d-prime: %2.3f,  percent correct: %2.3f\n', c, dprimeAnalytical(c==expParams.contrastLevels),percentCorrect(c==expParams.contrastLevels))

    else
        numerator = sum( (templateCCW-templateCW).*log(templateCCW./templateCW) );
        denominator = sqrt( 0.5* sum( (templateCW+templateCCW) .* (log(templateCCW./templateCW)).^2 ));
    
        dprime = numerator/denominator;

        dprimeAnalytical(c==expParams.contrastLevels) = dprime;
        percentCorrect(c==expParams.contrastLevels) = normcdf(dprime/2);
        
        fprintf('Contrast %1.4f \t d-prime: %2.3f,  percent correct: %2.3f\n', c, dprimeAnalytical(c==expParams.contrastLevels),percentCorrect(c==expParams.contrastLevels))

    end
end

saveFolderClassification = fullfile(baseFolder, 'data', 'classification', 'rgc', expName, subFolder);
if ~exist('saveFolderClassification', 'dir'); mkdir(saveFolderClassification); end;

fnameClassify = sprintf('ideal_Classify_coneOutputs_contrast%1.4f_pa0_eye00_eccen4.50_defocus0.00_noise-none_sf4.00_lms-1.00.00.0_RGC%d', max(expParams.contrastLevels), cone2RGCRatio);
accuracy = percentCorrect.*100;
parsave(fullfile(saveFolderClassification, sprintf('%s.mat', fnameClassify)),'accuracy',accuracy);

if expParams.verbose
    logzero = 4e-5;
    figure; plot(expParams.contrastLevels(2:end), percentCorrect(2:end), 'o-'); hold on;
    plot(logzero, percentCorrect(1), 'o')
    set(gca, 'XScale', 'log', 'YLim', [.4, 1])
end

return
%% Simulation

% If absorption data are the same, it doesn't matter..
if all(templateCW == templateCCW)
    percentCorrect = 0.5;
    dPrime = 0;
    
    
else %% classify the template data
    index = find(templateCW ~= templateCCW);
    
    expNameData = 'defaultnophaseshift';
    expParams = loadExpParams(expNameData, false);
    subFolderName = 'run1';
    
    if saveFigs
        figurePth = fullfile(ogRootPath, 'figs', expNameData, subFolderName);
        if ~exist('figurePth','dir'); mkdir(figurePth); end;
    end
    
    for c = expParams.contrastLevels
        
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
            ZGivenCW(ii)= sum(CWData(ii,:)' .* log(templateCCW./templateCW));
            ZGivenCCW(ii)= sum(CCWData(ii,:)' .* log(templateCCW./templateCW));
        end
        
        % Compute mean and varianec across trials
        meanZGivenCCW = mean(ZGivenCCW);
        meanZGivenCW = mean(ZGivenCW);
        varZGivenCCW = var(ZGivenCCW);
        varZGivenCW = var(ZGivenCW);
        
        % Compute d-prime
        dprime = (meanZGivenCCW-meanZGivenCW) / sqrt(0.5*(varZGivenCCW+varZGivenCW));
        
        
        % Plot loglikelihood distributions of the two stimulus classes
        figure(1); clf; set(gcf, 'Color', 'w', 'Position', [508, 820, 1052, 518]);
        subplot(1,2,1); hold all;
        [nA,xA] = hist(ZGivenCW,40);
        bar(xA,nA, 'EdgeColor', 'k', 'FaceColor', 'w')
        predA = length(ZGivenCW)*normpdf(xA,meanZGivenCW,sqrt(varZGivenCW))*(xA(2)-xA(1));
        plot(xA,predA,'g','LineWidth',4);
        set(gca, 'TickDir', 'out', 'FontSize', 14);
        
        [nB,xB] = hist(ZGivenCCW,40);
        bar(xB,nB, 'EdgeColor', 'k', 'FaceColor', 'k');
        predB = length(ZGivenCCW)*normpdf(xB,meanZGivenCCW,sqrt(varZGivenCCW))*(xB(2)-xB(1));
        plot(xB,predB,'r','LineWidth',4);
        title(sprintf('Probability of absorptions for Alpha vs Beta stim, contrast: %1.4f',c))
        legend({'Alpha samples', 'Alpha fit', 'Beta samples', 'Beta fit'}, 'FontSize', 14, 'Location', 'best'); legend boxoff;
        
        % Get ROC curve from Hit and False Alarm rates and sample along many criteria
        criteria = linspace(min([ZGivenCW,ZGivenCCW]),max([ZGivenCW,ZGivenCCW]),1000);
        for cc = 1:length(criteria)
            HitRate(cc) = length(find(ZGivenCCW > criteria(cc))) / length(ZGivenCCW);
            FARate(cc) = length(find(ZGivenCW > criteria(cc))) / length(ZGivenCCW);
        end
        
        % Get the percent correct
        percentCorrect(c==expParams.contrastLevels) = -trapz([1 FARate 0],[1 HitRate 0]);
        subplot(1,2,2);
        plot([1 FARate 0],[1 HitRate 0],'r','LineWidth',4);
        xlabel('False Alarm Rate');
        ylabel('Hit Rate');
        title(sprintf('ROC Curve, d-prime: %2.2f, percent corrent: %2.2f',dprime,percentCorrect(c==expParams.contrastLevels)*100));
        axis('square'); axis([0 1 0 1]); set(gca, 'TickDir', 'out', 'FontSize', 14);
        drawnow;
        
        
        
        % Report percent correct
        fprintf('For contrast %1.4f: d-prime is %2.2f  percent correct %2.2f%%\n',c, dprime,percentCorrect(c==expParams.contrastLevels)*100);
        
        if saveFigs
            savefig(fullfile(figurePth,sprintf('Geisler_Ideal_dprime_%s_c%1.4f.fig',expName,c)))
            hgexport(gcf,fullfile(figurePth,sprintf('Geisler_Ideal_dprime_%s_c%1.4f.eps',expName,c)))
        end
    end
    
    saveFolderClassification = fullfile(ogRootPath, 'data', 'classification', expNameData, subFolderName);
    fnameClassify = sprintf('ideal_Classify_coneOutputs_contrast%1.3f_pa0_eye00_eccen4.50_defocus0.00_noise-random_sf4.00_lms-0.60.30.1', max(expParams.contrastLevels));
    accuracy = percentCorrect.*100;
    parsave(fullfile(saveFolderClassification, sprintf('%s.mat', fnameClassify)),'accuracy',accuracy);
    
end

return