%% s_makePsychometricFunctions


% Script to compute psychometric functions based on the computational
% observer model

%% 0. Set general experiment parameters
expName                  = 'defocus';
expParams                = loadExpParams(expName, false);
[xUnits, colors, labels, M] = loadWeibullPlottingParams(expName);

% Use cone current (flag = true) or cone absorptions (flag = false)
currentFlag = false;
polarAngles = expParams.polarAngle;
FFTflag     = true;
saveFig     = true;

% Where to find data and save figures
subFolderName = 'average';
dataPth     = fullfile(ogRootPath,'data','PF_data_alias','classification',expName,subFolderName);
figurePth   = fullfile(ogRootPath,'figs', expName, [subFolderName]);

% Number of total trials in computational observer model (50 clockwise, 50 counterclockwise)
nTotal      = expParams.nTrials;

%% 1. Set Weibull fitting parameters

% Prepare fit variables
fit = [];

% Predefine cell arrays
fit.ctrpred = cell(size(colors,1),1);
fit.ctrvar  = cell(size(colors,1),1);
fit.ctrr2   = cell(size(colors,1),1);
fit.data    = cell(size(colors,1),1);

% Set inital slope, threshold for first stage fitting
fit.init   = [2, 0.02]; % slope, threshold at ~80%
fit.thresh = 0.75;

% Get nr of conditions
nrEyemovTypes    = size(expParams.eyemovement,2);
nrEccen          = length(expParams.eccentricities);
nrSpatFreq       = length(expParams.spatFreq);
nrDefocusLevels  = length(expParams.defocusLevels);

% Check if different fitting accuracy for different cone types is requested
% (not implemented yet)
fn = fieldnames(expParams);
if any(strcmp(fn(:),'cparams')); nrConeTypes = size(expParams.cparams.spatialDensity,2);
    error('s_makePsychometricFunctions does not allow multiple cone type conditions (yet)'); end


count = 1;
for em = 1:nrEyemovTypes
    for eccen = 1:nrEccen
        for df = 1:nrDefocusLevels
            
            
            %% 2. Load results
            
            fName   = sprintf('Classify_coneOutputs_contrast%1.3f_pa%d_eye%s_eccen%1.2f_defocus%1.2f_noise-random_sf%1.2f_AVERAGE0.mat', ...
                max(expParams.contrastLevels),polarAngles,sprintf('%i',expParams.eyemovement(:,em)'),expParams.eccentricities(eccen),expParams.defocusLevels(df),expParams.spatFreq);
            if currentFlag; fName = ['current_' fName]; end;
            
             if subFolderName == 'average'
                 fNameSE   = sprintf('Classify_coneOutputs_contrast%1.3f_pa%d_eye%s_eccen%1.2f_defocus%1.2f_noise-random_sf%1.2f_SE0.mat', ...
                     max(expParams.contrastLevels),polarAngles,sprintf('%i',expParams.eyemovement(:,em)'),expParams.eccentricities(eccen),expParams.defocusLevels(df),expParams.spatFreq);
                 SE{count} = load(fullfile(dataPth, fNameSE));
             end
            
            accuracy = load(fullfile(dataPth, fName));
            accuracy.P = squeeze(accuracy.P);
            if size(accuracy.P,1)<size(accuracy.P,2)
                accuracy.P = accuracy.P';
            end
            
            %% 3. Fit Weibull
            % Make a Weibull function first with contrast levels and then search for
            % the best fit with the classifier data
            fit.ctrvar{count} = fminsearch(@(x) ogFitWeibull(x, expParams.contrastLevels, accuracy.P, nTotal), fit.init);
            
            % Then fit a Weibull function again, but now with the best fit parameters
            % from the previous step.
            fit.ctrpred{count} = ogWeibull(fit.ctrvar{count}, xUnits);
            
            %% 3. Find contrast threshold
%             diff   = abs(fit.ctrpred{count} - fit.thresh);
%             minval = find(diff == min(diff));
            %fit.ctrthresh{count} = xUnits(minval(1));
            fit.ctrthresh{count} = fit.ctrvar{count}(2);
            fit.data{count} = accuracy.P;
            
            count = count +1;
        end
        
    end
end


%% 4. Visualize psychometric curves

figure(3); clf; set(gcf,'Color','w', 'Position',  [1000, 850, 986, 488]); hold all;

idx = ~cellfun(@isempty, fit.ctrpred);
fit.ctrpred = fit.ctrpred(idx);

for ii = 1:length(fit.ctrpred)
    dataToFit = fit.data{ii};
    plot(xUnits(2:end), fit.ctrpred{ii}(2:end)*100, 'Color', colors(ii,:), 'LineWidth',2);
    scatter(expParams.contrastLevels(2:end), dataToFit(2:end), 80, colors(ii,:), 'filled');

    plot(3e-3,dataToFit(1),'o','Color',colors(ii,:), 'MarkerSize', 8, 'MarkerFaceColor',colors(ii,:))
    
    if strcmp(subFolderName, 'average')
        errorbar([3e-3, expParams.contrastLevels(2:end)], dataToFit, SE{ii}.P_SE,'Color', colors(ii,:), 'LineStyle','none');
    end
end

set(gca, 'XScale','log','XLim',[3e-3, max(expParams.contrastLevels)],'YLim', [40 100], 'TickDir','out','TickLength',[.015 .015],'FontSize',17, 'LineWidth',2);
set(gca, 'XTick', [3e-3, expParams.contrastLevels(2:2:end)], 'XTickLabel',sprintfc('%1.1f',[0 expParams.contrastLevels(2:2:end)]*100))

ylabel('Classifier Accuracy (% Correct)', 'FontSize',17)
xlabel('Stimulus Contrast (%)', 'FontSize',17);

h = findobj(gca,'Type','line');
legend([h(end:-2:2)],labels, 'Location','bestoutside'); legend boxoff

if saveFig
    if ~exist(figurePth,'dir'); mkdir(figurePth); end
    savefig(fullfile(figurePth,sprintf('WeibullFit_contrastVSperformance_fftFlag%d_%s_currentFlag%d',FFTflag,expName, currentFlag)))
    hgexport(gcf,fullfile(figurePth,sprintf('WeibullFit_contrastVSperformance_fftFlag%d_%s_currentFlag%d.eps',FFTflag,expName, currentFlag)))
end

%% Plot density thresholds
if strcmp('coneDensity',expName) || strcmp('eccbasedcoverage',expName)
    
    thresh = cell2mat(fit.ctrthresh);
    lm = fitlm(log10(M),thresh);
    
    figure(2); clf; set(gcf, 'Color', 'w', 'Position', [1318, 696, 836, 649])
    plot(lm, 'LineWidth', 3, 'MarkerSize',10, 'Marker','o','Color',[0 0 0]); box off;
    set(gca, 'TickDir', 'out','TickLength',[0.015 0.015], 'LineWidth',1,'Fontsize',25,'XScale','linear')
    xlabel('Cone Density (cones/deg^2)','FontSize',25); ylabel('Contrast sensitivity threshold','FontSize',25)
    set(gca, 'XTick',[2, 3, 4],'XTickLabel',[100 1000 10000], 'XLim', [1.99 5],'YLim', [0 0.04]),
    legend off; title('Contrast threshold versus Cone density')
    
    if saveFig
        if ~exist(figurePth,'dir'); mkdir(figurePth); end
        savefig(fullfile(figurePth,sprintf('contrastThresholdVS%s_fftFlag%d_currentFlag%d',expName,FFTflag,currentFlag)))
        hgexport(gcf,fullfile(figurePth,sprintf('contrastThresholdVS%s_fftFlag%d_currentFlag%d',expName,FFTflag,currentFlag)))
    end
    
    b_intcpt = lm.Coefficients.Estimate(1);
    a_coeff  = lm.Coefficients.Estimate(2);
    
    cThreshold = @(x) (a_coeff*x) + b_intcpt;
    diopters = @(y) (y-b_intcpt)./a_coeff;
    reportedBehavior = [b_intcpt; b_intcpt+0.015]; % contrast thresholds reported in behavior, corresponding to horizontal versus upper vertical meridian
    modelPredictionForPF = diopters(reportedBehavior);
    
    totalVariance.dioptersPredictedByModel = diff(modelPredictionForPF);
    totalVariance.dioptersReportedInLiterature = 4.1111e+04; % From Song's paper
    totalVariance.contributionOfDefocusPercent = (totalVariance.dioptersReportedInLiterature / totalVariance.dioptersPredictedByModel) * 100;
    
    figure(10); clf; set(gcf, 'Color', 'w', 'Position', [300, 982, 289, 363])
    bar([0.2 0.3], [totalVariance.dioptersPredictedByModel, totalVariance.dioptersReportedInLiterature], 'barWidth', 0.3, 'FaceColor', [0.5 0.5 0.5], 'EdgeColor', 'k');
    set(gca, 'TickDir', 'out', 'XTick', [0.2,0.3], 'XTickLabel', {'Model prediction for PF', 'Max reported difference literature'}, 'XTickLabelRotation', 45, 'FontSize', 9);
    ylim([0 max(totalVariance.dioptersPredictedByModel)+ 0.1*totalVariance.dioptersPredictedByModel]); 
    xlim([0.15 0.37]); box off;
    ylabel('Defocus (Diopters)')
    
    fprintf('Total contribution of cone density according to computational observer model: %1.1f percent\n', totalVariance.contributionOfDefocusPercent)
    
     if saveFig
        if ~exist(figurePth,'dir'); mkdir(figurePth); end
        savefig(fullfile(figurePth,sprintf('expVar_modelVSLiterature%s_fftFlag%d_currentFlag%d',expName,FFTflag, currentFlag)))
        hgexport(gcf,fullfile(figurePth,sprintf('expVar_modelVSLiterature%s_fftFlag%d_currentFlag%d',expName,FFTflag, currentFlag)))
    end
    
elseif strcmp(expName,'defocus')
    
    thresh = cell2mat(fit.ctrthresh);
    lm = fitlm(M(idx),thresh);
    
    figure(2); clf; set(gcf, 'Color', 'w', 'Position', [1318, 696, 836, 649])
    plot(lm, 'LineWidth', 3, 'MarkerSize',10, 'Marker','o','Color',[0 0 0]); box off;
    set(gca, 'TickDir', 'out','TickLength',[0.015 0.015], 'LineWidth',1,'Fontsize',25,'XScale','linear')
    xlabel('Defocus (Diopters)','FontSize',25); ylabel('Contrast threshold','FontSize',25)
    legend off; title(sprintf('Contrast threshold vs level of defocus - R2: %1.2f', lm.Rsquared.ordinary))
    
    if saveFig
        if ~exist(figurePth,'dir'); mkdir(figurePth); end
        savefig(fullfile(figurePth,sprintf('contrastThresholdVS%s_fftFlag%d_currentFlag%d',expName,FFTflag, currentFlag)))
        hgexport(gcf,fullfile(figurePth,sprintf('contrastThresholdVS%s_fftFlag%d_currentFlag%d',expName,FFTflag, currentFlag)))
    end
    
    b_intcpt = lm.Coefficients.Estimate(1);
    a_coeff  = lm.Coefficients.Estimate(2);
    
    cThreshold = @(x) (a_coeff*x) + b_intcpt;
    diopters = @(y) (y-b_intcpt)./a_coeff;
    reportedBehavior = [b_intcpt; b_intcpt+0.015]; % contrast thresholds reported in behavior, corresponding to horizontal versus upper vertical meridian
    modelPredictionForPF = diopters(reportedBehavior);
    
    totalVariance.dioptersPredictedByModel = diff(modelPredictionForPF);
    totalVariance.dioptersReportedInLiterature = 0.2; % From Artal's papers
    totalVariance.contributionOfDefocusPercent = (totalVariance.dioptersReportedInLiterature / totalVariance.dioptersPredictedByModel) * 100;
    
    figure(10); clf; set(gcf, 'Color', 'w', 'Position', [300, 982, 289, 363])
    bar([0.2 0.3], [totalVariance.dioptersPredictedByModel, totalVariance.dioptersReportedInLiterature], 'barWidth', 0.3, 'FaceColor', [0.5 0.5 0.5], 'EdgeColor', 'k');
    set(gca, 'TickDir', 'out', 'XTick', [0.2,0.3], 'XTickLabel', {'Model prediction for PF', 'Max reported difference literature'}, 'XTickLabelRotation', 45, 'FontSize', 9);
    ylim([0 7]); xlim([0.15 0.37]); box off;
    ylabel('Defocus (Diopters)')
    
    fprintf('Total contribution of defocus according to computational observer model: %1.1f percent\n', totalVariance.contributionOfDefocusPercent)
    
     if saveFig
        if ~exist(figurePth,'dir'); mkdir(figurePth); end
        savefig(fullfile(figurePth,sprintf('expVar_modelVSLiterature%s_fftFlag%d_currentFlag%d',expName,FFTflag, currentFlag)))
        hgexport(gcf,fullfile(figurePth,sprintf('expVar_modelVSLiterature%s_fftFlag%d_currentFlag%d',expName,FFTflag, currentFlag)))
    end
    
end