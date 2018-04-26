%% s_makePsychometricFunctions


% Script to compute psychometric functions based on the computational
% observer model

%% 0. Set general experiment parameters
expName                  = 'defocus';
expParams                = loadExpParams(expName, false);
[xUnits, colors, labels, allDensity] = loadWeibullPlottingParams(expName);

% Use cone current (flag = true) or cone absorptions (flag = false)
currentFlag = false;
polarAngles = expParams.polarAngle;
FFTflag     = true;
saveFig     = true;

% Where to find data and save figures
dataPth     = fullfile(ogRootPath,'data','classification','HPC',expName,'100trials');
figurePth   = fullfile(ogRootPath,'figs','HPC');

% Number of total trials in computational observer model (50 clockwise, 50 counterclockwise)
nTotal      = 100;%expParams.nTrials*4;

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
            
            fName   = sprintf('Classify_coneOutputs_contrast%1.3f_pa%d_eye%s_eccen%1.2f_defocus%1.2f_noise-random_sf%1.2f.mat', ...
                max(expParams.contrastLevels),polarAngles,sprintf('%i',expParams.eyemovement(:,em)'),expParams.eccentricities(eccen),expParams.defocusLevels(df),expParams.spatFreq);
            if currentFlag; fName = ['current_' fName]; end;
            
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
%             fit.ctrthresh{count} = xUnits(minval(1));
            fit.ctrthresh{count} = fit.ctrvar{count}(2);
            fit.data{count} = accuracy.P;
            
            count = count +1;
        end
        
    end
end


%% 4. Visualize psychometric curves

figure(3); clf; set(gcf,'Color','w', 'Position',  [1000, 850, 986, 488]); hold all;

for ii = 1:length(fit.ctrpred)
    dataToFit = fit.data{ii};
    plot(xUnits(2:end), fit.ctrpred{ii}(2:end)*100, 'Color', colors(ii,:), 'LineWidth',2);
    scatter(expParams.contrastLevels(2:end), dataToFit(2:end), 80, colors(ii,:), 'filled');
    plot(3e-3,dataToFit(1),'o','Color',colors(ii,:), 'MarkerSize', 8, 'MarkerFaceColor',colors(ii,:))
end

set(gca, 'XScale','log','XLim',[3e-3, max(expParams.contrastLevels)],'YLim', [min(dataToFit)-10 100], 'TickDir','out','TickLength',[.015 .015],'FontSize',17, 'LineWidth',2);
set(gca, 'XTick', [3e-3, expParams.contrastLevels(2:2:end)], 'XTickLabel',sprintfc('%1.1f',[0 expParams.contrastLevels(2:2:end)]*100))

ylabel('Classifier Accuracy (% Correct)', 'FontSize',17)
xlabel('Stimulus Contrast (%)', 'FontSize',17);

h = findobj(gca,'Type','line');
legend([h(end:-2:2)],labels, 'Location','bestoutside'); legend boxoff

if saveFig
    if ~exist(figurePth,'dir'); mkdir(figurePth); end
    savefig(fullfile(figurePth,sprintf('WeibullFit_contrastVSperformance_fft%d_%s_100trials',FFTflag,expName)))
    hgexport(gcf,fullfile(figurePth,sprintf('WeibullFit_contrastVSperformance_fft%d_%s_100trials.eps',FFTflag,expName)))
end

%% Plot density thresholds
if all(ismember('coneDensity',expName)) || strcmp('eccbasedcoverage',expName)
    
    thresh = cell2mat(fit.ctrthresh);
%     M  = allDensity/11.111; % Convert mm2 to deg2 [Note: not needed
%     anymore, weibull fit will be on cones/deg2 data]
    lm = fitlm(log10(M),thresh);
    
    figure(2); clf; set(gcf, 'Color', 'w', 'Position', [1318, 696, 836, 649])
    plot(lm, 'LineWidth', 3, 'MarkerSize',10, 'Marker','o','Color',[0 0 0]); box off;
    set(gca, 'TickDir', 'out','TickLength',[0.015 0.015], 'LineWidth',1,'Fontsize',25,'XScale','linear')
    xlabel('Cone Density (cones/deg^2)','FontSize',25); ylabel('Contrast sensitivity threshold','FontSize',25)
    set(gca, 'XTick',[2, 3, 4],'XTickLabel',[100 1000 10000], 'XLim', [1.99 5],'YLim', [0 0.04]),
    legend off; title('Contrast threshold versus Cone density')
    
    if saveFig
        if ~exist(figurePth,'dir'); mkdir(figurePth); end
        savefig(fullfile(figurePth,sprintf('contrastThresholdVS%s_fft%d_100trials',expName,FFTflag)))
        hgexport(gcf,fullfile(figurePth,sprintf('contrastThresholdVS%s_fft%d_100trials',expName,FFTflag)))
    end
    
elseif strcmp(expName,'defocus')
    
    thresh = cell2mat(fit.ctrthresh);
    lm = fitlm(M,thresh);
    
    figure(2); clf; set(gcf, 'Color', 'w', 'Position', [1318, 696, 836, 649])
    plot(lm, 'LineWidth', 3, 'MarkerSize',10, 'Marker','o','Color',[0 0 0]); box off;
    set(gca, 'TickDir', 'out','TickLength',[0.015 0.015], 'LineWidth',1,'Fontsize',25,'XScale','linear')
    xlabel('Defocus (diopters)','FontSize',25); ylabel('Contrast sensitivity threshold','FontSize',25)
    %         set(gca, 'XTick',[2, 3, 4],'XTickLabel',[100 1000 10000], 'XLim', [1.99 5],'YLim', [-0.01 0.08]),
    legend off; title('Contrast threshold versus Defocus')
    
    if saveFig
        if ~exist(figurePth,'dir'); mkdir(figurePth); end
        savefig(fullfile(figurePth,sprintf('contrastThresholdVS%s_fft%d_100trials',expName,FFTflag)))
        hgexport(gcf,fullfile(figurePth,sprintf('contrastThresholdVS%s_fft%d_100trials',expName,FFTflag)))
    end
    
end
