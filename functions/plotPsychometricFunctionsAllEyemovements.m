function [] = plotPsychometricFunctionsAllEyemovements(expName, varargin)


% Function to compute psychometric functions based on the computational
% observer model.

% INPUTS:
% expName         : string defining the condition you want to plot.
%                   (See load expParams for possible conditions)
% [subFolderNames] : string defining the sub folder you want to plot from.

% [saveFig]       : boolean defining to save figures or not
% [plotAvg]      : boolean defining to plot average across experiments runs or not


% Example:
% subFolderNames{1} = 'eyemovnophaseshift';
% subFolderNames{2} = 'eyemov';
% expName = 'eyemov';
% plotPsychometricFunctionsIdeal(expName, 'subFolderNames', subFolderNames, 'saveFig', true)

%% 0. Set general experiment parameters
p = inputParser;
p.KeepUnmatched = true;
p.addRequired('expName', @ischar);
p.addParameter('subFolderNames', [], @iscell)
p.addParameter('saveFig', false, @islogical);
p.parse(expName, varargin{:});

% Rename variables
expName        = p.Results.expName;
subFolderNames = p.Results.subFolderNames;
saveFig        = p.Results.saveFig;

% Load specific experiment parameters
expParams                   = loadExpParams(expName, false);
[xUnits, colors, labels, ~, lineStyles] = loadWeibullPlottingParams(expName);

% Define figurePth
figurePth   = fullfile(ogRootPath,'figs', expName, subFolderNames{1});

nTotal = expParams.nTrials;
%% 1. Set Weibull fitting parameters

% Prepare fit variables
fit = [];

% Predefine cell arrays
fit.ctrpred = cell(size(colors,1),1);
fit.ctrvar  = cell(size(colors,1),1);
fit.ctrr2   = cell(size(colors,1),1);
fit.data    = cell(size(colors,1),1);

% Set inital slope, threshold for first stage fitting
fit.thresh = 0.75;


%% 2. Load classification accuracy

count = 1;

for ii = 1:length(subFolderNames)

    % Load first folder: eye movement experiment without phase shifts
    subFolderName = subFolderNames{ii};
    expParams     = loadExpParams(subFolderName, false);
    xUnits =  linspace(min(expParams.contrastLevels),max(expParams.contrastLevels),200);

    for em = 1:size(expParams.eyemovement,2)
    
        if ii ==1
            fName   = sprintf('Classify_coneOutputs_contrast%1.3f_pa0_eye%s_eccen%1.2f_defocus%1.2f_noise-random_sf%1.2f_lms-%1.1f%1.1f%1.1f.mat', ...
                max(expParams.contrastLevels), ...
                sprintf('%i',expParams.eyemovement(:,em)'), ...
                expParams.eccentricities, ...
                expParams.defocusLevels, ...
                expParams.spatFreq, ...
                expParams.cparams.spatialDensity(1,2), ...
                expParams.cparams.spatialDensity(1,3), ...
                expParams.cparams.spatialDensity(1,4));
            
            dataPth     = fullfile(ogRootPath,'data','PF_data_alias','classification',subFolderName, 'run1');
            fit.init   = [2, 0.001]; % slope, threshold at ~80%

            else

            fName   = sprintf('Classify_coneOutputs_contrast%1.3f_pa0_eye%s_eccen%1.2f_defocus%1.2f_noise-random_sf%1.2f_lms-%1.1f%1.1f%1.1f_AVERAGE.mat', ...
                max(expParams.contrastLevels), ...
                sprintf('%i',expParams.eyemovement(:,em)'), ...
                expParams.eccentricities, ...
                expParams.defocusLevels, ...
                expParams.spatFreq, ...
                expParams.cparams.spatialDensity(1,2), ...
                expParams.cparams.spatialDensity(1,3), ...
                expParams.cparams.spatialDensity(1,4));

            fNameSE   = sprintf('Classify_coneOutputs_contrast%1.3f_pa0_eye%s_eccen%1.2f_defocus%1.2f_noise-random_sf%1.2f_lms-%1.1f%1.1f%1.1f_SE.mat', ...
                max(expParams.contrastLevels), ...
                sprintf('%i',expParams.eyemovement(:,em)'), ...
                expParams.eccentricities, ...
                expParams.defocusLevels, ...
                expParams.spatFreq, ...
                expParams.cparams.spatialDensity(1,2), ...
                expParams.cparams.spatialDensity(1,3), ...
                expParams.cparams.spatialDensity(1,4));

            dataPth     = fullfile(ogRootPath,'data','PF_data_alias','classification',expName, 'toPlot', 'average');
            SE{em}      = load(fullfile(dataPth, fNameSE));
            
            fit.init   = [2, 0.01]; % slope, threshold at ~80%

        end
    
      
    % load model performance
    accuracy = load(fullfile(dataPth, fName));
    fn = fieldnames(accuracy);
    accuracy.P = squeeze(accuracy.(fn{1}));
    
    % Transpose matrix if necessary
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
    
    %% 4. Find contrast threshold
    %             diff   = abs(fit.ctrpred{count} - fit.thresh);
    %             minval = find(diff == min(diff));
    fit.ctrthresh{count} = fit.ctrvar{count}(2);
    fit.data{count} = accuracy.P;
    
    count = count+1;
    end
    
end


%% 5. Visualize psychometric curves

figure(3); clf; set(gcf,'Color','w', 'Position',  [1000, 850, 986, 488], 'NumberTitle', 'off', 'Name', sprintf('Psychometric function condition: %s', expName)); hold all;

% Remove empty cells
idx = ~cellfun(@isempty, fit.ctrpred);
fit.ctrpred = fit.ctrpred(idx);

% Define a zero point (just a very small number), to plot the 0 contrast,
% since a log-linear plot does not define 0.
logzero = 4e-4;
colors = repmat(colors, 2,1);
markerColors = colors;
markerColors(4,:) = [1 1 1];
markerColors(5,:) = [1 1 1];
markerColors(6,:) = [1 1 1];
lineStyles =  {'-', '-', '-', ':',':', ':'};
labels = {'No eye movement,no phase shifts', ...
    'Drift, no phase shifts', ...
    'Drift and MS, no phase shifts', ...
    'No eye movement, with phase shifts', ...
    'Drift, with phase shifts', ...
    'Drift and MS, with phase shifts'};

% Loop over all functions to plot
for ii = 1:length(fit.data)
    
    % What to plot?
    dataToPlot = fit.data{ii};
    fitToPlot  = fit.ctrpred{ii}*100;
    
    if ii <= 3
        expParams.contrastLevels = [0:0.001:0.01, 0.015, 0.02:0.01:0.1];;
    else
        expParams.contrastLevels = [0:0.005:0.04, 0.05:0.01:0.1];      
    end
    
    % Redefine xUnits
    xUnits =  linspace(min(expParams.contrastLevels),max(expParams.contrastLevels),200);

    % Plot curve
    plot(xUnits(2:end), fitToPlot(2:end), 'Color', colors(ii,:), 'LineWidth',3, 'LineStyle', lineStyles{ii});
    scatter(expParams.contrastLevels(2:end), dataToPlot(2:end), 80, colors(ii,:),'MarkerEdgeColor',colors(ii,:), 'MarkerFaceColor', markerColors(ii,:));
    
    % Plot zero point on arbitrary small 'logzero' point, since log axis
    % does not have an actual origin point
    plot(logzero,dataToPlot(1),'o','Color',colors(ii,:), 'MarkerSize', 8, 'MarkerEdgeColor',colors(ii,:), 'MarkerFaceColor', markerColors(ii,:))
    
    if ii > 3
        errorbar([logzero, expParams.contrastLevels(2:end)], dataToPlot, SE{ii-3}.P_SE,'Color', colors(ii,:), 'LineStyle','none');
    end
end

xmax = 0.1;

set(gca, 'XScale','log','XLim',[logzero, xmax],'YLim', [40 100], 'TickDir','out','TickLength',[.015 .015],'FontSize',17, 'LineWidth',2);
set(gca, 'XTick', [logzero, 0.001, 0.01, 0.05, 0.1], 'XTickLabel',sprintfc('%1.2f',[0, 0.001, 0.01, 0.05, 0.1]*100))

ylabel('Classifier Accuracy (% Correct)', 'FontSize',17)
xlabel('Stimulus Contrast (%)', 'FontSize',17);

h = findobj(gca,'Type','line');
legend([h(end:-2:2)],labels, 'Location','bestoutside'); legend boxoff

if saveFig
    if ~exist(figurePth,'dir'); mkdir(figurePth); end
    savefig(fullfile(figurePth,sprintf('WeibullFit_contrastVSperformance_%s_eyemovVSphaseshift',expName)))
    hgexport(gcf,fullfile(figurePth,sprintf('WeibullFit_contrastVSperformance_%s_eyemovVSphaseshift.eps',expName)))
end

