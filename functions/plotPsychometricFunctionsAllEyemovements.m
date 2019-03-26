function [] = plotPsychometricFunctionsAllEyemovements(expName, varargin)
% Function to compute psychometric functions based on the computational
% observer model.
%
% INPUTS:
% expName         : string defining the condition you want to plot.
%                   (See load expParams for possible conditions)
% [subFolderNames] : string defining the sub folder you want to plot from.
%
% [saveFig]       : boolean defining to save figures or not
% [plotAvg]      : boolean defining to plot average across experiments runs or not
%
%
% Example:
% subFolderNames{1} = 'eyemovnophaseshift';
% subFolderNames{2} = 'eyemov';
% expName = 'eyemovnophaseshift';
% plotPsychometricFunctionsAllEyemovements(expName, 'subFolderNames', subFolderNames, 'saveFig', true)

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
[xUnits, colors, labels, ~, ~] = loadWeibullPlottingParams(subFolderNames{1});

% Only plot with all eyemovements or none (skip drift only)
conditionsToPlot = [1, size(expParams.eyemovement,2)];
[~, ~, labels2, ~, ~] = loadWeibullPlottingParams(subFolderNames{2});
labels = {labels{conditionsToPlot}, labels2{conditionsToPlot}};

startIdx=0.0005;

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

% Define a zero point (just a very small number), to plot the 0 contrast,
% since a log-linear plot does not define 0.
logzero = 4e-4;
xmax    = 0.1;
xt = find(expParams.contrastLevels==startIdx);
xticks  = [logzero, expParams.contrastLevels(xt), 0.001, 0.01, 0.1]; 
markerFaceColors{1} = colors(conditionsToPlot,:);
markerFaceColors{1}(2,:) = [1 1 1];
markerFaceColors{2} = markerFaceColors{1};
markerFaceColors{2}(1,:) = [1 0 0];
lineStyles = {'-', ':', ':'};

% Set up figure
figure(3); clf; set(gcf,'Color','w', 'Position',  [1000, 850, 986, 488], 'NumberTitle', 'off', 'Name', sprintf('Psychometric function condition: %s', expName)); hold all;

%% 2. Load classification accuracy data
for ii = 1:length(subFolderNames)
    
    % Define data path
    dataPth     = fullfile(ogRootPath,'data','PF_data_alias','classification',subFolderNames{ii}, 'average');
    
    % Reload exp params
    switch subFolderNames{ii}
        case 'eyemov'
            expParams = loadExpParams(subFolderNames{ii},false);
            fit.init{1}   = [2, 0.02]; % slope, threshold at ~80% 
            fit.init{2}   = [2, 0.02]; % slope, threshold at ~80% 
        case 'eyemovnophaseshift'
            expParams = loadExpParams(subFolderNames{ii},false);
            fit.init{1}   = [2, 0.0005]; % slope, threshold at ~80% 
            fit.init{2}   = [1, 0.01]; % slope, threshold at ~80% 
    end
    
    for emIdx = conditionsToPlot % only plot with  all or without any eye movements

        % Get file name
        fName{ii, emIdx}   = sprintf('Classify_coneOutputs_contrast%1.4f_pa0_eye%s_eccen%1.2f_defocus%1.2f_noise-%s_sf%1.2f_lms-%1.1f%1.1f%1.1f', ...
            max(expParams.contrastLevels), ...
            sprintf('%i',expParams.eyemovement(:,emIdx)'), ...
            expParams.eccentricities, ...
            expParams.defocusLevels, ...
            expParams.cparams.noise, ...
            expParams.spatFreq, ...
            expParams.cparams.spatialDensity(1,2), ...
            expParams.cparams.spatialDensity(1,3), ...
            expParams.cparams.spatialDensity(1,4));

        fNameAVERAGE{ii, emIdx} = [fName{ii, emIdx} '_AVERAGE.mat'];
        fNameSE{ii, emIdx} = [fName{ii, emIdx} '_SE.mat'];

        % Load model performance
        accuracy = load(fullfile(dataPth, fNameAVERAGE{ii, emIdx}));
        fn = fieldnames(accuracy);
        accuracy.P = squeeze(accuracy.(fn{1}));

        % Transpose matrix if necessary
        if size(accuracy.P,1)<size(accuracy.P,2)
            accuracy.P = accuracy.P';
        end
    
        %% 4. Fit Weibull
        % Make a Weibull function first with contrast levels and then search for
        % the best fit with the classifier data
        fit.ctrvar{ii, emIdx} = fminsearch(@(x) ogFitWeibull(x, expParams.contrastLevels, accuracy.P, nTotal), fit.init{emIdx==conditionsToPlot});

        % Then fit a Weibull function again, but now with the best fit parameters
        % from the previous step.
        fit.ctrpred{ii, emIdx} = ogWeibull(fit.ctrvar{ii, emIdx}, xUnits);

        % Find contrast threshold
        fit.ctrthresh{ii, emIdx} = fit.ctrvar{ii, emIdx}(2);
        fit.data{ii, emIdx} = accuracy.P;

        %% 5. Visualize psychometric curves

        % What to plot?
        switch subFolderNames{ii}
            case 'eyemovnophaseshift'
                idxStartData = find(expParams.contrastLevels==startIdx);
                [~, idxStartFit] = min(abs(xUnits-startIdx));
            case 'eyemov'
                idxStartData = 2;
        end
        
        dataToPlot = fit.data{ii, emIdx};
        fitToPlot  = fit.ctrpred{ii, emIdx}*100;

        % Plot curve
        plot(xUnits(idxStartFit:end), fitToPlot(idxStartFit:end), 'Color', colors(ii,:), 'LineWidth',2, 'LineStyle', lineStyles{emIdx==conditionsToPlot});
        scatter(expParams.contrastLevels(idxStartData:end), dataToPlot(idxStartData:end), 80, colors(ii,:),'MarkerEdgeColor',colors(ii,:), 'MarkerFaceColor', markerFaceColors{ii}(emIdx==conditionsToPlot,:));

        % Plot zero point on arbitrary small 'logzero' point, since log axis
        % does not have an actual origin point
        plot(logzero,dataToPlot(1),'o','Color',colors(ii,:), 'MarkerSize', 8, 'MarkerEdgeColor',colors(ii,:), 'MarkerFaceColor', markerFaceColors{ii}(emIdx==conditionsToPlot,:))

        if ~isempty(fNameSE{ii})
            SE{ii} = load(fullfile(dataPth, fNameSE{ii}));
            errorbar([logzero, expParams.contrastLevels(idxStartData:end)], dataToPlot([1, (idxStartData:end)]), SE{ii}.P_SE([1, (idxStartData:end)]) ,'Color', colors(ii,:), 'LineStyle','none');
        end
    end
end


set(gca, 'XScale','log','XLim',[logzero, xmax],'YLim', [40 100], 'TickDir','out','TickLength',[.015 .015],'FontSize',17, 'LineWidth',2);
set(gca, 'XTick', xticks, 'XTickLabel',sprintfc('%1.2f',xticks*100))

ylabel('Classifier Accuracy (% Correct)', 'FontSize',17)
xlabel('Stimulus Contrast (%)', 'FontSize',17);

h = findobj(gca,'Type','line');
legend([h(end:-2:2)],labels, 'Location','bestoutside'); legend boxoff

if saveFig
    if ~exist(figurePth,'dir'); mkdir(figurePth); end
    savefig(fullfile(figurePth,sprintf('WeibullFit_contrastVSperformance_%s_eyemovVSphaseshift',expName)))
    hgexport(gcf,fullfile(figurePth,sprintf('WeibullFit_contrastVSperformance_%s_eyemovVSphaseshift.eps',expName)))
end

