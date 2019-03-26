function [] = plotPsychometricFunctionsIdeal(expName, varargin)
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
% subFolderNames{1} = 'idealtemplate';
% subFolderNames{2} = 'idealsimulation';
% subFolderNames{3} = 'svm';
% subFolderNames{4} = 'svm400trials';
%
% expName = 'idealobserver';
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
fit.thresh = 0.75;

%% 2. Set up figure
figure(3); clf; set(gcf,'Color','w', 'Position',  [1000, 850, 986, 488], 'NumberTitle', 'off', 'Name', sprintf('Psychometric function condition: %s', expName)); hold all;

% Define a zero point (just a very small number), to plot the 0 contrast,
% since a log-linear plot does not define 0.
logzero = 1e-5;
markerColors = colors;
markerColors(2,:) = [1 1 1];

for ii = 1:length(subFolderNames)
    
    %% 3. Load data
    % Define data path
    dataPth     = fullfile(ogRootPath,'data','PF_data_alias','classification',expName, subFolderNames{ii});
    
    % Reload exp params
    switch subFolderNames{ii}
        case 'idealtemplate'
            expParams = loadExpParams('idealobserver',false);
            fit.init   = [2, 0.0005]; % slope, threshold at ~80%
        case 'idealsimulation'
            expParams = loadExpParams('defaultnophaseshift',false);
            fit.init   = [2, 0.0005]; % slope, threshold at ~80%
        case 'svm'
            expParams = loadExpParams('defaultnophaseshift',false);
            fit.init   = [2, 0.01]; % slope, threshold at ~80% 
        case 'svm400trials'
           expParams = loadExpParams('defaultnophaseshift',false);
           fit.init   = [2, 0.01]; % slope, threshold at ~80% 
           expParams.contrastLevels = [0:0.001:0.01, 0.015, 0.02, 0.03, 0.04];
    end
    
   

    % Get file name
    fName{ii}   = sprintf('Classify_coneOutputs_contrast%1.4f_pa0_eye%s_eccen%1.2f_defocus%1.2f_noise-%s_sf%1.2f_lms-%1.1f%1.1f%1.1f', ...
        max(expParams.contrastLevels), ...
        sprintf('%i',expParams.eyemovement'), ...
        expParams.eccentricities, ...
        expParams.defocusLevels, ...
        expParams.cparams.noise, ...
        expParams.spatFreq, ...
        expParams.cparams.spatialDensity(1,2), ...
        expParams.cparams.spatialDensity(1,3), ...
        expParams.cparams.spatialDensity(1,4));

    % Check if ideal observer model or computational model was used to
    % classify
    if regexp(subFolderNames{ii}, 'ideal', 'ONCE')
        % change file name
        fName{ii} = ['ideal_' fName{ii}];
    end
    
    % Check if average or not was used    
    d = dir(fullfile(dataPth, [fName{ii} '_AVERAGE.mat']));
    if ~isempty(d)
        if exist(fullfile(d.folder, d.name), 'file')
            fNameSE{ii} = [fName{ii} '_SE.mat'];
            fName{ii}   = d.name;
        end
        
    else
        fNameSE{ii} = [];
        fName{ii} = [fName{ii} '.mat'];
    end
        
    % Load model performance
    accuracy = load(fullfile(dataPth, fName{ii}));
    fn = fieldnames(accuracy);
    accuracy.P = squeeze(accuracy.(fn{1}));
    
    % Transpose matrix if necessary
    if size(accuracy.P,1)<size(accuracy.P,2)
        accuracy.P = accuracy.P';
    end
    
    %% 4. Fit Weibull
    % Make a Weibull function first with contrast levels and then search for
    % the best fit with the classifier data
    fit.ctrvar{ii} = fminsearch(@(x) ogFitWeibull(x, expParams.contrastLevels, accuracy.P, nTotal), fit.init);
    
    % Then fit a Weibull function again, but now with the best fit parameters
    % from the previous step.
    fit.ctrpred{ii} = ogWeibull(fit.ctrvar{ii}, xUnits);
    
    % Find contrast threshold
    fit.ctrthresh{ii} = fit.ctrvar{ii}(2);
    fit.data{ii} = accuracy.P;
    
    %% 5. Visualize psychometric curves

    % What to plot?
    dataToPlot = fit.data{ii};
    fitToPlot  = fit.ctrpred{ii}*100;

    % Plot curve
    plot(xUnits(2:end), fitToPlot(2:end), 'Color', colors(ii,:), 'LineWidth',2, 'LineStyle', lineStyles{ii});
    scatter(expParams.contrastLevels(2:end), dataToPlot(2:end), 80, colors(ii,:),'MarkerEdgeColor',colors(ii,:), 'MarkerFaceColor', markerColors(ii,:));
    
    % Plot zero point on arbitrary small 'logzero' point, since log axis
    % does not have an actual origin point
    plot(logzero,dataToPlot(1),'o','Color',colors(ii,:), 'MarkerSize', 8, 'MarkerEdgeColor',colors(ii,:), 'MarkerFaceColor', markerColors(ii,:))
    
    if ~isempty(fNameSE{ii})
        SE{ii} = load(fullfile(dataPth, fNameSE{ii}));
        errorbar([logzero, expParams.contrastLevels(2:end)], dataToPlot, SE{ii}.P_SE,'Color', colors(ii,:), 'LineStyle','none');
    end
end

xmax = 0.04;
xticks = [logzero, expParams.contrastLevels(2:9:end)]; 

set(gca, 'XScale','log','XLim',[logzero, xmax],'YLim', [40 100], 'TickDir','out','TickLength',[.015 .015],'FontSize',17, 'LineWidth',2);
set(gca, 'XTick', xticks, 'XTickLabel',sprintfc('%1.2f',xticks.*100))

ylabel('Classifier Accuracy (% Correct)', 'FontSize',17)
xlabel('Stimulus Contrast (%)', 'FontSize',17);

h = findobj(gca,'Type','line');
legend([h(end:-2:2)],labels, 'Location','bestoutside'); legend boxoff

if saveFig
    if ~exist(figurePth,'dir'); mkdir(figurePth); end
    savefig(fullfile(figurePth,sprintf('WeibullFit_contrastVSperformance_%s_idealvsSVM',expName)))
    hgexport(gcf,fullfile(figurePth,sprintf('WeibullFit_contrastVSperformance_%s_idealvsSVM.eps',expName)))
end

