function plotConeTypesVSThreshold(expName, fit, xThresh, varargin)
% Function to plot L:M mixture levels versus stimulus contrast thresholds.
%
% INPUTS:
% expName         : string defining the condition you want to plot.
%                   (See load expParams for possible conditions)
% fit             : struct with fit data
%
% xThresh         : vector with x units for plot
% [saveFig]       : boolean defining to save figure or not
% [figurePth]     : boolean defining directory where to save figure

%% 0. Parse input parameters
p = inputParser;
p.KeepUnmatched = true;
p.addRequired('expName', @ischar);
p.addRequired('fit', @isstruct);
p.addRequired('xThresh', @isvector);
p.addParameter('saveFig', false, @islogical);
p.addParameter('figurePth', fullfile(ogRootpath, 'figs'), @isdir);
p.parse(expName, varargin{:});

% Rename variables
expName       = p.Results.expName;
fit           = p.Results.fit;
xThresh       = p.Results.xThresh;
saveFig       = p.Results.saveFig;
figurePth     = p.Results.figurePth;

%% 1. Fit a quadratic function to thresholds as a function of defocus levels

yThresh = cell2mat(fit.ctrthresh);
lm = fitlm(xThresh,yThresh, 'quadratic');
    
%% 2. Visualize plot
figure(2); clf; set(gcf, 'Color', 'w', 'Position', [1318, 696, 836, 649])
plot(lm, 'LineWidth', 3, 'MarkerSize',10, 'Marker','o','Color',[0 0 0]); box off;
set(gca, 'TickDir', 'out','TickLength',[0.015 0.015], 'LineWidth',1,'Fontsize',25,'XScale','linear', 'XLim', [-10,110], 'YLim', [0,.016])
xlabel('Probability of L-cones in L:M cone ratio ','FontSize',25); ylabel('Contrast threshold','FontSize',25)
legend off; title(sprintf('Contrast threshold vs probability of L-cones - r2: %1.2f', lm.Rsquared.ordinary))

% Save figure if requested
if saveFig
    if ~exist(figurePth,'dir'); mkdir(figurePth); end
    savefig(fullfile(figurePth,sprintf('contrastThresholdVS%s',expName)))
    hgexport(gcf,fullfile(figurePth,sprintf('contrastThresholdVS%s.eps',expName)))
end