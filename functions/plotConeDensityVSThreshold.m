function [] = plotConeDensityVSThreshold(expName, fit, xThresh, varargin)
% Function to plot cone density levels versus stimulus contrast thresholds.
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


%% 1. Fit a log-linear function to thresholds vs density
yThresh = cell2mat(fit.ctrthresh);
lm = fitlm(log10(xThresh),yThresh);

%% 2. Visualize plot
xrange = unique(round(log10(xThresh)));
xticks = [xrange(1)-1; xrange; xrange(end)+1];
for ii = 1:length(xticks); xticklbls{ii} = sprintf('10^%d', xticks(ii)'); end% get x axis range

% Plot it!
figure(2); clf; set(gcf, 'Color', 'w', 'Position', [1318, 696, 836, 649])
plot(lm, 'LineWidth', 3, 'MarkerSize',10, 'Marker','o','Color',[0 0 0]); box off;
set(gca, 'TickDir', 'out','TickLength',[0.015 0.015], 'LineWidth',1,'Fontsize',25,'XScale','linear')
xlabel('Cone Density (cones/deg^2)','FontSize',25); ylabel('Contrast threshold (%)','FontSize',25)
set(gca, 'XTick',xticks,'XTickLabel',xticklbls, 'XLim', [2 5],'YLim', [0 0.04]),
set(gca, 'YTick',[ 0 .01 .02 .03 .04],'YTickLabel',[0 1 2 3 4]),
legend off; title(sprintf('Contrast threshold vs level of cone density - R2: %1.2f', lm.Rsquared.ordinary))

% Save fig if requested
if saveFig
    if ~exist(figurePth,'dir'); mkdir(figurePth); end
    savefig(fullfile(figurePth,sprintf('contrastThresholdVS%s',expName)))
    hgexport(gcf,fullfile(figurePth,sprintf('contrastThresholdVS%s',expName)))
end

%% 3. Quantify effect of cone density

% Get intercept and slope of log-linear fit
b_intcpt = lm.Coefficients.Estimate(1);
a_coeff  = lm.Coefficients.Estimate(2);

% Reconstruct log-linear function
cThreshold = @(x) (a_coeff* log10(x)) + b_intcpt;
predictedDensity = @(y) 10.^((y-b_intcpt)./a_coeff);

% contrast thresholds reported in behavior, corresponding to horizontal versus upper vertical meridian
reportedBehavior = [0.02+0.015, 0.02]; % (thresholds in % contrast)
modelPredictionForPF = predictedDensity(reportedBehavior);

ang = [0, 90, 180, 270]; % polar angles: nasal (HM), superior (LVM), temporal (HM), inferior (UVM) (radians)
for jj = 1:length(ang)
    densityMM2 = coneDensityReadData('coneDensitySource', 'Song2011Young','eccentricity',4.5,'angle',ang(jj),'eccentricityUnits', 'deg','whichEye','left');
    densityDeg2(jj) = densityMM2/((1/.3)^2);
end

totalVariance.densityPredictedByModel = diff(modelPredictionForPF);
totalVariance.densityReportedInLiterature =  diff([densityDeg2(4), mean(densityDeg2([1 3]))]); % From Song's paper
totalVariance.contributionOfDensityPercent = (totalVariance.densityReportedInLiterature / totalVariance.densityPredictedByModel) * 100;

fprintf('Total contribution of cone density according to computational observer model: %1.1f percent\n', totalVariance.contributionOfDensityPercent)

end
