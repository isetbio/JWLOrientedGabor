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
%
% Example: 
% plotConeDensityVSThreshold('conedensity', 'linear', x 

%% 0. Parse input parameters
p = inputParser;
p.KeepUnmatched = true;
p.addRequired('expName', @ischar);
p.addRequired('fit', @isstruct);
p.addRequired('xThresh', @isvector);
p.addParameter('varThresh',[],  @isvector);
p.addParameter('fitTypeName', 'linear', @ischar);
p.addParameter('yScale', 'log', @ischar);
p.addParameter('saveFig', false, @islogical);
p.addParameter('figurePth', fullfile(ogRootPath, 'figs'), @isdir);
p.parse(expName, fit, xThresh, varargin{:});

% Rename variables
expName       = p.Results.expName;
fit           = p.Results.fit;
xThresh       = p.Results.xThresh;
varThresh     = p.Results.varThresh;
fitTypeName   = p.Results.fitTypeName;
yScale        = p.Results.yScale;
saveFig       = p.Results.saveFig;
figurePth     = p.Results.figurePth;


%% 1. Fit a log-linear function to thresholds vs density

yThresh = cell2mat(fit.ctrthresh);

if strcmp(yScale,'log') % linear in log-log space
    yData = log10(yThresh)';
else
    yData = yThresh';
end

if strcmp(fitTypeName, 'linear')
    [fitResult, err]  = polyfit(log10(xThresh),yData,1);
    [y, delta] = polyval(fitResult,log10(xThresh), err);
    R2 = 1 - (err.normr/norm(yData - mean(yData)))^2;
elseif strcmp(fitTypeName, 'poly2')
    [fitResult, err]  = polyfit(log10(xThresh),yData,2);
    [y, delta] = polyval(fitResult,log10(xThresh), err);
    R2 = 1 - (err.normr/norm(yData - mean(yData)))^2;
end


    
%% 2. Visualize plot
xrange = unique(round(log10(xThresh)));
xticks = [xrange(1)-1; xrange; xrange(end)+1];

for ii = 1:length(xticks); xticklbls{ii} = sprintf('10^%d', xticks(ii)'); end% get x axis range

% Plot it!
figure(2); clf; set(gcf, 'Color', 'w', 'Position', [ 394   156   836   649])
if strcmp(fitTypeName, 'linear')
    plot(xThresh, 10.^y, 'r-', 'LineWidth', 3); hold all;
    if ~isempty(varThresh)
        errorbar(xThresh, yThresh, varThresh, 'Color', 'k', 'LineStyle','none', 'LineWidth', 2);
    end
   scatter(xThresh, yThresh, 80, 'MarkerFaceColor', 'w', 'MarkerEdgeColor','k', 'LineWidth',2);

end


box off;
set(gca, 'TickDir', 'out','TickLength',[0.015 0.015], 'LineWidth',1,'Fontsize',20,'XScale','log', 'YScale', yScale)
xlabel('Cone density (cells/deg^2)','FontSize',20); ylabel('Contrast threshold (%)','FontSize',25)
set(gca, 'XTick',10.^xticks,'XTickLabel',xticklbls, 'XLim', 10.^[2 5]);
if strcmp(yScale, 'log')
   yrange = [0.01 0.1];
   set(gca,'YLim', [0.005 0.2], 'YTick',yrange,'YTickLabel',sprintfc('%1.2f',yrange));
   set(gca,'XGrid','on', 'YGrid','on','XMinorGrid','on','YMinorGrid','on', ...
       'GridAlpha',0.25, 'LineWidth',0.5); drawnow;
else
    set(gca, 'YTick',[0, 0.05, 0.1],'YTickLabel',[0 5 10], 'YLim', [0 0.1]);
end
legend off; title(sprintf('Contrast threshold vs level of cone density - R^2: %1.2f', R2))

% Save fig if requested
if saveFig
    if ~exist(figurePth,'dir'); mkdir(figurePth); end
    savefig(fullfile(figurePth,sprintf('contrastThresholdVS%s_%s_%s',expName, fitTypeName, yScale)))
    hgexport(gcf,fullfile(figurePth,sprintf('contrastThresholdVS%s_%s_%s',expName,fitTypeName, yScale)))
    print(fullfile(figurePth,sprintf('contrastThresholdVS%s_%s_%s',expName,fitTypeName, yScale)), '-dpng')
end

%% 3. Quantify effect of cone density

% contrast thresholds reported in behavior, corresponding to horizontal versus upper vertical meridian
% reportedBehavior = [0.02+0.015, 0.02]; % (thresholds in % contrast)

ang = [0, 90, 180, 270]; % polar angles: nasal (HM), superior (LVM), temporal (HM), inferior (UVM) (radians)
for jj = 1:length(ang)
    coneDensityMM2 = coneDensityReadData('coneDensitySource', 'Curcio1990','eccentricity',4.5,'angle',ang(jj),'eccentricityUnits', 'deg','whichEye','left');
    coneDensityDeg2(jj) = coneDensityMM2/((1/.3)^2);
end


if strcmp(fitTypeName, 'linear')
    % Get intercept and slope of log-linear fit
    b_intcpt = fitResult(2);
    a_coeff  = fitResult(1);
    
    % Reconstruct log-linear function
    cThreshold = @(x) 10.^(a_coeff* log10(x) + b_intcpt);
    predictedDensity = @(y) 10.^((y-b_intcpt)./a_coeff);
    modelPredictionForPF = cThreshold(coneDensityDeg2);
    
    errorRatioConeDensity = 2*abs(coneDensityDeg2-mean(coneDensityDeg2));

    % Nasal retina
    predictedError = [];
    predictedError(1,1) = cThreshold(coneDensityDeg2(1)-errorRatioConeDensity(1));
    predictedError(1,2) = cThreshold(coneDensityDeg2(1)+errorRatioConeDensity(1));

    % Superior retina
    predictedError(2,1) = cThreshold(coneDensityDeg2(2)-errorRatioConeDensity(2));
    predictedError(2,2) = cThreshold(coneDensityDeg2(2)+errorRatioConeDensity(2));

    % Temporal retina
    predictedError(3,1) = cThreshold(coneDensityDeg2(3)-errorRatioConeDensity(3));
    predictedError(3,2) = cThreshold(coneDensityDeg2(3)+errorRatioConeDensity(3));

    % Inferior retina
    predictedError(4,1) = cThreshold(coneDensityDeg2(4)-errorRatioConeDensity(4));
    predictedError(4,2) = cThreshold(coneDensityDeg2(4)+errorRatioConeDensity(4));

    save(fullfile(ogRootPath, 'data', 'conedensity','coneOnly_predictedMeanAndError_stimeccen'),'modelPredictionForPF','predictedError', 'coneDensityDeg2', 'ang')
elseif strcmp(fitTypeName, 'poly2')
    
    modelPredictionForPF = polyval(fitResult,reportedBehavior);

end


totalVariance.densityPredictedByModelHVA = 100.*diff([mean(modelPredictionForPF([1 3])),mean(modelPredictionForPF([2 4]))]) ./ mean(modelPredictionForPF);
totalVariance.densityPredictedByModelVMA = 100.*diff([modelPredictionForPF(4), modelPredictionForPF(2)]) ./ mean(modelPredictionForPF([2 4]));

fprintf('HVA contrast thresholds predicted by computational observer model given cone density: %2.1f percent\n', totalVariance.densityPredictedByModelHVA)
fprintf('VMA contrast thresholds predicted by computational observer model given cone density: %2.1f percent\n', totalVariance.densityPredictedByModelVMA)

end
