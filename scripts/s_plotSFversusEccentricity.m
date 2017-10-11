%% s_plotSFversusEccentricity

eyemovement         = {'110'}; % No eyemovement
defocusZ            = 0;      % No defocus

colors              = copper(11);

eccentricities      = [0 2 5 10 20 40];      % deg
contrastLevels      = [0:0.01:0.1]; % Contrast levels of stimulus used in simulation

spatFreqs           = [0.25, 0.4, 0.65, 1, 1.6, 2.6, 4, 8, 10, 16, 26];

usedLabels          = contrastLevels;

labels              = {'0.25','', '',...
    '0.4', '', '',...
    '0.65', '', '',...
    '1', '', '',...
    '1.6', '', '',...
    '2.6','', '',...
    '4', '', '',...
    '8', '', '',...
    '10', '', '',...
    '16', '', '',...
    '26'}; % labels when plotting psychometric functions

polarAngles = 0;
FFTflag     = true;

% Where to find data
dataPth     = fullfile(ogRootPath,'figs');

% Number of total trials in computational observer model (50 clockwise, 50 counterclockwise)
nTotal      = 100;

% Prepare fit variables
fit = [];

fit.ctrpred = cell(size(colors,1),1);
fit.ctrvar  = cell(size(colors,1),1);
fit.ctrr2   = cell(size(colors,1),1);
fit.data    = cell(size(colors,1),1);

fit.init = [0.5, 0.1];
fit.thresh = 0.75;

count = 1;


for eccen = eccentricities
    
    %% 1. Load results
    fName   = sprintf('Classify_coneOutputs_contrast0.10_pa0_eye110_eccen%1.2f_defocus0.00_noise-random_phasescrambled_fft%d.mat', eccen,FFTflag);
    
    accuracy = load(fullfile(dataPth, 'Eccen_SF', fName));
    accuracy.P = squeeze(accuracy.P);
    if size(accuracy.P,1)<size(accuracy.P,2)
        accuracy.P = accuracy.P';
    end
    
    for ii = 1:size(accuracy.P)
        
        %% 2. Fit Weibull
        % Make a Weibull function first with contrast levels and then search for
        % the best fit with the classifier data
        fit.ctrvar{ii} = fminsearch(@(x) ogFitWeibull(x, usedLabels, accuracy.P(:,ii), nTotal), fit.init);
        
        % Then fit a Weibull function again, but now with the best fit parameters
        % from the previous step.
        fit.ctrpred{ii} = ogWeibull(fit.ctrvar{ii}, usedLabels);
        
        % Not sure what this line is for..
        % fit.ctrr2 = corr(fit.ctrpred', dec.ctr).^2;
        
        %% 3. Find contrast threshold
        diff   = abs(fit.ctrpred{ii} - fit.thresh);
        minval = find(diff == min(diff));
        fit.ctrthresh{ii} = usedLabels(minval(1));
        fit.data{ii} = accuracy.P(:,ii);
        
        thresholds(eccen==eccentricities,ii) = usedLabels(minval(1));
        
    end
    
    %% 4. Visualize
    
    figure; set(gcf,'Color','w'); hold all;
    
    for ii = 1:length(fit.ctrpred)
        dataToFit = squeeze(fit.data{ii})';
        plot(usedLabels, fit.ctrpred{ii}*100, 'Color', colors(ii,:), 'LineWidth',2);
        scatter(usedLabels, dataToFit, 80, colors(ii,:), 'filled');
        plot(10.^-2.1,dataToFit(1),'o','Color',colors(ii,:), 'MarkerSize', 8, 'MarkerFaceColor',colors(ii,:))
    end
    
    set(gca, 'XScale','log', 'XLim',[min(usedLabels) max(usedLabels)],'YLim', [min(dataToFit)-10 100], 'TickDir','out','TickLength',[.015 .015],'FontSize',17, 'LineWidth',2);
    set(gca, 'XTick', [0.01,0.03,0.05, 0.1], 'XTickLabel',{'1','3','5','10'})
    
    ylabel('Classifier Accuracy (% Correct)', 'FontSize',17)
    % if strcmp(whatPlot,'ConeDensity')
    %     xlabel('Eccentricity (deg)', 'FontSize',17);
    %     set(gca,'XScale','linear')
    xlabel('Stimulus Contrast (%)', 'FontSize',17);
    
    legend(labels, 'Location','bestoutside'); legend boxoff
    title(sprintf('Eccentricity: %1.2f',eccen))
    
    
end


figure; hold on;
for ii = 1:length(eccentricities)
    plot(spatFreqs,1./thresholds(ii,:), 'o-','LineWidth',4)
end

xlim([0.1 100]); ylim([1 100])
set(gca,'XScale','log','YScale','Log','TickDir','out', 'FontSize',17)

legend('Eccen: 0 deg','Eccen: 2 deg','Eccen: 5 deg','Eccen: 10 deg','Eccen: 20 deg','Eccen: 40 deg', 'Location','Best'); legend boxoff; 
xlabel('Spatial frequency (cpd)','FontSize',17)
ylabel('Contrast sensitivity (1/threshold)','FontSize',17)




%
% savefig(fullfile(dataPth,sprintf('WeibullFit_contrastVSperformance_all_pa%d_fft%d_%s_noPCA',polarAngles,FFTflag,whatPlot)))
% hgexport(gcf,fullfile(dataPth,sprintf('WeibullFit_contrastVSperformance_all_pa%d_fft%d_%s_noPCA.eps',polarAngles,FFTflag,whatPlot)))
