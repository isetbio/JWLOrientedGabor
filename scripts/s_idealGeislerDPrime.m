%% s_idealGeislerDPrime

% save figures?
saveFigs = false;
figure(1); clf; hold all;

% Set arbitrary zero point for log-linear plot
logzero = 5e-5;

 %% 1. Analytical solution
 
 % Load parameters
expNameTemplate = 'idealobserver';
expParams       = loadExpParams(expNameTemplate, false);

% Get percent correct
percentCorrectAnalytic = geislerIdealAnalytical(expParams);

% Plot it
figure(1); plot(expParams.contrastLevels(2:end), percentCorrectAnalytic(2:end), 'ko-'); hold on;
plot(logzero, percentCorrectAnalytic(1), 'ko')
set(gca, 'XScale', 'log', 'YLim', [40, 100])

%% 2. Simulation of analytic solution

% Load data parameters
expNameData     = 'defaultnophaseshift';
expParams       = loadExpParams(expNameData, false);
subFolderName   = 'run1';

% Get percent correct
percentCorrectSVMDiffTemplate = geislerIdealSimulation(expParams, expNameTemplate, expNameData, subFolderName);

% Plot it
figure(1); plot(expParams.contrastLevels(2:end), percentCorrectSVMDiffTemplate(2:end), 'ko-'); hold on;
plot(logzero, percentCorrectSVMDiffTemplate(1), 'ko')
set(gca, 'XScale', 'log', 'YLim', [.4, 1])

%% 3. SVM with data sample match to template

% Load data parameters
expNameData     = 'defaultnophaseshift';
expParams       = loadExpParams(expNameData, false);
expNameTemplate = 'idealobserver';
subFolderName   = 'run1';

% Get percent correct
percentCorrectSVMDiffTemplate = svmTemplateSimulation(expParams, expNameTemplate, expNameData, subFolderName);

% Plot simulation results
figure(1); plot(expParams.contrastLevels(2:end), percentCorrectSVMDiffTemplate(2:end), 'ro-'); hold on;
plot(logzero, percentCorrectSVMDiffTemplate(1), 'ro')
set(gca, 'XScale', 'log', 'YLim', [.40 1])





% OLD CODE 
%     % Compute mean and varianec across trials
%     meanZGivenCCW = mean(ZGivenCCW);
%     meanZGivenCW = mean(ZGivenCW);
%     varZGivenCCW = var(ZGivenCCW);
%     varZGivenCW = var(ZGivenCW);
%
%     % Compute d-prime
%     dprime = (meanZGivenCCW-meanZGivenCW) / sqrt(0.5*(varZGivenCCW+varZGivenCW));
%
%
%     % Plot loglikelihood distributions of the two stimulus classes
%     figure(2); clf; set(gcf, 'Color', 'w', 'Position', [508, 820, 1052, 518]);
%     subplot(1,2,1); hold all;
%     [nA,xA] = hist(ZGivenCW,40);
%     bar(xA,nA, 'EdgeColor', 'k', 'FaceColor', 'w')
%     predA = length(ZGivenCW)*normpdf(xA,meanZGivenCW,sqrt(varZGivenCW))*(xA(2)-xA(1));
%     plot(xA,predA,'g','LineWidth',4);
%     set(gca, 'TickDir', 'out', 'FontSize', 14);
%
%     [nB,xB] = hist(ZGivenCCW,40);
%     bar(xB,nB, 'EdgeColor', 'k', 'FaceColor', 'k');
%     predB = length(ZGivenCCW)*normpdf(xB,meanZGivenCCW,sqrt(varZGivenCCW))*(xB(2)-xB(1));
%     plot(xB,predB,'r','LineWidth',4);
%     title(sprintf('Probability of absorptions for Alpha vs Beta stim, contrast: %1.4f',c))
%     legend({'Alpha samples', 'Alpha fit', 'Beta samples', 'Beta fit'}, 'FontSize', 14, 'Location', 'best'); legend boxoff;
%
%     % Get ROC curve from Hit and False Alarm rates and sample along many criteria
%     criteria = linspace(min([ZGivenCW,ZGivenCCW]),max([ZGivenCW,ZGivenCCW]),1000);
%     for cc = 1:length(criteria)
%         HitRate(cc) = length(find(ZGivenCCW > criteria(cc))) / length(ZGivenCCW);
%         FARate(cc) = length(find(ZGivenCW > criteria(cc))) / length(ZGivenCCW);
%     end
%
%     % Get the percent correct
%     percentCorrectSim2(c==expParams.contrastLevels) = -trapz([1 FARate 0],[1 HitRate 0]);
%     subplot(1,2,2);
%     plot([1 FARate 0],[1 HitRate 0],'r','LineWidth',4);
%     xlabel('False Alarm Rate');
%     ylabel('Hit Rate');
%     title(sprintf('ROC Curve, d-prime: %2.2f, percent corrent: %2.2f',dprime,percentCorrectSim2(c==expParams.contrastLevels)*100));
%     axis('square'); axis([0 1 0 1]); set(gca, 'TickDir', 'out', 'FontSize', 14);
%     drawnow;
%     if saveFigs
%         savefig(fullfile(figurePth,sprintf('Geisler_Ideal_dprime_%s_c%1.4f.fig',expNameTemplate,c)))
%         hgexport(gcf,fullfile(figurePth,sprintf('Geisler_Ideal_dprime_%s_c%1.4f.eps',expNameTemplate,c)))
%     end