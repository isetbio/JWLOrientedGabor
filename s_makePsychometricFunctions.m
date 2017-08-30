%% s_makePsychometricFunctions


% Script to compute psychometric functions based on the computational
% observer model 

%% 0. Load results

eyemovement = {'100'};
polarAngles = 0;
FFTflag = false;

dataPth = fullfile(ogRootPath,'data');
fName   = sprintf('contrastVSperformance_eye%s_pa%d_fft%d.mat',cell2mat(eyemovement),polarAngles,FFTflag);

accuracy = load(fullfile(dataPth, fName));

%% 1. Set parameters

contrastLevels    = [0.01:0.01:0.09, 0.1:0.1:1.0]; % Contrast levels of stimulus used in simulation
% nCLevels           = length(Contrastlevels);     % Number of contrast levels
% nAccuracy         = length(accuracy.P);          % Actual data (accuracy of linear classifier)
nTotal            = 100;                           % Number of total trials in computational observer model (50 clockwise, 50 counterclockwise)

%% 2. Fit Weibull

% Prepare fit variables
fit = [];

fit.ctrpred = []; 
fit.ctrvar  = []; 
fit.ctrr2   = [];

fit.init = [0.5, 0.1];
fit.thresh = 0.75;

% Make a Weibull function first with contrast levels and then search for
% the best fit with the classifier data
fit.ctrvar = fminsearch(@(x) ogFitWeibull(x, contrastLevels, accuracy.P', nTotal), fit.init);

% Then fit a Weibull function again, but now with the best fit parameters
% from the previous step.
fit.ctrpred = ogWeibull(fit.ctrvar, contrastLevels);

% Not sure what this line is for..
% fit.ctrr2 = corr(fit.ctrpred', dec.ctr).^2;

%% 3. Find contrast threshold
diff   = abs(fit.ctrpred - fit.thresh);
minval = find(diff == min(diff));
fit.ctrthresh = contrastLevels(minval(1));

%% 4. Visualize

% colors = lines(length(eyemovement));

figure(1); clf; set(gcf,'Color','w'); 
plot(contrastLevels, fit.ctrpred*100, 'Color', 'k', 'LineWidth',2);
set(gca, 'XScale','log', 'YLim', [40 100], 'TickDir','out','TickLength',[.015 .015],'FontSize',12);
ylabel('Classifier Accuracy')
xlabel('Contrast level (Michelson)')
title('Linear SVM, fixed eye', 'FontSize',12)
box off

savefig(fullfile(ogRootPath,'figs',sprintf('WeibullFit_contrastVSperformance_eye%s_pa%d_fft%d',cell2mat(eyemovement),polarAngles,FFTflag)))
hgexport(gcf,fullfile(ogRootPath,'figs',sprintf('WeibullFit_contrastVSperformance_eye%s_pa%d_fft%d.eps',cell2mat(eyemovement),polarAngles,FFTflag)))

