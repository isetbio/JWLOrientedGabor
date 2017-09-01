%% s_makePsychometricFunctions


% Script to compute psychometric functions based on the computational
% observer model


%% 0. Set parameters
eyemovement = {'000','100','010','110'};
polarAngles = 0;
FFTflag = true;

dataPth = fullfile(ogRootPath,'figs');

contrastLevels    = [0.01:0.01:0.09, 0.1:0.1:1.0]; % Contrast levels of stimulus used in simulation
eyemovement       = {'110'};
% nCLevels           = length(Contrastlevels);     % Number of contrast levels
% nAccuracy         = length(accuracy.P);          % Actual data (accuracy of linear classifier)
nTotal            = 100;                           % Number of total trials in computational observer model (50 clockwise, 50 counterclockwise)
    
% Prepare fit variables
fit = [];

fit.ctrpred = cell(length(eyemovement),1);
fit.ctrvar  = cell(length(eyemovement),1);
fit.ctrr2   = cell(length(eyemovement),1);

fit.init = [0.5, 0.1];
fit.thresh = 0.75;

coneTypes = {'','_L','_M','_S'};

for em = 1:length(eyemovement)
    for ct = 1:length(coneTypes)     
    
        %% 1. Load results
    fName   = sprintf('contrastVSperformance_eye%s_pa%d_fft%d%s.mat',cell2mat(eyemovement(em)),polarAngles,FFTflag, coneTypes{ct});
    accuracy = load(fullfile(dataPth, fName));
    
    
    %% 2. Fit Weibull
    % Make a Weibull function first with contrast levels and then search for
    % the best fit with the classifier data
    fit.ctrvar{ct} = fminsearch(@(x) ogFitWeibull(x, contrastLevels, accuracy.P', nTotal), fit.init);
    
    % Then fit a Weibull function again, but now with the best fit parameters
    % from the previous step.
    fit.ctrpred{ct} = ogWeibull(fit.ctrvar{ct}, contrastLevels);
    
    % Not sure what this line is for..
    % fit.ctrr2 = corr(fit.ctrpred', dec.ctr).^2;
    
    %% 3. Find contrast threshold
    diff   = abs(fit.ctrpred{ct} - fit.thresh);
    minval = find(diff == min(diff));
    fit.ctrthresh{ct} = contrastLevels(minval(1));
    fit.data{ct} = accuracy.P;
    
    end
end

%% 4. Visualize

% colors = lines(length(eyemovement));

figure(1); clf; set(gcf,'Color','w'); hold all;

plot(contrastLevels, fit.ctrpred{1}*100, 'Color', 'k', 'LineWidth',2);
scatter(contrastLevels, fit.data{1}, [], 'k');
plot(contrastLevels, fit.ctrpred{2}*100, 'Color', 'b', 'LineWidth',2);
scatter(contrastLevels, fit.data{2},[], 'b');
plot(contrastLevels, fit.ctrpred{3}*100, 'Color', 'g', 'LineWidth',2);
scatter(contrastLevels, fit.data{3}, [], 'g');
plot(contrastLevels, fit.ctrpred{4}*100, 'Color', 'r', 'LineWidth',2);
scatter(contrastLevels, fit.data{4}, [], 'r');

set(gca, 'XScale','log', 'YLim', [0 100], 'TickDir','out','TickLength',[.015 .015],'FontSize',12);
ylabel('Classifier Accuracy')
xlabel('Contrast level (Michelson)')
title('Linear SVM, Ecc 6, Polar Angle 0, 6cpd, monochromatics versus LMS','FontSize',12)
% legend('Fixed Fit','Fixed Data','Tremor Fit','Tremor Data','Drift Fit','Drift Data','Tremor+Drift Fit', 'Tremor+Drift Data','Location','Best')
legend('LMS Fit','LMS Data','L Fit','L Data','M Fit','M Data','S Fit', 'S Data','Location','Best')

box off

savefig(fullfile(ogRootPath,'figs',sprintf('WeibullFit_contrastVSperformance_all_pa%d_fft%d_monochromatic',polarAngles,FFTflag)))
hgexport(gcf,fullfile(ogRootPath,'figs',sprintf('WeibullFit_contrastVSperformance_all_pa%d_fft%d_monochromatic.eps',polarAngles,FFTflag)))

