%% s_makePsychometricFunctions


% Script to compute psychometric functions based on the computational
% observer model

whatPlot = 'Defocus';

switch whatPlot
    case 'default'
        eyemovement         = {'000'}; % No eyemovement
        defocusZ            = [];      % No defocus
        coneType            = {''};    % All LMS cones
        colors              = {'k'};
        labels              = {'Fit','data'};
        eccentricities      = 6;      % deg
        contrastLevels      = [0.01:0.01:0.09, 0.1:0.1:1.0]; % Contrast levels of stimulus used in simulation
        usedLabels          = contrastLevels;
        
    case 'ConeTypes'
        
        eyemovement         = {'110'}; % Drift, tremor, no ms
        defocusZ            = [];
        coneType            = {'','_L', '_M','_S'};
        colors              = {'k', 'r','g','b'};
        labels              = {'LMS cones','', ...
            'L cones','', ...
            'M cones','', ...
            'S cones', ''};
        eccentricities      = 6;
        contrastLevels      = [0.01:0.01:0.09, 0.1:0.1:1.0]; % Contrast levels of stimulus used in simulation
        usedLabels          = contrastLevels;
        
        
        
    case 'EyeMovement'
        eyemovement         = {'000','100','010','110'}; % No, Drift, tremor, Drift + tremor
        defocusZ            = [];
        coneType            = {''};
        colors              = {'k', 'r','g','b'};
        labels              = {'Fixed','', ...
            'Tremor','', ...
            'Drift','', ...
            'Tremor+Drift', ''};
        eccentricities      = 6;
        contrastLevels      = [0.01:0.01:0.09, 0.1:0.1:1.0]; % Contrast levels of stimulus used in simulation
        usedLabels          = contrastLevels;
        
    case 'ConeDensity'
        eyemovement         = {'110'};
        defocusZ            = [];
        coneType            = {''};
        colors              = {'k'};
        labels              = {'Fit','Data'};
        eccentricities      = 40;
        contrastLevels      = 0.08;
        usedLabels          = 1:eccentricities;
        
    case 'Defocus'
        eyemovement         = {'110'};
        defocusZ            = [0:0.5:2];
        coneType            = {''};
        colors              = {'k','r','g','b','c'};
        labels              = {'0 Defocus','','0.5 Defocus','','1.0 Defocus','','1.5 Defocus','','2.0 Defocus',''};
        eccentricities      = 6;
        contrastLevels      = [0.01:0.01:0.09, 0.1:0.1:1.0];
        usedLabels          = contrastLevels;
        
        
end



%% 0. Set other parameters
polarAngles = 0;
FFTflag     = false;

% Where to find data
dataPth     = fullfile(ogRootPath,'figs');

% Number of total trials in computational observer model (50 clockwise, 50 counterclockwise)
nTotal      = 100;

% Prepare fit variables
fit = [];

fit.ctrpred = cell(size(colors));
fit.ctrvar  = cell(size(colors));
fit.ctrr2   = cell(size(colors));

fit.init = [0.5, 0.1];
fit.thresh = 0.75;

count = 1;
for em = 1:length(eyemovement)
    for ct = 1:length(coneType)
        for eccen = eccentricities
            for df = 1:length(defocusZ)
                
                %% 1. Load results
                fName   = sprintf('contrastVSperformance_eye%s_pa%d_fft%d_eccen%1.2f%s_defocus%1.2f.mat', ...
                    cell2mat(eyemovement(em)),polarAngles,FFTflag,eccen,coneType{ct,defocusZ(df));
                
                
                accuracy = load(fullfile(dataPth, fName));
                
                
                %% 2. Fit Weibull
                % Make a Weibull function first with contrast levels and then search for
                % the best fit with the classifier data
                fit.ctrvar{count} = fminsearch(@(x) ogFitWeibull(x, usedLabels, (accuracy.P(:,:,:,:,df))', nTotal), fit.init);
                
                % Then fit a Weibull function again, but now with the best fit parameters
                % from the previous step.
                fit.ctrpred{count} = ogWeibull(fit.ctrvar{count}, usedLabels);
                                
                % Not sure what this line is for..
                % fit.ctrr2 = corr(fit.ctrpred', dec.ctr).^2;
                
                %% 3. Find contrast threshold
                diff   = abs(fit.ctrpred{count} - fit.thresh);
                minval = find(diff == min(diff));
                fit.ctrthresh{count} = usedLabels(minval(1));
                fit.data{count} = accuracy.P;
                
                count = count +1;
            end
            
        end
    end
end

%% 4. Visualize

% colors = lines(length(eyemovement));

figure(1); clf; set(gcf,'Color','w'); hold all;

dataToFit = squeeze(fit.data{1});

for ii = 1:length(fit.ctrpred)
    plot(usedLabels, fit.ctrpred{ii}*100, 'Color', colors{ii}, 'LineWidth',2);
    scatter(usedLabels, dataToFit(:,ii), [], colors{ii});
    
end

set(gca, 'XScale','log', 'YLim', [40 100], 'TickDir','out','TickLength',[.015 .015],'FontSize',12);
ylabel('Classifier Accuracy')
if strcmp(whatPlot,'Density'); xlabel('Eccentricity (deg)'); else xlabel('Contrast (Michelson)'); end
title(sprintf('Linear SVM, Ecc %d, Polar Angle %d, 6cpd, FFT%d',max(eccentricities),polarAngles,FFTflag),'FontSize',12)
legend(labels, 'Location','Best')


box off

savefig(fullfile(ogRootPath,'figs',sprintf('WeibullFit_contrastVSperformance_all_pa%d_fft%d_%s',polarAngles,FFTflag,whatPlot)))
hgexport(gcf,fullfile(ogRootPath,'figs',sprintf('WeibullFit_contrastVSperformance_all_pa%d_fft%d_%s.eps',polarAngles,FFTflag,whatPlot)))

