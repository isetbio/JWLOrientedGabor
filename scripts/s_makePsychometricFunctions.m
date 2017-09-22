%% s_makePsychometricFunctions


% Script to compute psychometric functions based on the computational
% observer model

whatPlot = 'ConeDensity';

switch whatPlot
    case 'default'
<<<<<<< HEAD
        eyemovement         = {'000'}; % No eyemovement
=======
        eyemovement         = {'110'}; % No eyemovement
>>>>>>> 8723303378da234edfd736e95317319491a54833
        defocusZ            = 0;      % No defocus
        coneType            = {''};    % All LMS cones
        colors              = [0 0 0];
        labels              = {'Fit','data'};
        eccentricities      = 6;      % deg
        contrastLevels      = [0:0.01:0.09, 0.1:0.1:1.0]; % Contrast levels of stimulus used in simulation
        usedLabels          = contrastLevels;
        
    case 'ConeTypes'
        
        eyemovement         = {'110'}; % Drift, tremor, no ms
        defocusZ            = 0;
        coneType            = {'','_L', '_M','_S'};
        colors              = {'k', 'r','g','b'};
        labels              = {'LMS cones','', ...
                                'L cones','', ...
                                'M cones','', ...
                                'S cones',''};
        eccentricities      = 6;
        contrastLevels      = [0:0.01:0.09, 0.1:0.1:1.0]; % Contrast levels of stimulus used in simulation
        usedLabels          = contrastLevels;
        
        
        
    case 'EyeMovement'
        eyemovement         = {'000','100','010','110'}; % No, Drift, tremor, Drift + tremor
        defocusZ            = 0;
        coneType            = {''};
        colors              = copper(5); %{'k', 'r','g','b'};
        labels              = {'Fixed','','', ...
                                'Tremor','','', ...
                                'Drift','', '',...
                                'Tremor+Drift','',''};
        eccentricities      = 6;
<<<<<<< HEAD
        contrastLevels      = [0 0.01:0.01:0.09, 0.1:0.1:1.0]; % Contrast levels of stimulus used in simulation
=======
        contrastLevels      = [0:0.01:0.09, 0.1:0.1:1.0]; % Contrast levels of stimulus used in simulation
>>>>>>> 8723303378da234edfd736e95317319491a54833
        usedLabels          = contrastLevels;
        
    case 'EyeMovementEnhanced'
        eyemovement         = {'000','110','220','330'}; % No, Drift, tremor, Drift + tremor
        defocusZ            = 0;
        coneType            = {''};
        colors              = copper(4); %{'k', 'r','g','b'};
        labels              = {'Fixed','', '', ...
                                'Default Tremor + Drift','', '', ...
                                '2x Tremor + Drift','', '',...
                                '3x Tremor + Drift','',''};
        eccentricities      = 6;
<<<<<<< HEAD
        contrastLevels      = [0 0.01:0.01:0.09, 0.1:0.1:1.0]; % Contrast levels of stimulus used in simulation
=======
        contrastLevels      = [0:0.01:0.09, 0.1:0.1:1.0]; % Contrast levels of stimulus used in simulation
>>>>>>> 8723303378da234edfd736e95317319491a54833
        usedLabels          = contrastLevels;
        
    case 'ConeDensity'
        eyemovement         = {'110'};
        defocusZ            = 0;
        coneType            = {''};
        colors              = [0 0 0];
        labels              = {'Fit','Data'};
        eccentricities      = 40;
        contrastLevels      = 0.04;
        
        %% Change x labels to density
        whichEye          = 'left';
        cparams.cmFOV     =  2; % degrees

        % Convertion deg to m
        deg2m  = 1/3 * 0.001; % 3 deg per mm, .001 mm per meter

        % Predefine density vector
        allDensity = nan(eccentricities,1);

        for eccen = 1:eccentricities
            % Specify retinal location where stimulus is presented
            cparams.eccentricity = eccen;             % Visual angle of stimulus center, in deg
            cparams.polarAngle   = deg2rad(0);   % Polar angle (radians): 0 is right, pi/2 is superior, pi is left, 3*pi/2 inferior

            % Compute x,y position in m of center of retinal patch from ecc and angle
            [x, y] = pol2cart(cparams.polarAngle, cparams.eccentricity);
            x = x * deg2m;  y = y * deg2m;
            cMosaic = coneMosaic('center', [x, y], 'whichEye', whichEye);

            % Set the field of view (degrees)
            cMosaic.setSizeToFOV(cparams.cmFOV);
            allDensity(eccen,:) = eccen2density(cMosaic, 'm');
        end
        
<<<<<<< HEAD
        usedLabels          = 1:40;
=======
        usedLabels          = 1:eccentricities; % For now use eccentricities as labels, but we could plot it against cone density
>>>>>>> 8723303378da234edfd736e95317319491a54833
        
    case 'Defocus'
        eyemovement         = {'110'};
        defocusZ            = [0:0.5:2];
        coneType            = {''};
        colors              = copper(5);%{'k','r','g','b','c'};
<<<<<<< HEAD
        labels              = {'0 Defocus','','','0.5 Defocus','','','1.0 Defocus','','','1.5 Defocus','','','2.0 Defocus','',''};
        eccentricities      = 6;
        contrastLevels      = [0 0.01:0.01:0.09, 0.1:0.1:1.0];
=======
        labels              = {'0 Defocus','','', ...
                                '0.5 Defocus','','', ...
                                '1.0 Defocus','','',...
                                '1.5 Defocus','','',...
                                '2.0 Defocus','',''};
        eccentricities      = 6;
        contrastLevels      = [0:0.01:0.09, 0.1:0.1:1.0];
>>>>>>> 8723303378da234edfd736e95317319491a54833
        usedLabels          = contrastLevels;
        
        
end



%% 0. Set other parameters
polarAngles = 0;
FFTflag     = true;

% Where to find data
dataPth     = fullfile(ogRootPath,'figs');

% Number of total trials in computational observer model (50 clockwise, 50 counterclockwise)
nTotal      = 100;

% Prepare fit variables
fit = [];

fit.ctrpred = cell(size(colors,1));
fit.ctrvar  = cell(size(colors,1));
fit.ctrr2   = cell(size(colors,1));
fit.data    = cell(size(colors,1));

fit.init = [0.5, 0.1];
fit.thresh = 0.75;

count = 1;
for em = 1:length(eyemovement)
    for ct = 1:length(coneType)
        for eccen = eccentricities
            for df = 1:length(defocusZ)
                
                %% 1. Load results
                fName   = sprintf('contrastVSperformance_eye%s_pa%d_fft%d_eccen%1.2f%s_defocus%1.2f_pca0.mat', ...
                    cell2mat(eyemovement(em)),polarAngles,FFTflag,eccen,coneType{ct},defocusZ(df));
                
                
                accuracy = load(fullfile(dataPth, fName));
                accuracy.P = squeeze(accuracy.P);
                if size(accuracy.P,1)<size(accuracy.P,2)
                    accuracy.P = accuracy.P';
                end
                
                %% 2. Fit Weibull
                % Make a Weibull function first with contrast levels and then search for
                % the best fit with the classifier data
                fit.ctrvar{count} = fminsearch(@(x) ogFitWeibull(x, usedLabels, accuracy.P, nTotal), fit.init);
                
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

figure(1); clf; set(gcf,'Color','w'); hold all;

for ii = 1:length(fit.ctrpred)
    dataToFit = squeeze(fit.data{ii});
    if ~strcmp(whatPlot,'ConeDensity'); plot(usedLabels, fit.ctrpred{ii}*100, 'Color', colors(ii,:), 'LineWidth',2); end
    scatter(usedLabels, dataToFit, 80, colors(ii,:), 'filled');
    plot(10.^-2.1,dataToFit(1),'o','Color',colors(ii,:), 'MarkerSize', 8, 'MarkerFaceColor',colors(ii,:))
end

set(gca, 'XScale','log', 'XLim',[min(usedLabels) max(usedLabels)],'YLim', [40 100], 'TickDir','out','TickLength',[.015 .015],'FontSize',12, 'LineWidth',2);

ylabel('Classifier Accuracy', 'FontSize',12)
if strcmp(whatPlot,'ConeDensity'); 
    xlabel('Eccentricity (deg)', 'FontSize',12); 
    set(gca,'XScale','linear')
else xlabel('Contrast (Michelson)', 'FontSize',12); end

title(sprintf('Linear SVM, Ecc %d, Polar Angle %d, 6cpd, FFT%d',max(eccentricities),polarAngles,FFTflag),'FontSize',12)
legend(labels, 'Location','Best', 'box','off')

box off

savefig(fullfile(dataPth,sprintf('WeibullFit_contrastVSperformance_all_pa%d_fft%d_%s_noPCA',polarAngles,FFTflag,whatPlot)))
hgexport(gcf,fullfile(dataPth,sprintf('WeibullFit_contrastVSperformance_all_pa%d_fft%d_%s_noPCA.eps',polarAngles,FFTflag,whatPlot)))

