%% s_makePsychometricFunctions


% Script to compute psychometric functions based on the computational
% observer model

whatPlot = 'EyeMovement';
currentFlag = false;

switch whatPlot
    case 'default'
        eyemovement         = {'110'}; % No eyemovement
        defocusZ            = 0;      % No defocus
        coneType            = {''};    % All LMS cones
        colors              = [0 0 0];
        labels              = {'Fit','data'};
        eccentricities      = 4.5;      % deg
        contrastLevels      = [0:0.01:0.1];%, 0.1:0.1:1.0]; % Contrast levels of stimulus used in simulation
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
        eyemovement         = {'000', '100', '010', '001','110'}; % No, Drift, tremor, Drift + tremor, Drift + tremor + MS
        defocusZ            = 0;
        coneType            = {''};
        colors              = copper(5); %{'k', 'r','g','b'};
        labels              = {'Fixed','','', ...
                                'Tremor','','', ...
                                'Drift','', '',...
                                'MicroSaccades','', '',...
                                'Tremor+Drift','',''};
        eccentricities      = 4.5;
        contrastLevels      = [0:0.01:0.1 0.2]; % Contrast levels of stimulus used in simulation
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
        eccentricities      = 4.5;
        contrastLevels      = [0:0.01:0.1];%, 0.1:0.1:1.0]; % Contrast levels of stimulus used in simulation
        usedLabels          = contrastLevels;
        
    case 'ConeDensity'
        eyemovement         = {'110'};
        defocusZ            = 0;
        coneType            = {''};
        eccentricities      = [0,4.5,16,60];
        colors              = copper(length(eccentricities));

        contrastLevels      = [0:0.01:0.1];
        
        %% Change x labels to density
        whichEye          = 'left';
        cparams.cmFOV     =  1; % degrees

        % Convertion deg to m
        deg2m  = 1/3 * 0.001; % 3 deg per mm, .001 mm per meter

        % Predefine density vector
        allDensity = nan(length(eccentricities),1);

        for eccen = eccentricities
            % Specify retinal location where stimulus is presented
            cparams.eccentricity = eccen;             % Visual angle of stimulus center, in deg
            cparams.polarAngle   = deg2rad(0);   % Polar angle (radians): 0 is right, pi/2 is superior, pi is left, 3*pi/2 inferior

            % Compute x,y position in m of center of retinal patch from ecc and angle
            [x, y] = pol2cart(cparams.polarAngle, cparams.eccentricity);
            x = x * deg2m;  y = y * deg2m;
            cMosaic = coneMosaic('center', [x, y], 'whichEye', whichEye);

            % Set the field of view (degrees)
            cMosaic.setSizeToFOV(cparams.cmFOV);
            allDensity(eccen==eccentricities,:) = eccen2density(cMosaic, 'mm');
        end
        
        
        labels              = num2str(round(allDensity./10.^5,3));


        usedLabels          = contrastLevels; % For now use eccentricities as labels, but we could plot it against cone density

        
    case 'Defocus'
        eyemovement         = {'110'};
        defocusZ            = [0:0.5:2.0];
        coneType            = {''};
        colors              = copper(5);%{'k','r','g','b','c'};
        labels              = {'0 Diopters of Defocus','','', ...
                                '1.5 Diopters of Defocus','','', ...
                                '3.1 Diopters of Defocus','','',...
                                '4.6 Diopters of Defocus','','',...
                                '6.2 Diopters of Defocus','',''};
        eccentricities      = 4.5;
        contrastLevels      = [0:0.01:0.1];%, 0.1:0.1:1.0];
        usedLabels          = contrastLevels;
        
        
end



%% 0. Set other parameters
polarAngles = 0;
FFTflag     = true;

% Where to find data and save figures
dataPth     = fullfile(ogRootPath,'data','classification');
figurePth     = fullfile(ogRootPath,'figs');


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
for em = 1:length(eyemovement)
    for ct = 1:length(coneType)
        for eccen = eccentricities
            for df = 1:length(defocusZ)
                
                                
                %% 1. Load results
               
                fName   = sprintf('Classify_coneOutputs_contrast0.20_pa%d_eye%s_eccen%1.2f_defocus%1.2f_noise-random_phasescrambled.mat', ...
                    polarAngles,cell2mat(eyemovement(em)),eccen,defocusZ(df));
                 if currentFlag; fName = ['current_' fName]; end;
                
                accuracy = load(fullfile(dataPth, whatPlot, fName));
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

figure(3); clf; set(gcf,'Color','w'); hold all;

for ii = 1:length(fit.ctrpred)
    dataToFit = squeeze(fit.data{ii})';
    plot(usedLabels, fit.ctrpred{ii}*100, 'Color', colors(ii,:), 'LineWidth',2); 
    scatter(usedLabels, dataToFit, 80, colors(ii,:), 'filled');
    plot(10.^-2.1,dataToFit(1),'o','Color',colors(ii,:), 'MarkerSize', 8, 'MarkerFaceColor',colors(ii,:))
end

set(gca, 'XScale','log', 'XLim',[min(usedLabels) max(usedLabels)],'YLim', [min(dataToFit)-10 100], 'TickDir','out','TickLength',[.015 .015],'FontSize',17, 'LineWidth',2);
set(gca, 'XTick', [0.01,0.03,0.05, 0.1], 'XTickLabel',{'1','3','5','10'})

ylabel('Classifier Accuracy (% Correct)', 'FontSize',17)
xlabel('Stimulus Contrast (%)', 'FontSize',17);

legend(labels, 'Location','bestoutside'); legend boxoff


savefig(fullfile(figurePth,sprintf('WeibullFit_contrastVSperformance_all_pa%d_fft%d_%s_noPCA',polarAngles,FFTflag,whatPlot)))
hgexport(gcf,fullfile(figurePth,sprintf('WeibullFit_contrastVSperformance_all_pa%d_fft%d_%s_noPCA.eps',polarAngles,FFTflag,whatPlot)))

