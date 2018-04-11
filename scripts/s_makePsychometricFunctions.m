%% s_makePsychometricFunctions


% Script to compute psychometric functions based on the computational
% observer model


%% 0. Set general experiment parameters
expName = 'coneDensity';
expParams = loadExpParams(expName, false);

% Units on x axis
xUnits              = expParams.contrastLevels;

% Use cone current (flag = true) or cone absorptions (flag = false)
currentFlag = false;

polarAngles = 0;
FFTflag     = true;

% Where to find data and save figures
dataPth     = fullfile(ogRootPath,'data','classification');
classificationDir  = 'coneDensity'; % possibility to call a specific directory with the classification results, else, leave empty
figurePth   = fullfile(ogRootPath,'figs');

% Number of total trials in computational observer model (50 clockwise, 50 counterclockwise)
nTotal      = expParams.nTrials*4;


%% Get appropiate colors and labels

switch expName
    case 'default'
        colors              = [0 0 0];
        labels              = {'Fit','data'};
        
    case 'ConeTypes'
        colors              = {'k','r','g','b'};
        labels              = {'LMS cones','', ...
            'L cones','', ...
            'M cones','', ...
            'S cones',''};
        xUnits          = expParams.contrastLevels;
        
    case 'eyemov'
        colors              = copper(size(expParams.eyemovement,2));
        labels              = cell(1,length(colors));
        for emIdx = 1:length(colors)
            thisCondition = expParams.eyemovement(:,emIdx)';
            if all(thisCondition == [0 0 0])
                labels{:,emIdx} = 'No eyemovements';
            elseif  all(thisCondition == [1 0 0])
                labels{:,emIdx} = 'Tremor';
            elseif  all(thisCondition == [0 1 0])
                labels{:,emIdx} = 'Drift';
            elseif  all(thisCondition == [1 1 0])
                labels{:,emIdx} = 'Tremor+Drift';
            elseif  all(thisCondition == [0 0 1])
                labels{:,emIdx} = 'Microsaccades';
            end 
        end
        xUnits          = expParams.contrastLevels;
        
    case 'eyemovEnhanced'
        colors              = copper(size(expParams.eyemovement,2));
        labels              = cell(1,length(colors));
        for emIdx = 1:length(colors)
            thisCondition = expParams.eyemovement(:,emIdx)';
            if all(thisCondition == [0 0 0])
                labels{:,emIdx} = 'No eyemovements';
            elseif  all(thisCondition == [2 0 0])
                labels{:,emIdx} = '2xTremor';
            elseif  all(thisCondition == [0 2 0])
                labels{:,emIdx} = '2xDrift';
            elseif  all(thisCondition == [2 2 0])
                labels{:,emIdx} = '2xTremor+Drift';
            elseif  all(thisCondition == [0 0 2])
                labels{:,emIdx} = '2xMicrosaccades';
            end
        end
        
    case 'coneDensity'
      
        colors              = jet(length(expParams.eccentricities));
        
        % Change x labels to density
        whichEye          = 'left';
        cparams.cmFOV     =  1; % degrees
        
        % Convertion deg to m
        deg2m  = 0.3 * 0.001; % .3 deg per mm, .001 mm per meter
        
        % Predefine density vector
        allDensity = nan(length(expParams.eccentricities),1);
        labels = cell(length(expParams.eccentricities),1);
        
        for eccen = expParams.eccentricities
            % Specify retinal location where stimulus is presented
            cparams.eccentricity = eccen;             % Visual angle of stimulus center, in deg
            cparams.polarAngle   = deg2rad(0);   % Polar angle (radians): 0 is right, pi/2 is superior, pi is left, 3*pi/2 inferior
            
            % Compute x,y position in m of center of retinal patch from ecc and angle
            [x, y] = pol2cart(cparams.polarAngle, cparams.eccentricity);
            x = x * deg2m;  y = y * deg2m;
            cMosaic = coneMosaic('center', [x, y], 'whichEye', whichEye);
            
            % Set the field of view (degrees)
            cMosaic.setSizeToFOV(cparams.cmFOV);
            allDensity(eccen==expParams.eccentricities,:) = eccen2density(cMosaic, 'mm');
            
            labels{eccen==expParams.eccentricities} = sprintf('%1.3f x10^5 cells/mm2', round(allDensity(eccen==expParams.eccentricities)./10.^5,3));
        end
          
        xUnits              = expParams.contrastLevels; % For now use eccentricities as labels, but we could plot it against cone density
        
        
    case 'Defocus'
        colors              = copper(length(expParams.defocusLevels));
        labels              = {'0 Diopters of Defocus', ...
            '1.5 Diopters of Defocus', ...
            '3.1 Diopters of Defocus', ...
            '4.6 Diopters of Defocus', ...
            '6.2 Diopters of Defocus'};
        xUnits              = expParams.contrastLevels;
        
        
end



%% 0. Set other parameters


% Prepare fit variables
fit = [];

fit.ctrpred = cell(size(colors,1),1);
fit.ctrvar  = cell(size(colors,1),1);
fit.ctrr2   = cell(size(colors,1),1);
fit.data    = cell(size(colors,1),1);

fit.init = [0.5, 0.05];
fit.thresh = 0.75;

% Predefine matrix for predictions
nrEyemovTypes    = size(expParams.eyemovement,2);
nrEccen          = length(expParams.eccentricities);
nrSpatFreq       = length(expParams.spatFreq);
nrDefocusLevels  = length(expParams.defocusLevels);

fn = fieldnames(expParams);
if any(strcmp(fn(:),'cparams')); nrConeTypes = size(expParams.cparams.spatialDensity,2);
    error('s_makePsychometricFunctions does not allow multiple cone type conditions (yet)'); end

count = 1;
for em = 1:nrEyemovTypes
    for eccen = 1:nrEccen
        for df = 1:nrDefocusLevels
            
            
            %% 1. Load results
            
            fName   = sprintf('Classify_coneOutputs_contrast%1.3f_pa%d_eye%s_eccen%1.2f_defocus%1.2f_noise-random_sf%1.2f.mat', ...
                max(expParams.contrastLevels),polarAngles,sprintf('%i',expParams.eyemovement(:,em)'),expParams.eccentricities(eccen),expParams.defocusLevels(df),expParams.spatFreq);
            if currentFlag; fName = ['current_' fName]; end;
            
            accuracy = load(fullfile(dataPth, classificationDir, fName));
            accuracy.P = squeeze(accuracy.P);
            if size(accuracy.P,1)<size(accuracy.P,2)
                accuracy.P = accuracy.P';
            end
            
            %% 2. Fit Weibull
            % Make a Weibull function first with contrast levels and then search for
            % the best fit with the classifier data
            fit.ctrvar{count} = fminsearch(@(x) ogFitWeibull(x, xUnits, accuracy.P, nTotal), fit.init);
            
            % Then fit a Weibull function again, but now with the best fit parameters
            % from the previous step.
            fit.ctrpred{count} = ogWeibull(fit.ctrvar{count}, xUnits);
            
            % Not sure what this line is for..
            % fit.ctrr2 = corr(fit.ctrpred', dec.ctr).^2;
            
            %% 3. Find contrast threshold
            diff   = abs(fit.ctrpred{count} - fit.thresh);
            minval = find(diff == min(diff));
            fit.ctrthresh{count} = xUnits(minval(1));
            fit.data{count} = accuracy.P;
            
            count = count +1;
        end
        
    end
end


%% 4. Visualize

figure(3); clf; set(gcf,'Color','w'); hold all;

for ii = 1:length(fit.ctrpred)
    dataToFit = squeeze(fit.data{ii})';
    plot(xUnits, fit.ctrpred{ii}*100, 'Color', colors(ii,:), 'LineWidth',2);
    scatter(xUnits, dataToFit, 80, colors(ii,:), 'filled');
    plot(10.^-2.1,dataToFit(1),'o','Color',colors(ii,:), 'MarkerSize', 8, 'MarkerFaceColor',colors(ii,:))
end

set(gca, 'XScale','linear', 'XLim',[min(xUnits) max(xUnits)],'YLim', [min(dataToFit)-10 100], 'TickDir','out','TickLength',[.015 .015],'FontSize',17, 'LineWidth',2);
set(gca, 'XTick', [0, 0.01,0.03,0.05, 0.1], 'XTickLabel',{'0','1','3','5','10'})

ylabel('Classifier Accuracy (% Correct)', 'FontSize',17)
xlabel('Stimulus Contrast (%)', 'FontSize',17);

h = findobj(gca,'Type','line');
legend([h(2:2:end)],labels, 'Location','bestoutside'); legend boxoff


savefig(fullfile(figurePth,sprintf('WeibullFit_contrastVSperformance_fft%d_%s',FFTflag,expName)))
hgexport(gcf,fullfile(figurePth,sprintf('WeibullFit_contrastVSperformance_fft%d_%s.eps',FFTflag,expName)))

