%% s_makePsychometricFunctions


% Script to compute psychometric functions based on the computational
% observer model

%% 0. Set general experiment parameters
expName = 'defocus';
expParams = loadExpParams(expName, false);

% Units on x axis
xUnits              = expParams.contrastLevels;

% Use cone current (flag = true) or cone absorptions (flag = false)
currentFlag = false;

polarAngles = 0;
FFTflag     = true;

% Where to find data and save figures
dataPth     = fullfile(ogRootPath,'data','classification','HPC',expName);
figurePth   = fullfile(ogRootPath,'figs','HPC');

% Number of total trials in computational observer model (50 clockwise, 50 counterclockwise)
nTotal      = expParams.nTrials*4;


%% Get appropiate colors and labels

switch lower(expName)
    case 'default'
        colors              = [0 0 0];
        labels              = {'Fit','data'};
        
    case 'conetypes'
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
        
    case 'eyemovenhanced'
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
        
    case {'conedensity','coneDensitynoeyemov'}
      
        colors              = jet(length(expParams.eccentricities));
        
        % Change x labels to density
        whichEye          = 'left';
        cparams.cmFOV     =  1; % degrees
        
        % Convertion deg to m
        deg2m  = 0.3 * 0.001; % .3 deg per mm, .001 mm per meter
        
        % Predefine density vector
        allDensity = nan(length(expParams.eccentricities),1);
        labels = cell(length(expParams.eccentricities),1);
        
        for ec = expParams.eccentricities
            % Specify retinal location where stimulus is presented
            cparams.eccentricity = ec;                     % Visual angle of stimulus center, in deg
            cparams.polarAngle   = expParams.polarAngle;   % Polar angle (radians): 0 is right, pi/2 is superior, pi is left, 3*pi/2 inferior
            
            % Compute x,y position in m of center of retinal patch from ecc and angle
            [x, y] = pol2cart(cparams.polarAngle, cparams.eccentricity);
            x = x * deg2m;  y = y * deg2m;
            cMosaic = coneMosaic('center', [x, y], 'whichEye', whichEye);
            
            % Set the field of view (degrees)
            cMosaic.setSizeToFOV(cparams.cmFOV);
            allDensity(ec==expParams.eccentricities,:) = eccen2density(cMosaic, 'mm');
            
            labels{ec==expParams.eccentricities} = sprintf('%1.3f x10^5 cells/mm2', allDensity(ec==expParams.eccentricities)./10.^5);
        end
          
        xUnits              = expParams.contrastLevels; % For now use eccentricities as labels, but we could plot it against cone density
        
        
    case 'defocus'
        
        for df = expParams.defocusLevels
            
            % compute defocus
            pupilRadiusMM = 1.5; % mm
            M(df==expParams.defocusLevels) = 4*pi*sqrt(3) * df / (pi* pupilRadiusMM^2); % convert to diopters
            labels{df==expParams.defocusLevels} = sprintf('%2.2f Diopters of Defocus',M(df==expParams.defocusLevels));
        end
        
        colors              = copper(length(expParams.defocusLevels));
        xUnits              = expParams.contrastLevels;
 
end



%% 0. Set other parameters


% Prepare fit variables
fit = [];

fit.ctrpred = cell(size(colors,1),1);
fit.ctrvar  = cell(size(colors,1),1);
fit.ctrr2   = cell(size(colors,1),1);
fit.data    = cell(size(colors,1),1);

fit.init = [2, 0.02]; % slope, threshold at ~80%
fit.thresh = 0.75;

% Predefine matrix for predictions
nrEyemovTypes    = size(expParams.eyemovement,2);
nrEccen          = length(expParams.eccentricities);
nrSpatFreq       = length(expParams.spatFreq);
nrDefocusLevels  = length(expParams.defocusLevels);

fn = fieldnames(expParams);
if any(strcmp(fn(:),'cparams')); nrConeTypes = size(expParams.cparams.spatialDensity,2);
    error('s_makePsychometricFunctions does not allow multiple cone type conditions (yet)'); end

% conditions = [nrEyemovTypes, nrEccen, nrSpatFreq, nrDefocusLevels];
% condOfInterest = find(conditions>1);

count = 1;
for em = 1:nrEyemovTypes
    for eccen = 1:nrEccen
        for df = 1:nrDefocusLevels
            
            
            %% 1. Load results
            
            fName   = sprintf('Classify_coneOutputs_contrast%1.3f_pa%d_eye%s_eccen%1.2f_defocus%1.2f_noise-random_sf%1.2f.mat', ...
                max(expParams.contrastLevels),polarAngles,sprintf('%i',expParams.eyemovement(:,em)'),expParams.eccentricities(eccen),expParams.defocusLevels(df),expParams.spatFreq);
            if currentFlag; fName = ['current_' fName]; end;
            
            accuracy = load(fullfile(dataPth, fName));
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

figure(3); clf; set(gcf,'Color','w', 'Position',  [1000, 850, 986, 488]); hold all;

for ii = 1:length(fit.ctrpred)
    dataToFit = fit.data{ii};
    plot(xUnits(2:end), fit.ctrpred{ii}(2:end)*100, 'Color', colors(ii,:), 'LineWidth',2);
    scatter(xUnits(2:end), dataToFit(2:end), 80, colors(ii,:), 'filled');
    plot(3e-3,dataToFit(1),'o','Color',colors(ii,:), 'MarkerSize', 8, 'MarkerFaceColor',colors(ii,:))
end

set(gca, 'XScale','log','XLim',[3e-3, max(xUnits)],'YLim', [min(dataToFit)-10 100], 'TickDir','out','TickLength',[.015 .015],'FontSize',17, 'LineWidth',2);
set(gca, 'XTick', [3e-3, xUnits(2:2:end)], 'XTickLabel',sprintfc('%1.1f',[0 xUnits(2:2:end)]*100))

ylabel('Classifier Accuracy (% Correct)', 'FontSize',17)
xlabel('Stimulus Contrast (%)', 'FontSize',17);

h = findobj(gca,'Type','line');
legend([h(end:-2:2)],labels, 'Location','bestoutside'); legend boxoff

if ~exist(figurePth,'dir'); mkdir(figurePth); end
savefig(fullfile(figurePth,sprintf('WeibullFit_contrastVSperformance_fft%d_%s',FFTflag,expName)))
hgexport(gcf,fullfile(figurePth,sprintf('WeibullFit_contrastVSperformance_fft%d_%s.eps',FFTflag,expName)))


% Plot density thresholds
if all(ismember('coneDensity',expName))
    
    thresh = cell2mat(fit.ctrthresh);
    M  = allDensity/11.111; % Convert mm2 to deg2
    lm = fitlm(log10(M),thresh);


    figure(2); clf; set(gcf, 'Color', 'w', 'Position', [1318, 696, 836, 649])
    plot(lm, 'LineWidth', 3, 'MarkerSize',10, 'Marker','o','Color',[0 0 0]); box off;
    set(gca, 'TickDir', 'out','TickLength',[0.015 0.015], 'LineWidth',1,'Fontsize',25,'XScale','linear')
    xlabel('Cone Density (cones/deg^2)','FontSize',25); ylabel('Contrast sensitivity threshold','FontSize',25)
    set(gca, 'XTick',[2, 3, 4],'XTickLabel',[100 1000 10000], 'XLim', [1.99 5],'YLim', [-0.01 0.08]),
    legend off; title('Contrast threshold versus Cone density')
    
elseif strcmp(expName,'defocus')
    
     thresh = cell2mat(fit.ctrthresh);    
     lm = fitlm(M,thresh);
     
     figure(2); clf; set(gcf, 'Color', 'w', 'Position', [1318, 696, 836, 649])
        plot(lm, 'LineWidth', 3, 'MarkerSize',10, 'Marker','o','Color',[0 0 0]); box off;
        set(gca, 'TickDir', 'out','TickLength',[0.015 0.015], 'LineWidth',1,'Fontsize',25,'XScale','linear')
        xlabel('Defocus (diopters)','FontSize',25); ylabel('Contrast sensitivity threshold','FontSize',25)
%         set(gca, 'XTick',[2, 3, 4],'XTickLabel',[100 1000 10000], 'XLim', [1.99 5],'YLim', [-0.01 0.08]),
        legend off; title('Contrast threshold versus Defocus')

    
end
