%% s_ConeDensityClassifierAccuracy

% 0. Set other parameters
polarAngles = 0;
FFTflag     = true;

% Where to find data
dataPth     = fullfile(ogRootPath,'data','classification','ConeDensity');

% Number of total trials in computational observer model (50 clockwise, 50 counterclockwise)
nTotal      = 100;

eyemovement         = {'110'};
defocusZ            = 0;
coneType            = {''};
eccentricities      = [0,4.5,16,60];
colors              = copper(length(eccentricities));

contrastLevels      = [0:0.01:0.1];

%% Change x labels to density
whichEye          = 'left';
cparams.cmFOV     =  1; % use 1 degree of FOV for labels

% Convertion deg to m
deg2m  = .3 * 0.001; % 3 deg per mm, .001 mm per meter

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

usedLabels          = contrastLevels;

%% Prepare fit variables
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
                fName   = sprintf('Classify_coneOutputs_contrast0.10_pa%d_eye%s_eccen%1.2f_defocus%1.2f_noise-random_phasescrambled_fft%d.mat', ...
                    polarAngles,cell2mat(eyemovement(em)),eccen,defocusZ(df),FFTflag);
                
                
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




% thresh = [0.01,0.02,0.02,0.03,0.05];

% thresh = [0.01,0.02,0.03,0.03,0.02,0.04,0.03,0.03,0.05];

% thresh =  [0.0100,0.0200,0.0300,0.0500];
thresh = cell2mat(fit.ctrthresh);

%%

M  = allDensity/11.111; % Convert mm2 to deg2
lm = fitlm(log10(M),thresh);


figure(2); clf;
plot(lm, 'LineWidth', 3, 'MarkerSize',10, 'Marker','o','Color',[0 0 0]); box off;
set(gca, 'TickDir', 'out','TickLength',[0.015 0.015], 'LineWidth',1,'Fontsize',25,'XScale','linear')
xlabel('Cone Density (cones/deg^2)','FontSize',25); ylabel('Contrast sensitivity threshold','FontSize',25)
set(gca, 'XTick',[2, 3, 4],'XTickLabel',[100 1000 10000], 'XLim', [1.99 5],'YLim', [-0.01 0.08]),
legend off
% title('Contrast threshold versus Cone density')

% hgexport(gcf, '/Volumes/server/Projects/PerformancefieldsIsetBio/figs/contrastThreshVsConeDensity.eps')