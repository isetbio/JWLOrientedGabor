%% s_makeFigure4_CorrelationBehaviorAndConeDensity

% Plot the cone density (cones/deg2) for the 4 meridians (upper, lower,
% left and right) at 4.5 deg eccentricity (Song et al 2011 data).
% And plot these densities against the behavior of 3 observers from the
% paper Cameron, Tai, Carrasco (2002) showing a difference in thresholds
% for the four meridians.

% Unit Converters
deg2m   = 0.3 * 0.001;
m2deg   = 1/deg2m;

deg2mm  = 0.3; % Degrees of visual angle, 0.3 mm/deg.
mm2deg  = 1/deg2mm;

% Define eccentricity range
eccDeg = 4.5;
eccMM  = eccDeg*deg2mm;

% Define polar angles and their names
ang = [0 90 180 270];

% Preallocate matrices
densityMM2 = NaN(length(ang),length(eccMM));

cmap = [0 1 0; 0.5 0.5 0.5; 1 0 0; 0 0 1];

%% 1. Plot full cone density by eccentricity for different polar angles

% Get SE from Song et al. 2011 paper 
SEinMM2.nasal    = 1.5*10^3;
SEinMM2.superior = 0.8*10^3;
SEinMM2.temporal = 0.9*10^3; 
SEinMM2.inferior = 1.1*10^3;

% Create vectors with the length of eccentricities used for cone density, to
% index our limited numbers SE correctly
SEinDeg2.nasal    = SEinMM2.nasal/(mm2deg.^2);
SEinDeg2.superior = SEinMM2.superior/(mm2deg.^2);
SEinDeg2.temporal = SEinMM2.temporal/(mm2deg.^2);
SEinDeg2.inferior = SEinMM2.inferior/(mm2deg.^2);

% Loop over polar angles
for ii = 1:length(ang)
    
    % Get cone density in cones/mm2
    densityMM2(ii) = coneDensityReadData('coneDensitySource', 'Song2011Young','eccentricity',eccDeg,'angle',ang(ii),'eccentricityUnits', 'deg','whichEye','left');
    
    % Convert density from mm squared to degrees squared
    densityDeg2(ii) = densityMM2(ii,:)/(mm2deg.^2);
end


%% 2.Get thresholds at 4.5 degrees eccentricity

% behavior in percent correct for Nasal, Superior, Temporal, Inferior polar coords
behaviorJT = [2.3, 4.1, 2.5, 3.5];
behaviorJH = [4.2, 4.2, 4.0, 4.3];
behaviorLC = [3.8, 4.8, 3.8, 4.5];

meanBehavior4CPD = mean([behaviorJT;behaviorJH;behaviorLC],1);
seBehavior4CPD   = std([behaviorJT;behaviorJH;behaviorLC],[],1)/sqrt(3);

nasalConeDensityDeg2        = densityDeg2(1);
temporalConeDensityDeg2     = densityDeg2(3);
densityDeg2([1,3])  = mean([nasalConeDensityDeg2,temporalConeDensityDeg2]);

meanDensityDeg2     = densityDeg2;
seDensityDeg2(1)    = mean([SEinDeg2.nasal; SEinDeg2.temporal]);
seDensityDeg2(2)    = SEinDeg2.superior;
seDensityDeg2(3)    = mean([SEinDeg2.nasal; SEinDeg2.temporal]);
seDensityDeg2(4)    = SEinDeg2.inferior;

% Fit linear model
lm = fitlm(meanDensityDeg2,meanBehavior4CPD);

% Set up figure
figure(3); clf; set(3,'Color','w', 'Position',[276, 330, 658, 598]);

% Plot data and line
scatter(meanDensityDeg2, meanBehavior4CPD, [], cmap); hold on;
eb = [];
for ii = 1:4
    eb = errorbar(meanDensityDeg2(ii), meanBehavior4CPD(ii), seDensityDeg2(ii),'horizontal','LineStyle',':', 'Color', cmap(ii,:));
    eb.Bar.LineStyle = 'dashed';
    eb.Bar.LineWidth = 2;
    
    eb = errorbar(meanDensityDeg2(ii), meanBehavior4CPD(ii), seBehavior4CPD(ii),'vertical','LineStyle',':', 'Color', cmap(ii,:));
    eb.Bar.LineStyle = 'dashed';
    eb.Bar.LineWidth = 2;
end

y = (lm.Coefficients.Estimate(2)*meanDensityDeg2) + lm.Coefficients.Estimate(1);
plot(meanDensityDeg2,y, 'Color', [0 0 0], 'LineWidth',2); legend off;



% Make axes pretty
ylim([2.5 5]); xlim([1 1.6].*10^3); set(gca, 'FontSize', 16, 'TickDir', 'out', 'TickLength', [0.015 0.015], 'box', 'off');
xlabel('Cone density (cones/deg^2)'); ylabel('Threshold (% contrast)')
title(sprintf('Adjusted R-squared: %1.2f',lm.Rsquared.Adjusted))
axis square