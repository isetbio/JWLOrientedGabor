%% s_makeFigure2_ConeDensity

% Create figure 2 of manuscript, showing cone Density as a function of
% eccentricity and polar angle (Song et al. 2011), plus zoom, plus
% correlation with behavior of indivual observer JT in from Cameron, Tai,
% Carrasco (2002).

%% 0. Define parameters
% Unit Converters
deg2m   = 0.3 * 0.001;
m2deg   = 1/deg2m;

deg2mm  = 0.3; % Degrees of visual angle, 0.3 mm/deg.
mm2deg  = 1/deg2mm;

% Define eccentricity range
eccDeg = (0.5:.05:7);
eccMM  = eccDeg*deg2mm;

% Define polar angles and their names
ang = [0 90 180 270];
angNames = {'Nasal (HM)','', 'Superior (LVM)','', 'Temporal (HM)', '', 'Inferior (UVM)',''};

% Preallocate matrices
densityMM2 = NaN(length(ang),length(eccMM));
params = cell(length(ang),length(eccMM));
comment = cell(length(ang),length(eccMM));

cmap = [0 1 0; 0.5 0.5 0.5; 1 0 0; 0 0 1];

%% 1. Plot full cone density by eccentricity for different polar angles

% Get SE from Song et al. 2011 paper 
SEinMM.ecc      = [0.18, 0.27, 0.36, 0.45, 0.54, 0.72, 0.90, 1.08, 1.35, 1.62, 1.89, 2.16];
SEinMM.nasal    = [5.4, 2.8, 2.9, 2.7, 2.2, 1.8, 1.2, 1.5, 1.5, 1.4, 0.9, 0.6].*10^3;
SEinMM.superior = [ 2.3, 1.9, 2.0, 2.1, 1.9, 1.4, 1.2, 0.8, 0.8, 0.9, 0.8, 0.5].*10^3;
SEinMM.temporal = [7.5, 4.0, 2.4, 2.0, 1.7, 1.3, 1.5, 1.4, 0.9, 1.0, 0.5, 0.7].*10^3; 
SEinMM.inferior = [4.5, 1.4, 3.3, 2.8, 2.7, 2.1, 1.6, 1.7, 1.1, 0.9, 0.8, 0.7].*10^3;
SEinDeg.ecc      = SEinMM.ecc * mm2deg;

% Place SE eccentricies in vector that matches eccentricies (deg) vector
SEinDeg.eccMatched        =  NaN(size(eccDeg));
SEinDeg.nasalMatched      = SEinDeg.eccMatched;
SEinDeg.superiorMatched   = SEinDeg.eccMatched;
SEinDeg.temporalMatched   = SEinDeg.eccMatched;
SEinDeg.inferiorMatched   = SEinDeg.eccMatched;

% Loop over eccentricities to place the limited SE points at the correct 
% eccentricity in degrees.
for c = 1:length(SEinDeg.ecc)
    
    roundedEcc = round(SEinDeg.ecc(c)*100)./100; 
    
    idx = find(eccDeg == roundedEcc);
    SEinDeg.eccMatched(idx) = roundedEcc;
    
    SEinDeg.nasalMatched(idx)    = (sqrt(SEinMM.nasal(c))*mm2deg).^2;
    SEinDeg.superiorMatched(idx) = (sqrt(SEinMM.superior(c))*mm2deg).^2;
    SEinDeg.temporalMatched(idx) = (sqrt(SEinMM.temporal(c))*mm2deg).^2;
    SEinDeg.inferiorMatched(idx) = (sqrt(SEinMM.inferior(c))*mm2deg).^2;
end

% Get field names
fn = fieldnames(SEinDeg);

% Set up figure to plot cone densities
figure(1); clf; set(1,'Color','w', 'Position',[885, 327, 658, 598]); hold all;

% Loop over polar angles
for ii = 1:length(ang)
    
    % Loop over eccentricities
    for jj = 1:length(eccMM)
        densityMM2(ii,jj) = coneDensityReadData('coneDensitySource', 'Song2011Young','eccentricity',eccDeg(jj),'angle',ang(ii),'eccentricityUnits', 'deg','whichEye','left');
    end
    
    % Convert density from mm squared to degrees squared
    densityDeg(ii,:) = (sqrt(densityMM2(ii,:))*mm2deg).^2;
    
    % Plot it!
    plot(eccDeg, densityDeg(ii,:), 'Color', cmap(ii,:),  'LineWidth', 2);
    errorbar(eccDeg, densityDeg(ii,:),SEinDeg.(fn{2+ii}), 'Color', cmap(ii,:),  'LineWidth', 2);
    
    % Make axes pretty
    xlabel('Eccentricity (deg)'); ylabel('Density (cones/deg^2)')
    set(gca,'XScale', 'log', 'YScale', 'log','YLim', 10.^[4.7, 6], 'XLim', 10.^[-0.3, 1]);
    set(gca, 'FontSize', 16, 'TickDir', 'out', 'TickLength', [0.015 0.015], 'box', 'off');
    grid on; set(gca, 'GridLineStyle','-', 'MinorGridLineStyle', ':');
end

% Add legend
legend(angNames, 'Location', 'Best'); legend boxoff;
axis square


%% 2. Make inset of cone density plot for 3-7 degrees eccentricities

% % Find eccentricities
% idx3Deg = find(eccDeg==3);
% idx7Deg = find(eccDeg==7);
% 
% % Set up figure 
% figure(2); clf; set(2,'Color','w', 'Position',[1796, 456, 320, 342]); hold all;
% 
% % Loop over polar angles
% for ii = 1:4
%     plot(eccDeg, densityDeg(ii,:), 'Color', cmap(ii,:),  'LineWidth', 2);
%     errorbar(eccDeg, densityDeg(ii,:),SEinDeg.(fn{2+ii}), 'Color', cmap(ii,:),  'LineWidth', 2);
% end
% 
% % Make axes pretty
% set(gca,'XScale', 'log', 'YScale', 'log','YLim', 10.^[2.5, 3.5], 'XLim', 10.^[0.4771, 0.8451]);
% xlabel('Eccentricity (deg)'); ylabel('Density (cones/deg^2)')
% set(gca, 'FontSize', 16, 'TickDir', 'out', 'TickLength', [0.015 0.015], 'box', 'off');
% grid on; set(gca, 'GridLineStyle','-', 'MinorGridLineStyle', ':')

%% 3. Make correlation plot cone density at 4.5 degrees eccentricity & behavior

idx = find(eccDeg==4.5);

% behavior = [85, 69, 87, 65]; % percent correct for Nasal, Superior, Temporal, Inferior polar coords
behaviorJT = [2.3, 4.1, 2.5, 3.5];
behaviorJH = [4.2, 4.2, 4.0, 4.3];
behaviorLC = [3.8, 4.8, 3.8, 4.5];

meanBehavior4CPD = mean([behaviorJT;behaviorJH;behaviorLC],1);
seBehavior4CPD  = std([behaviorJT;behaviorJH;behaviorLC],[],1)/sqrt(4);


nasalConeDensity = densityDeg(1,idx);
temporalConeDensity = densityDeg(3,idx);
densityDeg([1,3],idx) = mean([nasalConeDensity,temporalConeDensity]);

coneDensMean  = densityDeg(:,idx);
coneDensSE(1) = mean([SEinDeg.nasalMatched(idx); SEinDeg.temporalMatched(idx)]);
coneDensSE(2) = SEinDeg.superiorMatched(idx);
coneDensSE(3) = mean([SEinDeg.nasalMatched(idx); SEinDeg.temporalMatched(idx)]);
coneDensSE(4) = SEinDeg.inferiorMatched(idx);

% Fit linear model
lm = fitlm(coneDensMean,meanBehavior4CPD);

% Set up figure
figure(3); clf; set(3,'Color','w', 'Position',[276, 330, 658, 598]);

% Plot data and line
scatter(coneDensMean, meanBehavior4CPD, [], cmap); hold on;
eb = [];
for ii = 1:4
    eb = errorbar(coneDensMean(ii), meanBehavior4CPD(ii), coneDensSE(ii),'horizontal','LineStyle',':', 'Color', cmap(ii,:));
    eb.Bar.LineStyle = 'dashed';
    eb.Bar.LineWidth = 2;
    
    eb = errorbar(coneDensMean(ii), meanBehavior4CPD(ii), seBehavior4CPD(ii),'vertical','LineStyle',':', 'Color', cmap(ii,:));
    eb.Bar.LineStyle = 'dashed';
    eb.Bar.LineWidth = 2;
end

y = (lm.Coefficients.Estimate(2)*coneDensMean) + lm.Coefficients.Estimate(1);
plot(coneDensMean,y, 'Color', [0 0 0], 'LineWidth',2); legend off;



% Make axes pretty
ylim([2.5 5]); xlim([1.2 2.1].*10^5); set(gca, 'FontSize', 16, 'TickDir', 'out', 'TickLength', [0.015 0.015], 'box', 'off');
xlabel('Cone density (cones/deg^2)'); ylabel('Threshold (% contrast)')
title(sprintf('Adjusted R-squared: %1.2f',lm.Rsquared.Adjusted))
axis square