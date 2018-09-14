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
eccInMMForSE    = [0.18, 0.27, 0.36, 0.45, 0.54, 0.72, 0.90, 1.08, 1.35, 1.62, 1.89, 2.16];
SEinMM.nasal    = [5.4, 2.8, 2.9, 2.7, 2.2, 1.8, 1.2, 1.5, 1.5, 1.4, 0.9, 0.6].*10^3;
SEinMM.superior = [ 2.3, 1.9, 2.0, 2.1, 1.9, 1.4, 1.2, 0.8, 0.8, 0.9, 0.8, 0.5].*10^3;
SEinMM.temporal = [7.5, 4.0, 2.4, 2.0, 1.7, 1.3, 1.5, 1.4, 0.9, 1.0, 0.5, 0.7].*10^3; 
SEinMM.inferior = [4.5, 1.4, 3.3, 2.8, 2.7, 2.1, 1.6, 1.7, 1.1, 0.9, 0.8, 0.7].*10^3;
eccInDegForSE   = eccInMMForSE * mm2deg;

% Create vectors with the length of eccentricities used for cone density, to
% index our limited numbers SE correctly
eccInDegToIndexSE            =  NaN(size(eccDeg));
SEinDeg.nasalMatchedIdx      = eccInDegToIndexSE;
SEinDeg.superiorMatchedIdx   = eccInDegToIndexSE;
SEinDeg.temporalMatchedIdx   = eccInDegToIndexSE;
SEinDeg.inferiorMatchedIdx   = eccInDegToIndexSE;

% Loop over eccentricities to place the limited SE points at the correct 
% eccentricity in degrees.
for c = 1:length(eccInDegForSE)
    
    % Round eccentricity
    roundedEcc = round(eccInDegForSE(c)*100)./100; 
    
    % Find index
    idx = find(eccDeg == roundedEcc);
    eccInDegToIndexSE(idx) = roundedEcc;
    
    % Place SE in correct index, and convert to cones/deg2
    SEinDeg.nasalMatchedIdx(idx)    = SEinMM.nasal(c)/(mm2deg.^2);
    SEinDeg.superiorMatchedIdx(idx) = SEinMM.superior(c)/(mm2deg.^2);
    SEinDeg.temporalMatchedIdx(idx) = SEinMM.temporal(c)/(mm2deg.^2);
    SEinDeg.inferiorMatchedIdx(idx) = SEinMM.inferior(c)/(mm2deg.^2);
end

% Get field names
fn = fieldnames(SEinDeg);

% Set up figure to plot cone densities
figure(1); clf; set(1,'Color','w', 'Position',[885, 327, 658, 598]); hold all;

% Loop over polar angles
for ii = 1:length(ang)
    
    % Loop over eccentricities
    for jj = 1:length(eccDeg)
        densityMM2(ii,jj) = coneDensityReadData('coneDensitySource', 'Song2011Young','eccentricity',eccDeg(jj),'angle',ang(ii),'eccentricityUnits', 'deg','whichEye','left');
    end
    
    % Convert density from mm squared to degrees squared
    densityDeg2(ii,:) = densityMM2(ii,:)/(mm2deg.^2);
    
    % Plot it!
    plot(eccDeg, densityDeg2(ii,:), 'Color', cmap(ii,:),  'LineWidth', 2);
    errorbar(eccDeg, densityDeg2(ii,:),SEinDeg.(fn{ii}), 'Color', cmap(ii,:),  'LineWidth', 2);
    
    % Make axes pretty
    xlabel('Eccentricity (deg)'); ylabel('Density (cones/deg^2)')
    set(gca,'XScale', 'log', 'YScale', 'log','YLim', 10.^[2.7, 4], 'XLim', 10.^[-0.3, 1]);
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

