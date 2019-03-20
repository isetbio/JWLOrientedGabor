%% s_plotDensitiesTutorial

% Tutorial script to plot cone Density as a function of eccentricity and polar
% angle. See also new ISETBIO function coneDensityReadData.m


% Note that there are multiple sources in the literature for Cone Density:
%   Curcio (default); Song etc.
% Compare to FOV, p.46, Figure 3.1


% Unit Converters
deg2m   = 0.3 * 0.001;
m2deg   = 1/deg2m;

% Define eccentricity range
ecc = (0.1:.05:40) * deg2m;

% Define polar angles and their names
ang = [0 90 180 270];
angNames = {'Nasal (HM)','Superior (LVM)','Temporal (HM)', 'Inferior (UVM)'};

% Preallocate matrices
spacing = NaN(length(ang),length(ecc));
aperture = spacing; 
density = spacing;
params = cell(length(ang),length(ecc));
comment = cell(length(ang),length(ecc));

% Plot cone densities
figure('Color','w', 'Position',[99, 45, 1230, 1300]); cmap = lines(4);
for ii = 1:length(ang)
    
    for jj = 1:length(ecc)
        [spacing(ii,jj), aperture(ii,jj), density(ii,jj), params{ii,jj}, comment{ii,jj}] = coneSizeReadData('eccentricity',ecc(jj),'angle',ang(ii),'eccentriticyUnits', 'm','whichEye','left');
    end
    
    subplot(411); hold all;
    plot([-flip(ecc) ecc]*m2deg, [flip(density(ii,:)) density(ii,:)], 'Color', cmap(ii,:), 'LineWidth', 2)
    xlabel('Eccentricity (deg)'); ylabel('Density (count/mm^2)')
    set(gca, 'FontSize', 16, 'TickDir', 'out', 'TickLength', [0.015 0.015], 'box', 'off')
    
    
    subplot(412); hold all;
    plot(ecc*m2deg, density(ii,:), 'Color', cmap(ii,:),  'LineWidth', 2);
    set(gca,'XScale', 'log', 'YScale', 'log');
    xlabel('Eccentricity (deg)'); ylabel('Density (count/mm^2)')
    set(gca, 'FontSize', 16, 'TickDir', 'out', 'TickLength', [0.015 0.015], 'box', 'off')
   
    
    subplot(413); hold all;
    plot(ecc*m2deg, spacing(ii,:), 'Color', cmap(ii,:),  'LineWidth', 2);
    set(gca,'XScale', 'log', 'YScale', 'log');
    xlabel('Eccentricity (deg)'); ylabel('Cone spacing (m)')
    set(gca, 'FontSize', 16, 'TickDir', 'out', 'TickLength', [0.015 0.015], 'box', 'off')
   
    
    subplot(414); hold all;
    plot(ecc*m2deg, aperture(ii,:), 'Color', cmap(ii,:),  'LineWidth', 2);
    set(gca,'XScale', 'log', 'YScale', 'log');
    xlabel('Eccentricity (deg)'); ylabel('Cone aperture (m)')
    set(gca, 'FontSize', 16, 'TickDir', 'out', 'TickLength', [0.015 0.015], 'box', 'off')
    

end


subplot(411); axis tight; legend(angNames, 'Location', 'Best')
subplot(412); axis tight; legend(angNames, 'Location', 'Best')
subplot(413); axis tight; legend(angNames, 'Location', 'Best')
subplot(414); axis tight; legend(angNames, 'Location', 'Best')

%% See where Banks data points fall and what we use in our Gabor experiment
eccen = [0, 2.5, 5, 10, 20, 40];
for ii = 1:length(eccen)
    out(ii,:) = getBanks1991ConeSpacing(eccen(ii), 'spacingUnit', 'm');  
end

toPlotBanks = [[0.09, eccen(2:end)]; out'];
subplot(413); scatter(toPlotBanks(1,:),toPlotBanks(2,:), 80, 'k', 'fill', 'k');
subplot(414); scatter(toPlotBanks(1,:),toPlotBanks(2,:), 80, 'k', 'fill', 'k');

% Get used parameters
p       = loadExpParams('ConeDensity', false);

for ii = 1:length(p.eccentricities)
    cm = coneMosaic('center',[p.eccentricities(ii)*deg2m,0],'whichEye','left');
    cm.setSizeToFOV(1);
    
    %% OPTION A We read out these parameters from the coneMosaic
    
    % S = Spacing between center cones in meters, defined as:
    %   1./conesPerM, where conesPerM = conesPerMM*1e3, where conePerMM =
    %   sqrt(density);
    % For a reference: cm.pigment.width - cm.pigment.pdWidth = cm.pigment.gapWidth
    % And:             cm.pigment.width = spacing and cm.pigment.pdWidth = aperture
    s1(:,ii) = cm.pigment.width;
    s_tmp(:,ii) = cm.pigment.gapWidth;
    
    % A = Inner cone aperture in meters, defined as:
    %  0.7*spacing
    a1(:,ii) = cm.pigment.pdWidth;
    
    % D = Density of cones within a patch, defined by:
    % density = eccen2density(cMosaic, unit)
    d1(:,ii) = (1./(prod(cm.patternSampleSize)))*10.^-6;
    
    %% OPTION B: Use coneSizeReadData, this function is used when creating a cMosaic, so should give the same answers  
    [s2(:,ii), a2(:,ii), d2(:,ii)] = coneSizeReadData('eccentricity',p.eccentricities(ii)*deg2m,'angle',p.polarAngle,'eccentricityUnits', 'm','whichEye','left');
end

toPlotPF_spacing = [[0.09, p.eccentricities(2:end)]; s1];
toPlotPF_aperture = [[0.09, p.eccentricities(2:end)]; a1];
toPlotPF_density = [[0.09, p.eccentricities(2:end)]; d1];

subplot(412); scatter(toPlotPF_density(1,:), toPlotPF_density(2,:), 80, 'r', 'fill', 'r'); legend({angNames{:},'PF Density Data'});
subplot(413); scatter(toPlotPF_spacing(1,:), toPlotPF_spacing(2,:), 80, 'r', 'fill', 'r'); legend({angNames{:}, 'Banks Spacing Data ', 'PF Spacing Data'});
subplot(414); scatter(toPlotPF_aperture(1,:), toPlotPF_aperture(2,:), 80, 'r', 'fill', 'r'); legend({angNames{:}, 'Banks Spacing Data ', 'PF Aperture Data'});


