%% s_plotDensities

% Cone Density, see new function coneDensityReadData.m


% Note that there are multiple sources in the literature for Cone Density:
%   Curcio (default); Song etc.
% Compare to FOV, p.46, Figure 3.1


% Unit Converters
m2deg = 1000*3;
deg2m = 1/m2deg;

% Define eccentricity range
ecc = (0.1:.05:20) * deg2m;

% Define polar angles
ang = [0 90 180 270];

angNames = {'Nasal (HM)','Superior (LVM)','Temporal (HM)', 'Inferior (UVM)'};

% Plot cone densities
figure('Color','w'); cmap = lines(4);
for ii = 1:length(ang)
    
    coneDensity = coneDensityReadData('eccentricity',ecc,'angle',ones(size(ecc))*ang(ii),'whichEye','left');

    subplot(211); hold all;
    plot([-flip(ecc) ecc]*m2deg, [flip(coneDensity) coneDensity], 'Color', cmap(ii,:), 'LineWidth', 2)
    xlabel('Eccentricity (deg)'); ylabel('Density (count/mm^2)')
    set(gca, 'FontSize', 16, 'TickDir', 'out', 'TickLength', [0.015 0.015], 'box', 'off')
    
    subplot(212); hold all;
    plot(ecc*m2deg, coneDensity, 'Color', cmap(ii,:),  'LineWidth', 2);
    set(gca,'XScale', 'log', 'YScale', 'log');
    xlabel('Eccentricity (deg)'); ylabel('Density (count/mm^2)')
    set(gca, 'FontSize', 16, 'TickDir', 'out', 'TickLength', [0.015 0.015], 'box', 'off')

end
subplot(211), axis tight; legend(angNames, 'Location', 'Best')
subplot(212), axis tight; legend(angNames, 'Location', 'Best')

% Plot RGC densities