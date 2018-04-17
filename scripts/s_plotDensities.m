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
    
    for jj = 1:length(ecc)
        [spacing(:,jj), aperture(:,jj), density(:,jj), params(:,jj), comment(:,jj)] = coneSizeReadData('eccentricity',ecc(jj),'angle',ang(ii),'eccentriticyUnits', 'm','whichEye','left');
    end
    
    subplot(411); hold all;
    plot([-flip(ecc) ecc]*m2deg, [flip(density) density], 'Color', cmap(ii,:), 'LineWidth', 2)
    xlabel('Eccentricity (deg)'); ylabel('Density (count/mm^2)')
    set(gca, 'FontSize', 16, 'TickDir', 'out', 'TickLength', [0.015 0.015], 'box', 'off')
    
    
    subplot(412); hold all;
    plot(ecc*m2deg, density, 'Color', cmap(ii,:),  'LineWidth', 2);
    set(gca,'XScale', 'log', 'YScale', 'log');
    xlabel('Eccentricity (deg)'); ylabel('Density (count/mm^2)')
    set(gca, 'FontSize', 16, 'TickDir', 'out', 'TickLength', [0.015 0.015], 'box', 'off')
   
    
    subplot(413); hold all;
    plot(ecc*m2deg, spacing, 'Color', cmap(ii,:),  'LineWidth', 2);
    set(gca,'XScale', 'log', 'YScale', 'log');
    xlabel('Eccentricity (deg)'); ylabel('Cone spacing (m)')
    set(gca, 'FontSize', 16, 'TickDir', 'out', 'TickLength', [0.015 0.015], 'box', 'off')
   
    
    subplot(414); hold all;
    plot(ecc*m2deg, aperture, 'Color', cmap(ii,:),  'LineWidth', 2);
    set(gca,'XScale', 'log', 'YScale', 'log');
    xlabel('Eccentricity (deg)'); ylabel('Cone aperture (m)')
    set(gca, 'FontSize', 16, 'TickDir', 'out', 'TickLength', [0.015 0.015], 'box', 'off')
    

end


subplot(411); axis tight; legend(angNames, 'Location', 'Best')
subplot(412); axis tight; legend(angNames, 'Location', 'Best')
subplot(413); axis tight; legend(angNames, 'Location', 'Best')
subplot(414); axis tight; legend(angNames, 'Location', 'Best')

% See where Banks data points fall
eccen = [0, 2.5, 5, 10, 20, 40];
for ii = 1:length(eccen)
    out(ii,:) = getBanks1991ConeSpacing(eccen(ii), 'unitOut', 'm');  
end

subplot(413); scatter(0.09, out(1), 80,'k'); scatter(eccen(2:end), out(2:end), 80, 'k');
subplot(414); scatter(0.09, out(1), 80,'k'); scatter(eccen(2:end), out(2:end), 80, 'k');

