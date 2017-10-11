
thresh = [0.01,0.02,0.02,0.03,0.05];

% thresh = [0.01,0.02,0.03,0.03,0.02,0.04,0.03,0.03,0.05];

thresh =  [0.0100,0.0200,0.0300,0.0500];


%%

M =    allDensity/11.111
lm = fitlm(log10(M),thresh);


figure(2); clf;
plot(lm, 'LineWidth', 3, 'MarkerSize',10, 'Marker','o','Color',[0 0 0]); box off;
set(gca, 'TickDir', 'out','TickLength',[0.015 0.015], 'LineWidth',1,'Fontsize',25,'XScale','linear')
xlabel('Cone Density (cones/deg^2)','FontSize',25); ylabel('Contrast sensitivity threshold','FontSize',25)
set(gca, 'XTick',[2, 3, 4],'XTickLabel',[100 1000 10000], 'XLim', [1.99 5],'YLim', [-0.01 0.08]),
legend off
% title('Contrast threshold versus Cone density')

% hgexport(gcf, '/Volumes/server/Projects/PerformancefieldsIsetBio/figs/contrastThreshVsConeDensity.eps')