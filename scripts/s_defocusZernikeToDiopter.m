%% s_defocusZernikeToDiopter

defocusZ = 0:0.5:2; % zernike coeffs for Defocus

thresh = [0.0600,0.0900,0.0700,0.1000,0.1000];  % stim contrast thresholds at ~80% for each defocus condition

M = 4*pi*sqrt(3) * defocusZ / (pi* 1.5^2); % convert to diopters

lm = fitlm(M,thresh);

%%
figure(1); clf;
plot(lm, 'LineWidth', 3, 'MarkerSize',10, 'Marker','o'); box off;
set(gca, 'TickDir', 'out','TickLength',[0.015 0.015], 'LineWidth',1,'Fontsize',12, 'XLim',[-1 7])
xlabel('Defocus (Diopters)'); ylabel('Contrast sensitivity threshold')
legend off
title('Contrast threshold versus defocus')

hgexport(gcf, '/Volumes/server/Projects/PerformancefieldsIsetBio/figs/contrastThreshVsDiopter.eps')

return