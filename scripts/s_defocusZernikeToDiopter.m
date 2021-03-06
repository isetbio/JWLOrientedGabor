%% s_defocusZernikeToDiopter

defocusZ = 0:0.5:2; % zernike coeffs for Defocus

thresh = [0.0200,0.0200,0.0300,0.0300,0.0400];  % stim contrast thresholds at ~80% for each defocus condition

% Where does this formula come from? Polans et al? 
%   see here for starters: https://www.researchgate.net/post/How_can_I_convert_Zernike_modes_in_micrometer_to_Diopter
pupilRadiusMM = 1.5; % mm
M = 4*pi*sqrt(3) * defocusZ / (pi* pupilRadiusMM^2); % convert to diopters

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




