%% s_coneDensityVSClassifierAccuracy

eccentricities = 0:60;

load('/Users/winawerlab/matlab/git/toolboxes/JWLOrientedGabor/figs/Classify_coneOutputs_contrast0.08_pa0_eye110_eccen60.00_defocus0.00_noise-random_phasescrambled_fft1.mat')

P = squeeze(P);

lm = fitlm(eccentricities,P);


%% Visualize

figure(1); clf;
plot(lm, 'LineWidth', 4, 'MarkerFace',[0 0 0],'MarkerSize',10, 'Marker','o', 'Color',[0 0 0], 'LineWidth',4); box off;
set(gca, 'TickDir', 'out','TickLength',[0.015 0.015], 'LineWidth',2,'Fontsize',25, 'XLim',[0 max(eccentricities)])
xlabel('Cone Density x10^11 (deg)','FontSize',25); ylabel('Classifier Accuracy', 'FontSize',25)
legend off; 
title('Eccentricity versus Cone Density versus Classifier Accuracy for 0.08 stimulus contrast','FontSize',25)


%% Change x labels to density
whichEye          = 'left';
cparams.cmFOV     =  2; % degrees

% Convertion deg to m
deg2m  = 1/3 * 0.001; % 3 deg per mm, .001 mm per meter

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
    allDensity(eccen==eccentricities,:) = eccen2density(cMosaic, 'm');
end


ylabel('Accuracy (% Correct))','Fontsize',25); 
set(gca,'XTick', 0:10:60)

set(gca,'XTickLabel',round((allDensity(1:10:61)/10.^5),3)); 
%%

% hgexport(1,'/Volumes/server/Projects/PerformancefieldsIsetBio/figs/eccenVSconeDensityVSPerformance')

return