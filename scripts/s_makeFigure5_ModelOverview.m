%% s_makeFigure5_ModelOverview.m

% Script to plot from stimulus (radiance, to photon absorptions).

% Requires s_ogRGC.m to be ran, and have save absorptions.

%% Set parameters

% Save figures?
saveFigs = false;

% Select a trial
trialNum = 4;
stimulusNum = 4;

timeWindowA = [1 2 3];  % select time points for absorptions (ms)
climsA = [0 220];  % color bar / ylims for absorptions (photon count)
ylA = [0 500];  % y limit for time series (absorptions, photon count)

% File names
fNames = {'absorptionsPlotted', 'stimulus'};
contrasts = [100, 10];

trialsToPlot = ceil(100*rand(5,1));

cmap = [1 0 0; 0 1 0; 0 0 1];
varySatValues = 1-linspace(0.1,1,5);
trialColors = varysat(cmap,varySatValues);
contrastColors = [0 0 0; 241 101 33]./255;

figPath = fullfile(ogRootPath, 'figs', 'overviewFigure5');

%% Plot the cone mosaic separate because of the different colormap
load(fullfile(figPath,'cMosaicPlotted'))

% Get size
m2deg = 10/3;
xPosInDeg = cMosaic.size .* m2deg;

figure(100); set(gcf, 'Color', 'w','Position', [1363, 916, 1112, 297], 'NumberTitle','off', 'Name', 'Fig 5 - Cone Mosaic');

% Plot mosaic pattern
% subplot(131); 
title('Cone Mosaic')
imagesc(cMosaic.pattern)
set(gca,'CLim', [2 4]); axis image; axis off
colormap(cmap)

% subplot(132); hold all; title('Normalized quanta absorbance')
% for ii = 1:3
%     plot(cMosaic.wave, cMosaic.pigment.absorbance(:,ii), 'Color', cmap(ii,:))
% end
% box off; ylabel('Normalized quanta absorbance', 'FontSize',12);
% xlabel('Wavelength (nm)', 'FontSize',12); set(gca, 'TickDir', 'out', 'FontSize', 12);
% 
% subplot(133); hold on; title('Normalized macular transmittance')
% plot(cMosaic.wave, cMosaic.macular.transmittance, 'Color', 'k')
% box off; ylabel('Normalized macular transmittance', 'FontSize',12);
% xlabel('Wavelength (nm)', 'FontSize',12); set(gca, 'TickDir', 'out', 'FontSize', 12);
% 
if saveFigs
    hgexport(100,fullfile(ogRootPath, 'figs','Fig2_coneMosaic.eps'))
end


% Set up figures
figure(1); clf; set(1, 'Color','w', 'Position', [513, 668, 1462, 648], 'NumberTitle','off', 'Name', 'Fig 5 - Radiance'); 
figure(2); clf; set(2, 'Color','w', 'Position', [513, 668, 1462, 648], 'NumberTitle','off', 'Name', 'Fig 5 - Irradiance'); 
figure(3); clf; set(3, 'Color','w', 'Position', [513, 668, 1462, 648], 'NumberTitle','off', 'Name', 'Fig 5 - Absorptions');


% Loop over stimulus contrast
for c = 1:length(contrasts)
    
    % Load data
    for ii = 1:length(fNames)
        load(fullfile(figPath, [fNames{ii} sprintf('_%d.mat',contrasts(c))]));
    end
        
    %% 1. RADIANCE
    set(0, 'CurrentFigure', 1)
    
    radiance = scenes{2}.data.photons;
    midpoint = ceil(size(radiance,1)/2);
    
    % Plot scene radiance
    subplot(1,3,c); cla; hold all;
    imagesc(sum(radiance,3));
    hold on; plot([0 size(radiance,1)], [midpoint midpoint], 'k:');
    set(gca,'CLim', [0, 2.3].*10^17);
    colormap gray; colorbar; axis image; axis off
    title('Scene radiance summed over nm','Fontsize',12)
    
    % Plot 1D scene radiance
    subplot(1,3,3); hold on;
    plot(sum(radiance(midpoint,:, :),3), 'Color', contrastColors(c,:));
    title('1D radiance','Fontsize',12); axis image; axis square; box off;
    xlabel('X pos (deg)');
    ylabel('radiance (q/s/sr/nm/m^2)')    
    set(gca, 'TickDir', 'out', 'XTick', [1, midpoint, size(radiance,2)], ...
    'XTickLabel',{sprintf('%1.1f',-1*scenes{2}.wAngular/2), '0',sprintf('%1.1f',-1*scenes{2}.wAngular/2)}, ...
    'TickDir', 'out', 'FontSize', 12);

    xlim([1 size(scenes{2}.data.photons,1)]);
    ylim([0, 2.3].*10^17)
    
        
    if saveFigs
        if c==2
            print(fullfile(figPath, 'Fig2_radiance_2d'),'-depsc')
        end
    end
   
    
    %% 2. IRRADIANCE
    set(0, 'CurrentFigure', 2)
    
    % Plot the stimulus after optics
    irradiance = OG(4).oiModulated.data.photons;
    midpoint   = ceil(size(irradiance,1)/2);
    gamma      = 1;
    wList = (400:10:700); 
    
    if c == 1 % Create color bar;
       
        XYZ = ieXYZFromPhotons(irradiance, wList);
        XYZ_normalized = XYZ/max(XYZ(:));    
        rgbIrradianceData = xyz2srgb(XYZ_normalized);
        
        rgbColorMapMin = min(reshape(rgbIrradianceData, [size(rgbIrradianceData,1)*size(rgbIrradianceData,2), 3]));
        rgbColorMapMax = max(reshape(rgbIrradianceData, [size(rgbIrradianceData,1)*size(rgbIrradianceData,2), 3]));
        cInt = [linspace(rgbColorMapMin(1), rgbColorMapMax(1), 64); ...
            linspace(rgbColorMapMin(2), rgbColorMapMax(2), 64); ...
            linspace(rgbColorMapMin(3), rgbColorMapMax(3), 64)]';
        
        irradianceSumAllWV = sum(irradiance,3);
        ylI = [min(irradianceSumAllWV(:)), max(irradianceSumAllWV(:))];
    else
        rgbIrradianceData = xyz2srgb(ieXYZFromPhotons(irradiance, wList)./max(XYZ(:)));
    end
    
%     irradiance2D = reshape(irradiance, [size(irradiance,1)*size(irradiance,2),size(irradiance,3)]);
    
    subplot(1,3,c); hold all;
    imagesc(rgbIrradianceData);
    hold on; plot([0 size(rgbIrradianceData,1)], [midpoint midpoint], 'k:');
    colormap(cInt); axis image; axis off; cb = colorbar;
    set(gca,'CLim', [0 1]); cb.Ticks = linspace(0, 1, 4) ; %Create 8 ticks from zero to 1
    cb.TickLabels = num2cell([ylI(1), ylI(2)*(1/3),  ylI(2)*(2/3), ylI(2)]);
    title('Retinal irradiance','Fontsize',12)
    
    subplot(1,3,3); hold on;
    plot(sum(irradiance(midpoint,:,:),3), 'Color', contrastColors(c,:));
    title('1D irradiance','Fontsize',12); axis square; box off;
    xlabel('X pos (deg)');  
    xlim([1 size(OG(4).oiModulated.data.photons,1)])
    ylabel('Irradiance (q/s/m^2/nm)'); 
    ylim(ylI)
    set(gca, 'TickDir', 'out',... 
    'XTick', [1, midpoint, size(irradiance,2)], ...
    'XTickLabel',{sprintf('%1.1f',-1*OG(4).oiModulated.wAngular/2), '0',sprintf('%1.1f',-1*OG(4).oiModulated.wAngular/2)},...
    'FontSize', 12);
   
    
    if saveFigs
        if c==2
            print(fullfile(figPath, 'Fig2_irradiance_2d'),'-depsc')
        end
    end
   
%     ylim([min(sum(irradiance2D,2)) max(sum(irradiance2D,2))]);
    % set(gca, 'YTick', [0:0.333:1], 'YTickLabel', cb.TickLabels);
    
    
    
    %% 3. ABSORPTIONS

    set(0, 'CurrentFigure', 3)
    
    % get average absorptions across time
    theseAbsorptions = squeeze(absorptions(trialNum,:,:,timeWindowA,stimulusNum))./2;
    
    % Integration time was 2 ms, thus recalculate to absorption per 1 ms
    % theseAbsorptionsMS = theseAbsorptions; %./(cMosaic.integrationTime*1000);
    
    midpoint = ceil(size(theseAbsorptions,1)/2);
    
    if c == 1; subplotIdx = [1:3]; else; subplotIdx = [5:7]; end
    
    for ii = subplotIdx
        subplot(2,4,ii); hold on;
        imagesc(theseAbsorptions(:,:,mod(ii,4)));
        colormap gray; axis image; axis off;
        set(gca,'CLim',climsA); colorbar;
        title(sprintf('Absorptions (photons) at t=%d',mod(ii-1,3)+1),'Fontsize',12)
        if ii == subplotIdx(end); plot([0 size(theseAbsorptions,1)], [midpoint midpoint], 'k:', 'LineWidth',2); end
        colormap gray; axis image; axis off;
        set(gca,'CLim',climsA);
    end
    
    subplot(2,4,4); hold all;
    plot(theseAbsorptions(midpoint,:,3), 'k', 'Color', contrastColors(c,:));
    title('1D Absorptions','Fontsize',12); axis image; axis square; box off;
    xlabel('X position (deg)'); ylabel('Average absorption rate (photons/ms)')
    set(gca, 'XLim', [1 size(theseAbsorptions,1)], ...
        'TickDir', 'out', 'FontSize', 12,...
        'YLim',climsA, ...
        'XTick', [1, midpoint, size(theseAbsorptions,2)], ...
        'XTickLabel',{sprintf('%1.0f',(-1*xPosInDeg(1)/2)*1000), '0',sprintf('%1.0f',(xPosInDeg(1)/2)*1000)});
    
    if saveFigs
        if c ==2
            print(fullfile(figPath, 'Fig2_absorptions_2d'),'-depsc')
        end
    end
    

    
    %% save classifier weights
    %
    % data = absorptions(:,:,:,1:28,:);
    %
    % % Get the trials and samples (should be the data for all data sets though
    % nStimuli = size(data,5);
    % nTrials  = size(data,1) * nStimuli/2;
    % tSamples = size(data,4);
    % nrows    = size(data,2);
    % ncols    = size(data,3);
    %
    % % absorptions is trials x rows x cols x time points x stimuli
    %
    %
    % %   permute to trials x stimuli x rows x cols x time points
    % data = permute(data, [1 5 2:4]);
    %
    % %   reshape to (trials x stimuli) x rows x cols x time points
    % data = reshape(data, [], nrows, ncols, tSamples);
    %
    % % permute to rows x cols x (trials x stimuli) x time points
    % data  = permute(data, [2 3 1 4]);
    %
    % % Compute fourier transform the cone array outputs
    % data  = abs(fft2(data));
    %
    % % reshape to all trials x [rows x colums x time] for classification
    % data = permute(data, [3 1 2 4]);
    % data = reshape(data, nTrials*2, []);
    %
    % % permute the trial order within each of the two classes
    % idx = [randperm(nTrials) randperm(nTrials)+nTrials];
    %
    % data = data(idx, :);
    %
    % label = [ones(nTrials, 1); -ones(nTrials, 1)];
    %
    % % Fit the SVM model.
    % cvmdl = fitcsvm(data, label, 'Standardize', true, 'KernelFunction', 'linear', 'kFold', 10);
    %
    % betas = reshape(cvmdl.Trained{1}.Beta, [nrows, ncols, tSamples]);
    % save(fullfile(ogRootPath, 'figs',sprintf('fftBetas_%d.mat',contrasts(c))),'betas')
    
    
end % for contrasts


%% Plot the fft and the classification accuracy + weibull function


% Load experiment parameters
expParams = loadExpParams('default', false);
[xUnits, contrastColors, labels, M] = loadWeibullPlottingParams('default');

figure(8); clf; set(gcf, 'Color', 'w', 'Position', [974, 790, 735, 241], 'NumberTitle','off', 'Name', 'Fig 5 - Classifier Weights');

for c = [1,2]    
    % Load beta's
    load(fullfile(figPath, sprintf('fftBetas_%d.mat',contrasts(c))));
    midpoint = ceil(size(betas,1)/2);
    
    % Apply FFTshift, take average over time samples
    mn_betas = fftshift(squeeze(mean(betas,3)));
    
    % Visualize mean betas
    subplot(1,2,c);
    imagesc(fftshift(betas(:,:,10))); colormap gray; colorbar;
    title(sprintf('FFT of classifier weights %d%',contrasts(c)), 'FontSize',17);
    set(gca,'CLim', 1*10^-3.5*[-1 1], ...
        'FontSize', 17, ...
        'XTick', [midpoint*0.5, midpoint, midpoint*1.5], ...
        'XTickLabel',{'-10', '0','-10'});
        xlabel('cycles/deg'); box off;
end

if saveFigs
    print(fullfile(figPath, 'Fig2_classifierweights'),'-deps2')
end


%% Load accuracy for 0-10% contrasts
accuracy = load(fullfile(figPath, 'classification_accuracy_10.mat'));

% Set inital slope, threshold for first stage fitting
fit.init   = [2, 0.02]; % slope, threshold at ~80%
fit.thresh = 0.75;

% Make a Weibull function first with contrast levels and then search for
% the best fit with the classifier data
fit.ctrvar = fminsearch(@(x) ogFitWeibull(x, expParams.contrastLevels, accuracy.P, expParams.nTrials), fit.init);

% Then fit a Weibull function again, but now with the best fit parameters
% from the previous step.
fit.ctrpred = ogWeibull(fit.ctrvar, xUnits);

%Find contrast threshold
fit.ctrthresh = fit.ctrvar(2);
fit.data = accuracy.P;

% Visualize psychometric curves
figure(9); clf; set(gcf,'Color','w', 'Position',  [1282, 771, 621, 492], 'NumberTitle','off', 'Name', 'Fig 5 - Psychometric function'); hold all;

dataToFit = fit.data;
plot(xUnits(2:end), fit.ctrpred(2:end)*100, 'k', 'LineWidth',2);
scatter(expParams.contrastLevels(2:end), dataToFit(2:end), 80, 'k', 'filled');
plot(3e-3,dataToFit(1),'ko', 'MarkerSize', 8, 'MarkerFaceColor','k')
set(gca, 'XScale','log','XLim',[3e-3, max(expParams.contrastLevels)],'YLim', [40 100], 'TickDir','out','TickLength',[.015 .015],'FontSize',17, 'LineWidth',2);
set(gca, 'XTick', [3e-3, expParams.contrastLevels(2:2:end)], 'XTickLabel',sprintfc('%1.1f',[0 expParams.contrastLevels(2:2:end)]*100))
title('Model performance using cone absorptions', 'FontSize',17)
ylabel('Classifier Accuracy (% Correct)', 'FontSize',17)
xlabel('Stimulus Contrast (%)', 'FontSize',17);


if saveFigs
    print(fullfile(figPath, 'Fig2_weibull'),'-deps2')
end

