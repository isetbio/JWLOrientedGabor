%% s_plotFigure4_plotStimToCurrent.m

% Script to plot from stimulus (radiance, to photo current).

% Requires s_ogRGC.m to be ran, save absorptions, current,

%% Set parameters

% Save figures?
saveFigs = false;

% Wavelength to plot
wv = 550;

% Select a trial
trialNum = 4;
stimulusNum = 4;

wA = [1 2 3];  % select time points for absorptions (ms)
wC = [20 21 22]; % select average time window for current (ms)

climsA = [0 220];  % color bar / ylims for absorptions (photon count)
climsC = [-35 0]; % color bar / ylims for current (pA)


ylA = [0 500];  % y limit for time series (absorptions, photon count)
ylC = [-70 0];  % y limit for time series (current, pA)

% File names
fNames = {'absorptionsPlotted', 'currentPlotted','stimulus'};
contrasts = [100, 10];

trialsToPlot = ceil(100*rand(5,1));

cmap = [1 0 0; 0 1 0; 0 0 1];
varySatValues = 1-linspace(0.1,1,5);
trialColors = varysat(cmap,varySatValues);
contrastColors = [0 0 0; 241 101 33]./255;

figPath = fullfile(ogRootPath, 'figs', 'overviewFigure4');

%% Plot the cone mosaic separate because of the different colormap
load(fullfile(figPath,'cMosaicPlotted'))

% Get size
m2deg = 10/3;
xPosInDeg = cMosaic.size .* m2deg;

figure(100); set(gcf, 'Color', 'w','Position', [1363, 916, 1112, 297]);

% Plot mosaic pattern
subplot(131); title('Cone Mosaic')
imagesc(cMosaic.pattern)
set(gca,'CLim', [2 4]); axis image; axis off
colormap(cmap)

subplot(132); hold all; title('Normalized quanta absorbance')
for ii = 1:3
    plot(cMosaic.wave, cMosaic.pigment.absorbance(:,ii), 'Color', cmap(ii,:))
end
box off; ylabel('Normalized quanta absorbance', 'FontSize',12);
xlabel('Wavelength (nm)', 'FontSize',12); set(gca, 'TickDir', 'out', 'FontSize', 12);

subplot(133); hold on; title('Normalized macular transmittance')
plot(cMosaic.wave, cMosaic.macular.transmittance, 'Color', 'k')
box off; ylabel('Normalized macular transmittance', 'FontSize',12);
xlabel('Wavelength (nm)', 'FontSize',12); set(gca, 'TickDir', 'out', 'FontSize', 12);

if saveFigs
    hgexport(100,fullfile(ogRootPath, 'figs','Fig2_coneMosaic.eps'))
end



% Set up figures
figure(1); clf; set(1, 'Color','w', 'Position', [513, 668, 1462, 648], 'NumberTitle','off', 'Name', 'Fig 4 - Radiance'); 
figure(2); clf; set(2, 'Color','w', 'Position', [513, 668, 1462, 648], 'NumberTitle','off', 'Name', 'Fig 4 - Irradiance'); 
figure(3); clf; set(3, 'Color','w', 'Position', [513, 668, 1462, 648], 'NumberTitle','off', 'Name', 'Fig 4 - Absorptions');
figure(4); clf; set(4, 'Color','w', 'Position', [513, 668, 1462, 648], 'NumberTitle','off', 'Name', 'Fig 4 - Photocurrent');
figure(5); clf; set(5, 'Color','w', 'Position', [808, 98, 640, 647], 'NumberTitle','off', 'Name', 'Time series high contrast stim');
figure(6); clf; set(6, 'Color','w', 'Position', [808, 98, 640, 647], 'NumberTitle','off', 'Name', 'Time series low contrast stim');

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
    
    if c == 1; % Create color bar;
       
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
    theseAbsorptions = squeeze(absorptions(trialNum,:,:,wA,stimulusNum))./2;
    
    % Integration time was 2 ms, thus recalculate to absorption per 1 ms
    % theseAbsorptionsMS = theseAbsorptions; %./(cMosaic.integrationTime*1000);
    
    midpoint = ceil(size(theseAbsorptions,1)/2);
    
    if c == 1; subplotIdx = [1:3]; else; subplotIdx = [5:7]; end
    
    for ii = subplotIdx
        subplot(2,4,ii); hold on;
        imagesc(theseAbsorptions(:,:,mod(ii,4)));
        colormap gray; axis image; axis off;
        set(gca,'CLim',climsA); colorbar;
        title(sprintf('Absorptions (photons) at t=%d',ii),'Fontsize',12)
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
        'XTickLabel',{sprintf('%1.3f',(-1*xPosInDeg(1)/2)), '0',sprintf('%1.3f',(xPosInDeg(1)/2))});
    
    if saveFigs
        if c ==2
            print(fullfile(figPath, 'Fig2_absorptions_2d'),'-depsc')
        end
    end
    
    %% 4. CURRENT
    set(0, 'CurrentFigure', 4)
    
    % get average absorptions across time
    theseCurrents = squeeze(current(trialNum,:,:,wC,stimulusNum))./2;
    
    % Integration time was 2 ms, thus recalculate to absorption per 1 ms
    % theseAbsorptionsMS = theseAbsorptions; %./(cMosaic.integrationTime*1000);
    
    midpoint = ceil(size(theseCurrents,1)/2);
    
    if c ==1; subplotIdx = [1:3]; else; subplotIdx = [5:7]; end
    
    for ii = subplotIdx
        subplot(2,4,ii); hold on;
        imagesc(theseCurrents(:,:,mod(ii,4)));
        colormap gray; axis image; axis off
        set(gca,'CLim',climsC); colorbar;
        title(sprintf('Photocurrent (pA/ms) at t=%d',ii),'Fontsize',12)
        if ii == subplotIdx(end); plot([0 size(theseCurrents,1)], [midpoint midpoint],  'k:', 'LineWidth',2); end
        colormap gray; axis image; axis off;
        set(gca,'CLim',climsC);
    end
    
    subplot(2,4,4); hold all;
    plot(theseCurrents(midpoint,:,3), 'k', 'Color', contrastColors(c,:));
    title('1D Photocurrent','Fontsize',12); axis square; box off;
    xlabel('X pos (deg)'); ylabel('Average photocurrent rate (pA/ms)')
    set(gca, 'XLim', [1 size(theseCurrents,1)], ...
        'TickDir', 'out', 'FontSize', 12, ...
        'YLim',climsC, ...
        'XTick', [1, midpoint, size(theseCurrents,2)], ...
        'XTickLabel',{sprintf('%1.3f',(-1*xPosInDeg(1)/2)), '0',sprintf('%1.3f',xPosInDeg(1)/2)});
    
    if saveFigs
        if c ==2
            print(fullfile(figPath, 'Fig2_current_2d'),'-depsc')
        end
    end
    %% 5. PLOT INDIVIDUAL CONE TIMESERIES
    
    % Visualize the mean response of the best responding L-, M- and S-cone
    % Compute absorptions and current per cone class as cones x tPoints x trials
    absorptions1 = permute(squeeze(absorptions(:,:,:,:,stimulusNum)), [2 3 4 1]);
    current1     = permute(squeeze(current(:,:,:,:,stimulusNum)), [2 3 4 1]);
    [nRows, mCols, tBins, nTrials] = size(absorptions1);
    absorptions1 = reshape(absorptions1, [nRows*mCols tBins nTrials]);
    current1     = reshape(current1, [nRows*mCols tBins nTrials]);
    
    % Find indices of cones for each of the 3 cone types
    lcones = cMosaic.pattern == 2;
    mcones = cMosaic.pattern == 3;
    scones = cMosaic.pattern == 4;
    
    % Find indices of the best responding L, M ans S-cone
    l_absorptions = absorptions1(lcones,:,:)./2;
    m_absorptions = absorptions1(mcones,:,:)./2;
    s_absorptions = absorptions1(scones,:,:)./2;
    l_curr = current1(lcones,:,:)./2;
    m_curr = current1(mcones,:,:)./2;
    s_curr = current1(scones,:,:)./2;
    
    if c ==1
        [~,maxIndex] = max(l_absorptions(:));
        [maxLconeIndex, ~, ~] = ind2sub(size(l_absorptions), maxIndex);
        [~,maxIndex] = max(m_absorptions(:));
        [maxMconeIndex, ~, ~] = ind2sub(size(m_absorptions), maxIndex);
        [~,maxIndex] = max(s_absorptions(:));
        [maxSconeIndex, ~, ~] = ind2sub(size(s_absorptions), maxIndex);
    end
        % Compute mean over reps of best responding L-cone
        meanLabsorption = squeeze(mean(l_absorptions(maxLconeIndex,:, :),3));
        % stdLabsorption = squeeze(std(l_absorptions(maxLconeIndex,:, :),[], 3));

        % Compute mean over reps of best responding M-cone
        meanMabsorption = squeeze(mean(m_absorptions(maxMconeIndex,:, :),3));
        % stdMabsorption = squeeze(std(m_absorptions(maxMconeIndex,:, :),[], 3));

        % Compute mean over reps of best responding S-cone
        meanSabsorption = squeeze(mean(s_absorptions(maxSconeIndex,:, :),3));
        % stdSabsorption = squeeze(std(s_absorptions(maxSconeIndex,:, :),[], 3));
    

    %% plot the timeseries of an absorption rate
    set(0, 'CurrentFigure', 4+c)
    if c ==1  
        ylA = [50, 250; 25, 175; 0 30]; 
        ylC = [-20 0; -20 0; -40 -20]; 
    else
         ylA = [80, 150; 50, 120; 0 30]; 
         ylC = [-20 0; -20 0; -40 -20];
    end
    
    % One L Cone time series
    subplot(231); hold all;
    % patch([0 0.054 0.054 0], [ylA(1) ylA(1) ylA(2) ylA(2)], [0.7 0.7 0.7], 'EdgeColor','none')
    % x2 = [cMosaic.timeAxis, fliplr(cMosaic.timeAxis)];
    % y2 = [meanLabsorption+stdLabsorption, fliplr(meanLabsorption-stdLabsorption)];
    % h = fill(x2,y2 ,'r', 'edgecolor','none');
    % set(h,'facealpha',.5);
    
    for ii = 1:length(trialsToPlot)
            plot(cMosaic.timeAxis, squeeze(l_absorptions(maxLconeIndex,:, trialsToPlot(ii))),'Color',trialColors(:,ii,1),'LineWidth', 0.5);
    end
    plot(cMosaic.timeAxis, meanLabsorption, 'Color', 'k', 'LineWidth',2);
    plot([.054 .054], ylA(1,:), 'k--');
    xlim([0 0.216]); ylim(ylA(1,:));
    set(gca,'TickDir','out','Fontsize',12, 'LineWidth',1)
    title('L cone')
    xlabel('Time (s)'); ylabel('Absorption (photons/ms)')
    
    % One M Cone time series
    subplot(232); hold all;
    % patch([0 0.054 0.054 0], [ylA(1) ylA(1) ylA(2) ylA(2)], [0.7 0.7 0.7], 'EdgeColor','none')
    % x2 = [cMosaic.timeAxis, fliplr(cMosaic.timeAxis)];
    % y2 = [meanMabsorption+stdMabsorption, fliplr(meanMabsorption-stdMabsorption)];
    % h = fill(x2,y2 ,'g', 'edgecolor','none');
    % set(h,'facealpha',.5);
    for ii = 1:length(trialsToPlot)
        plot(cMosaic.timeAxis, squeeze(m_absorptions(maxMconeIndex,:, trialsToPlot(ii))),'Color',trialColors(:,ii,2),'LineWidth', 0.5);
    end
    plot(cMosaic.timeAxis, meanMabsorption, 'Color', 'k', 'LineWidth',2);
    plot([.054 .054], ylA(2,:), 'k--');
    xlim([0 0.216]); ylim(ylA(2,:));
    set(gca,'TickDir','out','Fontsize',12, 'LineWidth',1)
    title('M cone')
    xlabel('Time (s)'); ylabel('Absorption (photons/ms)')
    
    % One S Cone time series
    subplot(233); hold all;
    % patch([0 0.054 0.054 0], [ylA(1) ylA(1) ylA(2) ylA(2)], [0.7 0.7 0.7], 'EdgeColor','none')
    % x2 = [cMosaic.timeAxis, fliplr(cMosaic.timeAxis)];
    % y2 = [meanSabsorption+stdSabsorption, fliplr(meanSabsorption-stdSabsorption)];
    % h = fill(x2,y2 ,'b', 'edgecolor','none');
    % set(h,'facealpha',.5);
    for ii = 1:length(trialsToPlot)
        plot(cMosaic.timeAxis, squeeze(s_absorptions(maxSconeIndex,:, trialsToPlot(ii))),'Color',trialColors(:,ii,3),'LineWidth', 0.5);
    end
    plot(cMosaic.timeAxis, meanSabsorption, 'Color', 'k', 'LineWidth',2);
    plot([.054 .054], ylA(3,:), 'k--');
    xlim([0 0.216]); ylim(ylA(3,:));
    set(gca,'TickDir','out','Fontsize',12, 'LineWidth',1)
    title('S cone')
    xlabel('Time (s)'); ylabel('Absorption (photons/ms)')
    
    % Now for currents
    
    % Compute mean over reps of best responding L-cone
    meanLcurrent = squeeze(mean(l_curr(maxLconeIndex,:, :),3));
    stdLcurrent = squeeze(std(l_curr(maxLconeIndex,:, :),[], 3));
    
    % Compute mean over reps of best responding M-cone
    meanMcurrent = squeeze(mean(m_curr(maxMconeIndex,:, :),3));
    stdMcurrent = squeeze(std(m_curr(maxMconeIndex,:, :),[], 3));
    
    % Compute mean over reps of best responding S-cone
    meanScurrent = squeeze(mean(s_curr(maxSconeIndex,:, :),3));
    stdScurrent = squeeze(std(s_curr(maxSconeIndex,:, :),[], 3));
    
    % One L Cone time series
    subplot(234); cla; hold all;
    for ii = 1:length(trialsToPlot)
        plot(cMosaic.timeAxis, l_curr(maxLconeIndex,:,trialsToPlot(ii)), 'Color',trialColors(:,ii,1),'LineWidth', 0.5); end
    % patch(2*[wC(1) wC(2) wC(2) wC(1)]./1000, [ylC(1) ylC(1) ylC(2) ylC(2)], [0.7 0.7 0.7], 'EdgeColor','none')
    plot(cMosaic.timeAxis, meanLcurrent, 'Color','k', 'LineWidth',2);
    % x2 = [cMosaic.timeAxis, fliplr(cMosaic.timeAxis)];
    % y2 = [meanLcurrent+stdLcurrent, fliplr(meanLcurrent-stdLcurrent)];
    % h = fill(x2,y2 ,'r', 'edgecolor','none');
    % set(h,'facealpha',.5);
    plot([.054 .054], ylC(1,:), 'k--');
    xlim([0 0.216]);  ylim(ylC(1,:));
    set(gca,'TickDir','out','Fontsize',12, 'LineWidth',1)
    title('L cone')
    xlabel('Time (s)'); ylabel('Photocurrent (pA)')
    
    % One M Cone time series
    subplot(235); cla; hold all;
    for ii = 1:length(trialsToPlot)
        plot(cMosaic.timeAxis, m_curr(maxMconeIndex,:,trialsToPlot(ii)), 'Color',trialColors(:,ii,2),'LineWidth', 0.5); end
    % patch(2*[wC(1) wC(2) wC(2) wC(1)]./1000 , [ylC(1) ylC(1) ylC(2) ylC(2)], [0.7 0.7 0.7], 'EdgeColor','none')
    plot(cMosaic.timeAxis, meanMcurrent, 'Color','k', 'LineWidth',2);
    % x2 = [cMosaic.timeAxis, fliplr(cMosaic.timeAxis)];
    % y2 = [meanMcurrent+stdMcurrent, fliplr(meanMcurrent-stdMcurrent)];
    % h = fill(x2,y2 ,'g', 'edgecolor','none');
    % set(h,'facealpha',.5);
    plot([.054 .054], ylC(2,:), 'k--');
    xlim([0 0.216]); ylim(ylC(2,:));
    set(gca,'TickDir','out','Fontsize',12, 'LineWidth',1)
    title('M cone')
    xlabel('Time (s)'); ylabel('Photocurrent (pA)')
    
    % One S Cone time series
    subplot(236); cla; hold all;
    for ii = 1:length(trialsToPlot)
        plot(cMosaic.timeAxis, s_curr(maxSconeIndex,:,trialsToPlot(ii)), 'Color',trialColors(:,ii,3),'LineWidth', 0.5);
    end
    % patch(2*[wC(1) wC(2) wC(2) wC(1)]./1000, [ylC(1) ylC(1) ylC(2) ylC(2)], [0.7 0.7 0.7], 'EdgeColor','none')
    plot(cMosaic.timeAxis, meanScurrent, 'Color','k', 'LineWidth',2);
    % x2 = [cMosaic.timeAxis, fliplr(cMosaic.timeAxis)];
    % y2 = [meanScurrent+stdScurrent, fliplr(meanScurrent-stdScurrent)];
    % h = fill(x2,y2 ,'b', 'edgecolor','none');
    % set(h,'facealpha',.5);
    plot([.054 .054], ylC(3,:), 'k--');
    xlim([0 0.216]); ylim(ylC(3,:));
    set(gca,'TickDir','out','Fontsize',12, 'LineWidth',1)
    title('S cone')
    xlabel('Time (s)'); ylabel('Photocurrent (pA)')
    
    %
    if saveFigs
        print(fullfile(ogRootPath, 'figs', sprintf('Fig2_timeseries_%d',c)),'-depsc')
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

figure(8); clf; set(gcf, 'Color', 'w', 'Position', [974, 790, 735, 241]);

for c = [1,2]
    
    % Load beta's
    load(fullfile(figPath, sprintf('fftBetas_%d.mat',contrasts(c))));
    
    % Apply FFTshift, take average over time samples
    mn_betas = fftshift(squeeze(mean(betas,3)));
    
    % Visualize mean betas
    subplot(1,2,c);
    imagesc(fftshift(betas(:,:,10))); colormap gray; colorbar;
    title(sprintf('FFT of classifier weights %d%',contrasts(c)), 'FontSize',17);
    set(gca,'CLim', 1*10^-3.5*[-1 1], 'FontSize', 17); box off;
    
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
figure(9); clf; set(gcf,'Color','w', 'Position',  [1282, 771, 621, 492]); hold all;

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












