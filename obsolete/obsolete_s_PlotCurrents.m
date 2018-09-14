%% obsolete_s_PlotCurrents    


figure(4); clf; set(4, 'Color','w', 'Position', [513, 668, 1462, 648], 'NumberTitle','off', 'Name', 'Fig 4 - Photocurrent');
figure(5); clf; set(5, 'Color','w', 'Position', [808, 98, 640, 647], 'NumberTitle','off', 'Name', 'Time series high contrast stim');
figure(6); clf; set(6, 'Color','w', 'Position', [808, 98, 640, 647], 'NumberTitle','off', 'Name', 'Time series low contrast stim');


% Select a trial
trialNum = 4;
stimulusNum = 4;

timeWindowC     = [20 21 22]; % select average time window for current (ms)
climsC = [-35 0]; % color bar / ylims for current (pA)
ylC    = [-70 0];  % y limit for time series (current, pA)


%% 4. CURRENT
    set(0, 'CurrentFigure', 4)
    
    % get average absorptions across time
    theseCurrents = squeeze(current(trialNum,:,:,timeWindowC,stimulusNum))./2;
    
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