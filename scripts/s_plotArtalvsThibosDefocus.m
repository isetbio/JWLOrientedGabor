pupilMM = 3;
waveToPlot = 550;
[wvf, oi] = wvfLoadJaekenArtal2012Data('jIndex', 0:14, 'whichEye','left', 'eccentricity',4);

wvf = wvfSet(wvf,'calc pupil size',pupilMM);
wvf = wvfSet(wvf,'calc wave',waveToPlot);

uData = wvfPlot(wvf,'OTF','mm',waveToPlot);
set(gca,'xlim',[-250 250] ,'ylim',[-250 250], 'zlim', [0 1.2],'FontSize',16,'LineWidth',2);

fH = figure(1); set(fH,'Position',[0.0039,0.0257,0.5941,0.9076]); clf;
subplot(3,3,1); surf(uData.fx, uData.fy, uData.otf); view(2); colormap parula; title('Average of Artal data')
set(gca,'xlim',[-250 250] ,'ylim',[-250 250], 'zlim', [0 1.2],'FontSize',16,'LineWidth',2);
ylabel('freq/mm')
for ii = 0:1:7
    
    clear wvfP;
    zCoefs = wvfLoadThibosVirtualEyes(pupilMM);
    wave = 550;% wave = wave(:);
    
    wvfP = wvfCreate('wave',wave,'zcoeffs',zCoefs,'name',sprintf('human-%d',pupilMM));
    wvfP = wvfSet(wvfP,'calc pupil size',pupilMM);
    wvfP = wvfSet(wvfP,'calc wave',waveToPlot);
    wvfP = wvfSet(wvfP,'zcoeffs',ii,{'defocus'});
    wvfP = wvfComputePSF(wvfP);
    uData = wvfPlot(wvfP,'OTF','mm',waveToPlot);
    %
    %      view(2)
    %
    % %     otfSupport = wvfGet(wvf, 'otfsupport', 'mm', waveToPlot);
    %     otf = wvfGet(wvfP, 'otf', waveToPlot);
    figure(fH); subplot(3,3,ii+2); surf(uData.fx, uData.fy, uData.otf); view(2); colormap parula; 
    title(sprintf('Thibos with defocus %1.1f Diopters', zernickeDefocus2diopter(ii, pupilMM)));
    set(gca,'xlim',[-250 250] ,'ylim',[-250 250], 'zlim', [0 1.2],'FontSize',16,'LineWidth',2);
    
    if ismember(ii+2,[4,7])
        ylabel('freq (lines/mm)'); end
    if ismember(ii+2, [7,8,9])
        xlabel('freq (lines/mm)'); end
    
    
end
