% s_plotMTF_defocusLevels

pupilMM = 3;
defocusMicrons = [0, 1, 2];
% supportUnit = 'deg';
zCoefs = wvfLoadThibosVirtualEyes(pupilMM);
% wave = 400:10:700; wave = wave(:);
calc_wave = 550;


figure(1); clf; set(gcf, 'Color', 'w'); hold all;
for ii = defocusMicrons
    
    wvf = wvfCreate('wave',calc_wave);%,'zcoeffs',zCoefs,'name',sprintf('human-%d',pupilMM));
    wvf = wvfSet(wvf,'calc pupil size',pupilMM);


    wvf = wvfSet(wvf,'zcoeffs',ii,{'defocus'});
    wvf = wvfComputePSF(wvf);

    otf = fftshift(wvfGet(wvf, 'otf', calc_wave));

    mRow = wvfGet(wvf, 'middle row');
    supPSF = wvfGet(wvf, 'psf spatial sample', 'mm', 550);
    supOTF = wvfGet(wvf, 'otfsupport', 'mm' ,calc_wave);
    
    % mm2deg
    supOTF = supOTF*0.3;

    
    plot(supOTF(mRow:end),abs(otf(mRow:end,mRow)), 'LineWidth', 4);
     ylim([0 1]); xlim([0 100]);
    xlabel(sprintf('Spatial freq (cycle/%s)',supportUnit))
    ylabel('MTF (a.u.)');
    box off;
    set(gca,'TickDir', 'out', 'FontSize', 30, 'XScale', 'log', 'XTickLabel',[10 100])
    

end

cutOffFreq_degPerCycle = ((10^6)*pi*pupilMM) / (180*calc_wave);
disp(cutOffFreq_degPerCycle)