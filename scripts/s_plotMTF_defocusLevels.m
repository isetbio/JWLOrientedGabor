%% s_plotMTF_defocusLevels

% Script to plot the modulation transfer function (MTF) for different
% defocus levels

%% 0. Define parameters

pupilMM        = 3;         % pupil size (diameter) in mm
defocusMicrons = [0, 1, 2]; % defocus zernike coefficient in microns
calc_wave      = 550;       % wavelength used to calculate the MTF

cmap           = [0 0 0; 1 0 0; 1 0.7812 0.4975];

%% 1. Get zernike coefficients
zCoefs = wvfLoadThibosVirtualEyes(pupilMM);

%% 2. Plot MTFs

figure(100); clf; set(gcf, 'Color', 'w', 'NumberTitle', 'off', 'Name', 'Figure 8A - MTF for different defocus levels');
hold all;

for ii = defocusMicrons
    
    % Create a wavefront with pupil size and set wavelength of interest
    wvf = wvfCreate('wave',calc_wave); %,'zcoeffs',zCoefs,'name',sprintf('human-%d',pupilMM));
    wvf = wvfSet(wvf,'calc pupil size',pupilMM);

    % Set zernike coefficients and calculate point spread function (PSF)
    wvf = wvfSet(wvf,'zcoeffs',ii,{'defocus'});
    wvf = wvfComputePSF(wvf);

    % MTF (same as OTF, optical transfer function)
    otf = fftshift(wvfGet(wvf, 'otf', calc_wave));

    % Get support
    mRow = wvfGet(wvf, 'middle row');
    supPSF = wvfGet(wvf, 'psf spatial sample', 'mm', 550);
    supOTF = wvfGet(wvf, 'otfsupport', 'mm' ,calc_wave);
    
    % Convert mm 2 deg
    supOTF = supOTF*0.3;

    % Plot it!
    plot(supOTF(mRow:end),abs(otf(mRow:end,mRow)), 'LineWidth', 4, 'Color', cmap((ii==defocusMicrons),:));
    ylim([0 1]); xlim([0 100]);
    xlabel('Spatial freq (cycle/deg)')
    ylabel('Normalized modulation (a.u.)');
    box off;
    set(gca,'TickDir', 'out', 'FontSize', 30, 'XScale', 'log', 'XTickLabel',[10 100])

end

% cutOffFreq_degPerCycle = ((10^6)*pi*pupilMM) / (180*calc_wave);
% disp(cutOffFreq_degPerCycle)

%% 3. Stimuli (2D and 1D) under different defocus levels

figure(99); clf; set(gcf, 'Color', 'w', 'Position', [1000, 681, 1173, 657], 'NumberTitle', 'off', 'Name', 'Figure 8A - 1D and 2D stim for different defocus levels');
hold all;

plotIdx = 1;

for ii = defocusMicrons
    
    % Scene field of view
    sparams.sceneFOV  = 2;   % scene field of view in degrees (diameter)
    sparams.freqCPD   = 4;   % Gabor spatial frequency (cpd)
    sparams.gausSDdeg = sparams.sceneFOV/4; % Gabor SD in degrees of visual angle

    % Unit converters
    deg2fov = 1/sparams.sceneFOV;
    fov2deg = sparams.sceneFOV;
    deg2m   = 0.3 * 0.001;          % (we first choose 3 deg per mm, .001 mm per meter, but now adjusted to the default in isetbio)

    % Gabor parameters
    sparams.gabor           = harmonicP;                   % Standard Gabor
    sparams.gabor.ang       = (pi/180)* 15;                % Gabor orientation (radians) - question: what is 0??
    sparams.gabor.freq      = fov2deg*sparams.freqCPD;     % Spatial frequency (cycles/FOV)
    sparams.gabor.contrast  = 1;                           % Presumably michelson, [0 1]
    sparams.gabor.GaborFlag = sparams.gausSDdeg*deg2fov;   % Gaussian window

    % Add defocus
    sparams.oi = oiDefocus(ii);

    % Create Optical Image Sequence
    OG = ogStimuli(sparams);

    irradiance = OG(4).oiModulated.data.illuminance;
    midpoint = ceil(size(irradiance,1)/2);

    subplot(2, 3, plotIdx);
    imagesc(irradiance); colormap gray; axis square;
    set(gca, 'CLim', [0 4]); axis off;
    title(sprintf('Defocus level: %d microns', ii))
    
    subplot(2,3, plotIdx+2+1);
    plot(irradiance(midpoint,:), 'LineWidth', 2, 'Color', cmap((ii==defocusMicrons),:))
    set(gca, 'TickDir', 'out', 'XTick', [midpoint + [-30, -15, 0, 15, 30]], 'XTickLabel', {'-1', '-0.5', '0', '0.5', '1'})
    ylabel('Illuminance'); xlabel('X Position (deg)'); box off;
    ylim([0 4]);
    
    plotIdx = plotIdx + 1;
end