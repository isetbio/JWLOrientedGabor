% s_plotArtalvsThibosDefocus.m


% close all;

%% 0. Define OTF parameters
pupilMM     = 3; % diameter in millimeters
waveToPlot  = 550; % nm
eccen       = 0; % degrees
whichGroup  = 'emmetropes'; % can also be myopes

% Get wavefront aberrations from Jaeken and Artal 2012 dataset 
[wvf, oi] = wvfLoadJaekenArtal2012Data('jIndex', 0:14, 'whichEye','left', 'eccentricity',eccen, 'whichGroup', whichGroup);

% Compute OTF for given pupil and wavelength
wvf = wvfSet(wvf,'calc pupil size',pupilMM);
wvf = wvfSet(wvf,'calc wave',waveToPlot);

% Plot OTF data
uData1(1) = wvfPlot(wvf,'OTF','mm',waveToPlot);
set(gca, 'XLim', [-150 150], 'YLim', [-150 150], 'zlim', [0 1.2],'FontSize',16,'LineWidth',2);

%% 1. Loop over different levels of defocus aberrations (in um, for defocus)
defocusLevels = 0:0.25:1.75; % um
for ii = 1:length(defocusLevels)
    
    clear wvfP;
    
    % Load Thibos modeled zernike coefficients
    zCoefs = wvfLoadThibosVirtualEyes(pupilMM);
    
    % Create wavefront
    wvfP = wvfCreate('wave',550,'zcoeffs',zCoefs,'name', sprintf('Thibos human-%d',pupilMM));
    wvfP = wvfSet(wvfP,'calc pupil size',pupilMM);
    wvfP = wvfSet(wvfP,'calc wave',waveToPlot);
    wvfP = wvfSet(wvfP,'zcoeffs',ii,{'defocus'});
    
    % Compute PSF in order to get OTF
    wvfP = wvfComputePSF(wvfP);
    uData1(ii+1) = wvfPlot(wvfP,'OTF','mm',waveToPlot);
end

%% 2. Get Jaeken & Artal wavefront OTF at different eccentricities
eccentricities = -8:0; % degrees, negative numbers correspond to temporal retina in left eye and nasal retina in right eye.

for eccen = 1:length(eccentricities)
    
    % Get wavefront aberrations from Jaeken and Artal 2012 dataset 
    [wvf, oi] = wvfLoadJaekenArtal2012Data('jIndex', 0:14, 'whichEye','left', 'eccentricity',eccentricities(eccen), 'whichGroup', whichGroup);

    % Compute OTF for given pupil and wavelength
    wvf = wvfSet(wvf,'calc pupil size',pupilMM);
    wvf = wvfSet(wvf,'calc wave',waveToPlot);

    % Plot OTF data
    uData2(eccen) = wvfPlot(wvf,'OTF','mm',waveToPlot);
    
end


close all;

%% Visualize:
% 1. Plot Jaeken & Artal wavefront OTF, against different levels of defocus with Thibos model

% Set up figure
vcNewGraphWin([],'upperleftbig')

for ii = 1:length(defocusLevels)+1
    subplot(3,3,ii);
    surf(uData1(ii).fx, uData1(ii).fy, uData1(ii).otf); view(2); axis square
    colormap parula; 
    if ii == 1; title('Average of Artal data at fovea');
    else, title(sprintf('Thibos: Defocus %1.1f Diopters', zernikeDefocus2diopter(defocusLevels(ii-1), pupilMM))); end
    
    set(gca,'xlim',[-150 150] ,'ylim',[-150 150], 'zlim', [0 1.2],'FontSize',16,'LineWidth',2);
    xlabel('freq (lines/mm)'); ylabel('freq (lines/mm)')
end

% 2. Plot Jaeken & Artal wavefront OTF at different eccentricities
vcNewGraphWin([],'upperleftbig')
for ii = 1:9
    subplot(3,3,ii); hold all;
    surf(uData2(ii).fx, uData2(ii).fy, uData2(ii).otf); view(2); axis square
    colormap parula; 
    title(sprintf('Artal: %1.0f deg eccen',eccentricities(ii)))
    set(gca,'xlim',[-150 150] ,'ylim',[-150 150], 'zlim', [0 1.2],'FontSize',16,'LineWidth',2);
    xlabel('freq (lines/mm)'); ylabel('freq (lines/mm)');
end
