
%% 1. Create a letter image
fH = figure;
t = text(0.4,0.5,'E','FontName','Arial'); axis off;
s = t.FontSize; colormap(gray);
t.FontSize = 150;
im = getframe(fH);


%% 2. Define OTF parameters
pupilMM     = 3; % diameter in millimeters
waveToPlot  = 550; % nm
eccen       = 0; % degrees
whichGroup  = 'emmetropes'; % can also be myopes

mm2deg = 1000/(0.3*0.001);


% Get wavefront aberrations from Jaeken and Artal 2012 dataset 
[wvf, oi] = wvfLoadWavefrontOpticsData('jIndex', 0:14, 'whichEye','left', 'eccentricity',eccen, 'whichGroup', whichGroup);

% Compute OTF for given pupil and wavelength
wvf = wvfSet(wvf,'calc pupil size',pupilMM);
wvf = wvfSet(wvf,'calc wave',waveToPlot);

% Center OTF
% wvf.otf{1} = fftshift(wvf.otf{1});

% Get centered PSF
PSF_centered = wvfGet(wvf,'psf centered',550);
x_mm = wvfGet(wvf,'psf spatial samples', 'mm', 550); % in mm
x_deg = x_mm.*mm2deg;

centeredLetter = double(rgb2gray(im.cdata(1:402,100:501,:)));
centeredLetterDownsampled = centeredLetter(1:2:end,1:2:end);


% imagesc(data.x, data.y, centeredLetter); colormap gray


% Convolve PSF with letter
opticalQuality = conv2(PSF_centered,centeredLetter,'full');


ctr = ceil(0.5*size(opticalQuality,1));
halfwidth = floor(0.5*length(x_arcmin));
opticalQualityCropped = opticalQuality((ctr-halfwidth):(ctr+halfwidth),(ctr-halfwidth):(ctr+halfwidth),1);

figure(1); clf;
subplot(131)
imagesc(centeredLetterDownsampled);
xlabel('Pixels'); ylabel('Pixels'); 
set(gca, 'TickDir', 'out'); box off;


subplot(132);
imagesc(x_deg, x_deg,PSF_centered);
xlabel('deg'); ylabel('deg'); 
set(gca, 'TickDir', 'out'); box off;


subplot(133)
imagesc(x_deg, x_deg, opticalQualityCropped); 
colormap gray;
xlabel('deg'); ylabel('deg'); 
set(gca, 'TickDir', 'out'); box off;
