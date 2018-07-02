%% s_OpticalQualityLetterE

% Script to illustrate the effect of optical quality on an example stimulus
% (Letter 'E') for different eccentricities

%% 1. Define OTF parameters
pupilMM     = 4; % diameter in millimeters
waveToPlot  = 550; % nm
eccen       = [0, -2, -5]; % degrees

[centralRefraction, ~, ~, group] = wvfSortSubjectDataJaekenArtal2012;
whichGroup  = group(2).LE(randi(length(group(2).LE),1));% 'emmetropes'; % can also be myopes

%% Loop over eccentricities to define wavefront and convolve with letter image
for ec = eccen
    
    clear wvf;
    
    % Get wavefront aberrations from Jaeken and Artal 2012 dataset
    [wvf, oi] = wvfLoadWavefrontOpticsData('jIndex', 0:14, 'whichEye','left', 'eccentricity',ec, 'whichGroup', whichGroup);
    
    % Compute OTF for given pupil and wavelength
    wvf = wvfSet(wvf,'calc pupil size',pupilMM);
    wvf = wvfSet(wvf,'calc wave',waveToPlot);
    
    % Get centered PSF
%     PSF_centered = wvfGet(wvf,'psf centered',550);
%     PSF_centered = wvf.psf{1}([101:end, 1:100],:);
    
    psfFromOI = abs(ifft2((oi.optics.OTF.OTF)));
    PSF_centered(:,:,ec==eccen) = psfFromOI([101:end, 1:100],[101:end, 1:100]);
    
    
end
    
support_samples = wvfGet(wvf, 'psfangularsamples', 'min', 550);    
support_minPerSample = wvfGet(wvf, 'psf arcmin per sample', 550);    


%% 2. Create a letter image

% To simulate figure 5 from Jaeken & Artal:
% - Letter has size 30 arc min

% Point spread function 
psf.pix = length(support_samples); %201
psf.min = max(support_samples)*2; %33 arc min

minperpix = psf.min /  psf.pix;

E.min       = 100;
E.fraction  = .33;
E.pix       = round(E.min / minperpix);

figure(1); clf;set(gcf, 'Color', 'w')
axis([0 1 0 1]*E.min); axis square;
t = text( .5*E.min, .5*E.min, 'E');
t.FontUnits = 'normalized';
t.FontSize = E.fraction; % .84;
t.VerticalAlignment = 'middle';
t.HorizontalAlignment = 'center';
t.FontName = 'Arial';

axis off
grid off

E1 = getframe(gca);
E.im = rgb2gray(E1.cdata);
E.im = imresize(E.im, E.pix * [1 1]);

axis on
grid on

E2 = getframe(gca);
E2.im = rgb2gray(E2.cdata);
E2.im = imresize(E2.im, E.pix * [1 1]);


    
%% Convolve PSF with letter
figure(1); clf; hold all

for ecIdx = 1:length(eccen)
    
    opticalQuality = conv2(PSF_centered(:,:,ecIdx),E.im,'full');

    ctr = ceil(0.5*size(opticalQuality,1));
    halfwidth = floor(0.5*length(support_samples));
    opticalQualityCropped = opticalQuality((ctr-halfwidth):(ctr+halfwidth),(ctr-halfwidth):(ctr+halfwidth),1);

    subplot(241);
    imagesc(support_samples,support_samples,E.im); colormap gray
    xlabel('arc min'); ylabel('arc min');
    set(gca, 'TickDir', 'out'); box off;

    subplot(2,4,ecIdx+1)
    imagesc(support_samples, support_samples, opticalQualityCropped);
    colormap gray;
    xlabel('arc min'); ylabel('arc min'); title(sprintf('Eccentricity: %d',eccen(ecIdx)))
    set(gca, 'TickDir', 'out'); box off;

    subplot(2,4,ecIdx+5);
    imagesc(support_samples, support_samples, PSF_centered(:,:,ecIdx));
    xlabel('arc min'); ylabel('arc min'); title(sprintf('Eccentricity: %d',eccen(ecIdx)))
    set(gca, 'TickDir', 'out'); box off;

end

