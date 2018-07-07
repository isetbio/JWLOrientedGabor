


%% 2. Define OTF parameters
pupilMM     = 4; % diameter in millimeters
waveToPlot  = 550; % nm
eccen       = [-5 0 5]; % degrees

[centralRefraction, ~, ~, group] = wvfSortSubjectDataJaekenArtal2012;
whichGroup  = 38;%group(2).RE(randi(length(group(2).RE),1));% 'emmetropes'; % can also be myopes

horizontalDF = wvfDefocusDioptersToMicrons([-0.10, -0.02, -0.20], 4);
verticalDF   = wvfDefocusDioptersToMicrons([0.03, NaN, 0.01],4);

% pupilMM = 4.5;
% zCoefs = wvfLoadThibosVirtualEyes(pupilMM);
% wave = 550;


%% Loop over eccentricities to define wavefront and convolve with letter image
for ec = eccen
    
    clear wvf;
    
    % Get wavefront aberrations from Jaeken and Artal 2012 dataset
    [wvf, oi] = wvfLoadWavefrontOpticsData('jIndex', 0:14, 'whichEye','right', 'eccentricity',ec, 'whichGroup', whichGroup);
    
    %     wvf = wvfCreate('wave',wave,'zcoeffs',zCoefs,'name',sprintf('human-%d',pupilMM));
    %     wvf = wvfSet(wvf,'calc pupil size',pupilMM);
    %     wvf = wvfComputePSF(wvf);
    %
    % Compute OTF for given pupil and wavelength
    wvf = wvfSet(wvf,'calc pupil size',pupilMM);
    wvf = wvfSet(wvf,'calc wave',waveToPlot);
    
    if ec ~= 0
        wvf =  wvfSet(wvf,'zcoeffs', verticalDF(ec==eccen), {'defocus'});
        %
        %     % compute the PSF
        wvf = wvfComputePSF(wvf);
        otf = wvfGet(wvf,'otf');
    else
        otf = oi.optics.OTF.OTF;
        
    end
    
    %     this_psf = wvfGet(wvf, 'psf');
    %     PSF_centered(:,:,ec==eccen) = this_psf;%([101:end, 1:100],:);
    
    
    % Get centered PSF
    %     PSF_centered(:,:,ec==eccen) = wvfGet(wvf,'psf centered',550);
    %     PSF_centered(:,:,ec==eccen) = wvf.psf{1}([101:end, 1:100],:);
    
    psfFromOI = abs(ifft2((otf))); %abs(ifft2((oi.optics.OTF.OTF)));
    
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


E.min       = 50;
E.fraction  = .84;
E.pix       = round(E.min / support_minPerSample);

figure(1); clf;set(gcf, 'Color', 'w')
axis([0 1 0 1]*E.min); axis square;
t = text( .5*E.min, .5*E.min, 'E');
t.FontUnits = 'normalized';
t.FontSize = E.fraction; % .84;
t.VerticalAlignment = 'middle';
t.HorizontalAlignment = 'center';
t.FontName = 'Arial';
t.FontWeight = 'bold';

axis off
grid off

E1 = getframe(gca);
E.im = rgb2gray(E1.cdata);
E.im = imresize(E.im, E.pix * [1 1]);
E.support = linspace(-E.min/2, E.min/2, E.pix);

axis on
grid on

E2 = getframe(gca);
E2.im = rgb2gray(E2.cdata);
E2.im = imresize(E2.im, E.pix * [1 1]);

close(1)

%% Convolve PSF with letter
figure(2); clf; hold all
for ecIdx = 1:length(eccen)
    opticalQuality = conv2(double(E.im),double(PSF_centered(:,:,ecIdx)),'same');
    
    subplot(241);
    imagesc(E.support, E.support, E.im); colormap gray
    xlabel('arc min'); ylabel('arc min'); %title(sprintf('Subject: %d', whichGroup))
    set(gca, 'TickDir', 'out'); box off; axis equal
    axis([-1 1 -1 1]*E.min/2);
    
    subplot(2,4,ecIdx+1)
    imagesc(E.support, E.support, opticalQuality);
    colormap gray;
    xlabel('arc min'); ylabel('arc min'); title(sprintf('Eccentricity: %d',eccen(ecIdx)))
    set(gca, 'TickDir', 'out'); box off;
    axis equal
    axis([-1 1 -1 1]*E.min/2);
    
    subplot(2,4,ecIdx+5);
    imagesc(support_samples, support_samples, PSF_centered(:,:,ecIdx));
    xlabel('arc min'); ylabel('arc min'); title(sprintf('Eccentricity: %d',eccen(ecIdx)))
    set(gca, 'TickDir', 'out'); box off;
    axis equal
    axis([-1 1 -1 1]*max(support_samples));
    %     idx = idx+1;
    
    
end

