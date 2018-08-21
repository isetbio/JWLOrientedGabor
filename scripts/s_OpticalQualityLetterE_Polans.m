
clear all; close all;

%% 2. Define OTF parameters
pupilMM     = 4; % diameter in millimeters
waveToPlot  = 550; % nm
eccenHorz       = [5, 0, -5]; % degrees
eccenVert       = [5, 0, -5]; % degrees

locations = {[0,5],[-5,0],[0,0],[5,0],[0,-5]}; 

subject  = 10;




%% Loop over eccentricities to define wavefront and convolve with letter image

for ec = 1:numel(locations)
    
    clear wvf;
    
    % Get wavefront aberrations from Jaeken and Artal 2012 dataset
    [wvf, oi] = wvfLoadWavefrontOpticsData('wvfZcoefsSource', 'Polans2015', 'jIndex', 3:14, 'whichEye', 'right', 'eccentricity', locations{ec}, 'whichGroup', subject, 'verbose', false);
    disp(wvfDefocusMicronsToDiopters(wvf.zcoeffs,4))
    %     wvf = wvfCreate('wave',wave,'zcoeffs',zCoefs,'name',sprintf('human-%d',pupilMM));
    %     wvf = wvfSet(wvf,'calc pupil size',pupilMM);
    %     wvf = wvfComputePSF(wvf);
    %
    % Compute OTF for given pupil and wavelength
    wvf = wvfSet(wvf,'calc pupil size',pupilMM);
    wvf = wvfSet(wvf,'calc wave',waveToPlot);
    


    % Recompute the PSF
%     wvf = wvfComputePSF(wvf);
%      wvf.psf{1} = wvfGet(wvf,'psf centered',550);
%     if ec ==2;
%         otf = wvfGet(wvf,'otf');
%         wvf.otf{1} = fftshift(otf);
%         wvf.psf = {};
%         wvf = wvfComputePSF(wvf);
%         PSF_centered(:,:,ec) =  wvf.psf{1}(:,[101:end, 1:100]);
%     else
        PSF_centered(:,:,ec) =  wvf.psf{1};
%     end

    % Recompute the OTF
    %	otf = wvfGet(wvf,'otf');
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
fH2 = figure(2); clf; hold all

subplot(3,4,1);
imagesc(E.support, E.support, E.im); colormap gray
xlabel('arc min'); ylabel('arc min'); title(sprintf('Original image for subject: %d', subject))
set(gca, 'TickDir', 'out'); box off; axis equal
axis([-1 1 -1 1]*E.min/2);

fH3 = figure(3); clf; hold all
subplot(3,4,1);
imagesc(E.support, E.support, E.im); colormap gray
xlabel('arc min'); ylabel('arc min'); title(sprintf('Original image for subject: %d', subject))
set(gca, 'TickDir', 'out'); box off; axis equal
axis([-1 1 -1 1]*E.min/2);

subplotLoc = [3,6,7,8,11];

for ecIdx = 1:numel(locations)
    opticalQuality = conv2(double(E.im),double(PSF_centered(:,:,ecIdx)),'same');
    
    set(0,'currentfigure', fH2)
    subplot(3,4,subplotLoc(ecIdx))
    imagesc(E.support, E.support, opticalQuality);
    colormap gray;
    xlabel('arc min'); ylabel('arc min'); title(sprintf('Eccentricity: X: %d, Y: %d',locations{ecIdx}(1),locations{ecIdx}(2)))
    set(gca, 'TickDir', 'out'); box off;
    axis equal
    axis([-1 1 -1 1]*E.min/2);
    
    set(0,'currentfigure', fH3)
    subplot(3,4,subplotLoc(ecIdx))
    imagesc(support_samples, support_samples, PSF_centered(:,:,ecIdx));
    xlabel('arc min'); ylabel('arc min'); title(sprintf('Eccentricity: X: %d, Y: %d',locations{ecIdx}(1),locations{ecIdx}(2)))
    set(gca, 'TickDir', 'out'); box off;
    axis equal
    axis([-1 1 -1 1]*max(support_samples));
    %     idx = idx+1;
    
    
end

