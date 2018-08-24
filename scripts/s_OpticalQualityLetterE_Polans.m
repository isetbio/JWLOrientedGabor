%% s_OpticalQualityLetterE_Polans

% Demonstrate optical quality with convolving an image of the letter E with
% a wavefront from a single observer in the Polans et al (2015) dataset.

% To simulate figure 5 from Jaeken & Artal:
% - Letter has size 30 arc min

% Clean up
clear all; close all;

%% 1. Define parameters

pupilMM     	= 4; % diameter in millimeters
waveToPlot      = 550; % nm
eccenHorz       = [5, 0, -5]; % degrees
eccenVert       = [5, 0, -5]; % degrees

% Locations of the wavefronts (x,y)
locations = {[0,5],[-5,0],[0,0],[5,0],[0,-5]}; 

% What specific subject? (1-10, except 2)
for subject  = [1,3:10]


%% 2. Loop over eccentricities to define wavefront

for ec = 1:numel(locations)
    
    clear wvf oi;
    
    % Get wavefront aberrations from Jaeken and Artal 2012 dataset
    [wvf] = wvfLoadWavefrontOpticsData('wvfZcoefsSource', 'Polans2015', 'jIndex', 3:14, 'whichEye', 'right', 'eccentricity', locations{ec}, 'whichGroup', subject, 'relativeRefraction', true, 'verbose', false);
    
    fprintf('Subject %d, Z4: %1.2f (D)\n', subject, wvfDefocusMicronsToDiopters(wvf.zcoeffs(5),4))
   
    %     wvf = wvfCreate('wave',wave,'zcoeffs',zCoefs,'name',sprintf('human-%d',pupilMM));
    %     wvf = wvfSet(wvf,'calc pupil size',pupilMM);
    %     wvf = wvfComputePSF(wvf);
    %
    % Compute OTF for given pupil and wavelength
    wvf = wvfSet(wvf,'calc pupil size',pupilMM);
    wvf = wvfSet(wvf,'calc wave',waveToPlot);

    PSF_toPlot(:,:,ec) =  wvf.psf{1};

end



%% 3. Create a letter image

% Get support to plot PSFs 
support_samples = wvfGet(wvf, 'psfangularsamples', 'min', 550);
support_minPerSample = wvfGet(wvf, 'psf arcmin per sample', 550);

% Point spread function
psf.pix = length(support_samples); %201
psf.min = max(support_samples)*2; %33 arc min

% Define image size of letter E with the same samples
E.min       = 50;
E.fraction  = .84;
E.pix       = round(E.min / support_minPerSample);

% Plot figure with and without grid
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

%% 4. Plot PSF and image of E after convolving with PSF
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
subplotTitles = {'Superior retina', 'Nasal retina', 'Fovea', 'Temporal retina', 'Inferior retina'};

for ecIdx = 1:numel(locations)
    
    % Convolve PSF * letter E
    opticalQuality = conv2(double(E.im),double(PSF_toPlot(:,:,ecIdx)),'same');
    
    set(0,'currentfigure', fH2)
    subplot(3,4,subplotLoc(ecIdx))
    imagesc(E.support, E.support, opticalQuality);
    colormap gray;
    xlabel('arc min'); ylabel('arc min'); title(sprintf('%s (X:%d, Y:%d)',subplotTitles{ecIdx}, locations{ecIdx}(1),locations{ecIdx}(2)))
    set(gca, 'TickDir', 'out'); box off;
    axis equal
    axis([-1 1 -1 1]*E.min/2);
    
    set(0,'currentfigure', fH3)
    subplot(3,4,subplotLoc(ecIdx))
    imagesc(support_samples, support_samples, PSF_toPlot(:,:,ecIdx));
    xlabel('arc min'); ylabel('arc min'); title(sprintf('%s (X:%d, Y:%d)',subplotTitles{ecIdx}, locations{ecIdx}(1),locations{ecIdx}(2)))
    set(gca, 'TickDir', 'out'); box off;
    axis equal
    axis([-1 1 -1 1]*max(support_samples));
    %     idx = idx+1;
    
    
end

pth = '/Volumes/GoogleDrive/My Drive/prep_PF_ISETBIO/figures/Fig3_Introduction_OpticalQuality/Polans_subjects/';
hgexport(fH2, fullfile(pth, sprintf('letterE_subject%d.eps',subject)))
pth = '/Volumes/GoogleDrive/My Drive/prep_PF_ISETBIO/figures/Fig3_Introduction_OpticalQuality/Polans_subjects/';
hgexport(fH3, fullfile(pth, sprintf('psf_subject%d.eps',subject)))
end