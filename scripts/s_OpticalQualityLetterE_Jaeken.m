
clear all; close all;

%% 2. Define OTF parameters
pupilMM     = 4; % diameter in millimeters
waveToPlot  = 550; % nm
eccen       = [-5 0 5]; % degrees

[centralRefraction, ~, ~, group] = wvfSortSubjectDataJaekenArtal2012;
whichGroup  = 38;%group(2).RE(randi(length(group(2).RE),1));% 'emmetropes'; % can also be myopes

 meridian = 'horz';

% pupilMM = 4.5;
% zCoefs = wvfLoadThibosVirtualEyes(pupilMM);
% wave = 550;

% Z Coefficients for subject 38, for the three horizontal eccentricities

% Eccen -5, right eye, therefore nasal retina
% zcoefs(1,:) = [-0.2848, -5.3751, -0.7814, 0.1268, -0.0754, -0.0644, -0.0474, 0.1217, 0.0306, -0.0475, 0.0109, -0.0185, 0.0688, -0.0024, 0.0235];

% Eccen 0
% zcoefs(2,:) = [-0.4365, -5.4364, -0.8595, 0.0814, -0.1639, -0.1605, -0.0502, 0.1213, 0.0444, -0.0485, 0.0118, -0.0131, 0.0683, -0.0062, 0.0190];

% Eccen 5, right eye, therefore temporal retina
% zcoefs(3,:) = [-0.5445, -5.5653, -1.0369, 0.0779, -0.2309, -0.2028, -0.0532, 0.1199, 0.0491, -0.0454, 0.0131, -0.0132, 0.0635, -0.0072, 0.0192];

% We can eyeball the subject average of the horizontal and vertical level of Defocus (D) in the Polans
% paper, in order to adjust the wavefront for both vertical and horizontal axis
    % horizontalGroupAverageDF = wvfDefocusDioptersToMicrons([-0.10, -0.02, -0.20], 4);
    % verticalGroupAverageDF   = wvfDefocusDioptersToMicrons([0.03, NaN, 0.01],4); % we assume 4 mm pupil
% verticalGroupAverageDF   = [-0.1, NaN, -0.08];
%% Loop over eccentricities to define wavefront and convolve with letter image

% switch meridian
%     case 'horz'
%         % do nothing -- keep individual subject data, not group average
%         % data
% %         thisDF = horizontalGroupAverageDF;
%     case 'vert' 
%         % Change defocus component to group averaged defocus
%         thisDF = verticalGroupAverageDF;
% end
         

for ec = eccen
    
    clear wvf;
    
    % Get wavefront aberrations from Jaeken and Artal 2012 dataset
    [wvf, oi] = wvfLoadWavefrontOpticsData('jIndex', 0:14, 'whichEye','right', 'eccentricity',[ec 0], 'whichGroup', whichGroup, 'relativeRefraction', true, 'verbose', false);
    
    %     wvf = wvfCreate('wave',wave,'zcoeffs',zCoefs,'name',sprintf('human-%d',pupilMM));
    %     wvf = wvfSet(wvf,'calc pupil size',pupilMM);
    %     wvf = wvfComputePSF(wvf);
    %
    % Compute OTF for given pupil and wavelength
    wvf = wvfSet(wvf,'calc pupil size',pupilMM);
    wvf = wvfSet(wvf,'calc wave',waveToPlot);
    
%     if (ec ~= 0) && (strcmp(meridian,'vert'))
%         
%         % Set defocus zernike coefficient to group average, then
%         % recalculate the wave front
% %         zcoefs(ec==eccen, 5) = thisDF(ec==eccen);
%         wvf =  wvfSet(wvf,'zcoeffs', zcoefs(ec==eccen,:), 0:14);
% 
%         % Recompute the PSF
%         wvf = wvfComputePSF(wvf);
%         wvf.psf{1} = wvfGet(wvf,'psf centered',550);
        PSF_centered(:,:,ec==eccen) =  wvf.psf{1};
% 
%         % Recompute the OTF
%         %	otf = wvfGet(wvf,'otf');
% 
% %     else % If horizontal: We already got the OTF from wvfLoadWavefrontOpticsData function
%         otf = oi.optics.OTF.OTF;
%         psfFromOI = abs(ifft2((otf))); %abs(ifft2((oi.optics.OTF.OTF)));
%         PSF_centered(:,:,ec==eccen) = psfFromOI([101:end, 1:100],[101:end, 1:100]);    
%     end
    

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

% axis on
% grid on
% 
% E2 = getframe(gca);
% E2.im = rgb2gray(E2.cdata);
% E2.im = imresize(E2.im, E.pix * [1 1]);

close(1)

%% Convolve PSF with letter
figure(2); clf; hold all
for ecIdx = 1:length(eccen)
    opticalQuality = conv2(double(E.im),double(PSF_centered(:,:,ecIdx)),'same');
    
    subplot(241);
    imagesc(E.support, E.support, E.im); colormap gray
    xlabel('arc min'); ylabel('arc min'); title(sprintf('Original image for subject: %d', whichGroup))
    set(gca, 'TickDir', 'out'); box off; axis equal
    axis([-1 1 -1 1]*E.min/2);
    
    subplot(2,4,ecIdx+1)
    imagesc(E.support, E.support, opticalQuality);
    colormap gray;
    xlabel('arc min'); ylabel('arc min'); title(sprintf('Eccentricity: %d at %s',eccen(ecIdx), meridian))
    set(gca, 'TickDir', 'out'); box off;
    axis equal
    axis([-1 1 -1 1]*E.min/2);
    
    subplot(2,4,ecIdx+5);
    imagesc(support_samples, support_samples, PSF_centered(:,:,ecIdx));
    xlabel('arc min'); ylabel('arc min'); title(sprintf('Eccentricity: %d at %s',eccen(ecIdx), meridian))
    set(gca, 'TickDir', 'out'); box off;
    axis equal
    axis([-1 1 -1 1]*max(support_samples));
    %     idx = idx+1;
    
    
end

