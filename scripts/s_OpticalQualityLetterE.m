
%% 1. Create a letter image
fH = figure(99);
t = text(0.45, 0.5,'E','FontName','Arial', 'FontUnits', 'centimeters'); axis off;
s = t.FontSize; 
t.FontSize = 3;
im = getframe(fH);
close(99)

%% 2. Define OTF parameters
pupilMM     = 3; % diameter in millimeters
waveToPlot  = 550; % nm
eccen       = [0, -2, -5]; % degrees

[centralRefraction, ~, ~, group] = wvfSortSubjectDataJaekenArtal2012;
whichGroup  = group(2).idxRE(randi(length(group(2).idxRE),1));% 'emmetropes'; % can also be myopes



mm2deg = 1000/(0.3*0.001);

[wdth, hght, cl] = size(im.cdata);
half_wdth = wdth/2;
half_hght = hght/2;
x = (half_wdth-201):(half_wdth+201);
y = (half_hght-201):(half_hght+201);



% Convert letter to gray scale and double format, downsample
centeredLetter = double(rgb2gray(im.cdata(x,y,:)));
centeredLetterDownsampled = centeredLetter(1:2:end,1:2:end);

figure(1); clf;
subplot(241);
imagesc(centeredLetter); colormap gray
xlabel('Pixels'); ylabel('Pixels'); title(whichGroup)
set(gca, 'TickDir', 'out'); box off;

idx = 1;

for ec = eccen
    
    % Get wavefront aberrations from Jaeken and Artal 2012 dataset
    [wvf, oi] = wvfLoadWavefrontOpticsData('jIndex', 0:14, 'whichEye','right', 'eccentricity',ec, 'whichGroup', whichGroup);
    
    % Compute OTF for given pupil and wavelength
    wvf = wvfSet(wvf,'calc pupil size',pupilMM);
    wvf = wvfSet(wvf,'calc wave',waveToPlot);
    
    % Get centered PSF
%     PSF_centered = wvfGet(wvf,'psf centered',550);
    PSF_centered = wvf.psf{1}([101:end, 1:100],:);
    x_mm = wvfGet(wvf,'psf spatial samples', 'mm', 550); % in mm
    x_deg = x_mm.*mm2deg;
    
    % Convolve PSF with letter
    opticalQuality = conv2(PSF_centered,centeredLetter,'full');
    
    
    ctr = ceil(0.5*size(opticalQuality,1));
    halfwidth = floor(0.5*length(x_deg));
    opticalQualityCropped = opticalQuality((ctr-halfwidth):(ctr+halfwidth),(ctr-halfwidth):(ctr+halfwidth),1);
    
    subplot(2,4,idx+1)
    imagesc(x_deg, x_deg, opticalQualityCropped);
    colormap gray;
    xlabel('deg'); ylabel('deg'); title(ec)
    set(gca, 'TickDir', 'out'); box off;
    
    subplot(2,4,idx+5);
    imagesc(x_deg, x_deg,PSF_centered);
    xlabel('deg'); ylabel('deg');  title(ec)
    set(gca, 'TickDir', 'out'); box off;
    
    

    
    idx = idx+1;
    
    
end