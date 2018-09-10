%% s_numbersOnMeridonialDifferences

% Script to calculate the differences between cardinal meridians for cone
% density and midget retinal ganglion cell density

% Converter
deg2m   = 0.3 * 0.001;

% Polar angles
ang = [0 90 180 270]; % radians
angNames = {'Nasal (HM)','Superior (LVM)','Temporal (HM)', 'Inferior (UVM)'};

posVisualField_inDeg = 4.5;


%% Cone density (Sources in ISETBIO)
sourceNames = {'Song2011Young','Song2011Old','Curcio1990'};

for s = 1:numel(sourceNames)

    fprintf('Cone density ratio horizontal/vertical using data from %s:\n', sourceNames{s})
    
    for ii = 1:length(ang)
        [spacing(s, ii), aperture(s, ii), density(s, ii)] = coneSizeReadData('eccentricity',posVisualField_inDeg*deg2m, ...
            'angle',ang(ii),'eccentricityUnits', 'm','whichEye','left', 'coneDensitySource', sourceNames{s});
    end

    % display ratio horizontal / vertical
    disp((density(s, 1)+density(s, 3))/(density(s, 4)+density(s, 2)))
end 


%% RGC Midgets (Source from Barnet & Aguirre rgc Displacement code)

for pa = 1:length(ang)
    
    mRFDensitySqDegVisual(pa) = calcWatsonMidgetRFDensityByEccenDegVisual(posVisualField_inDeg, polarAngle(pa));

end

% Display ratio horizontal / vertival
fprintf('Midget RGC ratio horizontal/vertical using Watson (2014) equation 8:\n')
disp((mRFDensitySqDegVisual(1)+mRFDensitySqDegVisual(3)) / ...
    (mRFDensitySqDegVisual(2)+mRFDensitySqDegVisual(4)) )

fprintf('Midget RGC density horizontal and inferior retina:\n')
disp(mean([mRFDensitySqDegVisual(1),mRFDensitySqDegVisual(3)]))
disp(mRFDensitySqDegVisual(4))
fprintf('Midget RGC density ratio horizontal / inferior retina:\n')
disp(mean([mRFDensitySqDegVisual(1),mRFDensitySqDegVisual(3)]) / mRFDensitySqDegVisual(4))




return
%% OBSOLETE Getting RGC data (using Dacey data)

load(fullfile(isetbioDataPath,'rgc','midgetData.mat'))

% figure;
subplot(1,3,2);
scatter(midgetData(:,1),midgetData(:,2))

% Linear regression
midgetFit = ([ones(size(midgetData,1),1) midgetData(:,1)]\midgetData(:,2));
% Plot regression
hold on; plot(.1:.1:15,(.1:.1:15).*midgetFit(2)+midgetFit(1));

% Moving average bin
if movingAverageFlag
    clear xRange xRangeVal sizeAvg
    xRange = 1;
    for xRangeInd = 1:floor(max(midgetData(:,1))/xRange)
        xRangeVal(xRangeInd) = xRangeInd*xRange;
        xRangePts = find(abs(midgetData(:,1)-xRangeVal(xRangeInd))<xRange);
        if ~isempty(xRangePts)
            sizeAvg(xRangeInd) = mean(midgetData(xRangePts,2));
        end
        clear xRangePts
    end
    hold on;plot(xRangeVal,sizeAvg,'b')
end

% Format plot
title('Midget DF Size over Eccentricity, Dacey 2004');
xlabel('Eccentricity (mm)'); ylabel(sprintf('Dendritic field size (\\mum)'));
grid on;
set(gca,'fontsize',14);

legend('Data','Fit','Binned Average');
axis([0 18 0 450]);




%% OBSOLETE Getting RGC data (using isetbio mosaics)

imSize = 128; lineSpacing = 48; fov = 2; % deg
scene = sceneCreate('grid lines',imSize,lineSpacing);
scene = sceneSet(scene,'fov',fov);
oi = oiCreate;    % Standard human optics
oi = oiCompute(oi,scene);

nMovements = 25;
cMosaic = coneMosaic('center', [0 0.0015]);
cMosaic.setSizeToFOV(2);
cMosaic.emGenSequence(nMovements);
cMosaic.compute(oi, 'currentFlag', true);

bpL = bipolarLayer(cMosaic);
bpMosaicParams.spread  = 2;  % RF diameter w.r.t. input samples
bpMosaicParams.stride  = 2;  % RF diameter w.r.t. input samples
bpL.mosaic{1} = bipolarMosaic(cMosaic,'on midget',bpMosaicParams);
bpL.mosaic{1}.compute;

clear rgcL rgcParams

% Create RGC layer from bipolar layer
rgcL = rgcLayer(bpL);

% Spread and stride are not working
rgcParams.rfDiameter = 2;

% rgcL.mosaic{ii} = rgcGLM(rgcL, bpL.mosaic{1},'on midget');
rgcL.mosaic{1} = rgcGLM(rgcL, bpL.mosaic{1},'on midget',rgcParams);
rgcL.compute;

rgcL.window;