% s_numbersOnMerodionalDifferences

% Calculate the differences between cardinal meridians for cone density,
% retinal ganglion cell density

% Converter
deg2m   = 0.3 * 0.001;

% Polar angles
ang = [0 90 180 270]; % radians
angNames = {'Nasal (HM)','Superior (LVM)','Temporal (HM)', 'Inferior (UVM)'};

% Cone density source names
sourceNames = {'Song2011Young','Song2011Old','Curcio1990'};


for jj = 1:numel(sourceNames)

    fprintf('Ratio hor/ver for %s:\n', sourceNames{jj})
    
    for ii = 1:length(ang)

        [s(jj, ii), a(jj, ii), d(jj, ii)] = coneSizeReadData('eccentricity',4.5*deg2m, ...
            'angle',ang(ii),'eccentricityUnits', 'm','whichEye','left', 'coneDensitySource', sourceNames{jj});

    end

    % ratio: hor / ver
    disp((d(jj, 1)+d(jj, 3))/(d(jj, 4)+d(jj, 2)))
end 


%% RGC Midgets (using Dacey data)

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