% s_numbersOnMerodionalDifferences

% Calculate the differences between cardinal meridians for cone density,
% retinal ganglion cell density

deg2m   = 0.3 * 0.001;

 

ang = [0 90 180 270];

angNames = {'Nasal (HM)','Superior (LVM)','Temporal (HM)', 'Inferior (UVM)'};

 

for ii = 1:4

    [s(ii), a(ii), d(ii)] = coneSizeReadData('eccentricity',4.5*deg2m,...

        'angle',ang(ii),'eccentricityUnits', 'm','whichEye','left', 'coneDensitySource', 'Song2011Young');

end

 

% ratio: hor / ver

disp((d(1)+d(3))/(d(4)+d(2)))

 

for ii = 1:4

    [s(ii), a(ii), d(ii)] = coneSizeReadData('eccentricity',4.5*deg2m,...

        'angle',ang(ii),'eccentricityUnits', 'm','whichEye','left', 'coneDensitySource', 'Song2011Old');

end

 

% ratio: hor / ver

disp((d(1)+d(3))/(d(4)+d(2)))

 

for ii = 1:4

    [s(ii), a(ii), d(ii)] = coneSizeReadData('eccentricity',4.5*deg2m,...

        'angle',ang(ii),'eccentricityUnits', 'm','whichEye','left', 'coneDensitySource', 'Curcio1990');

end

 

% ratio: hor / ver

disp((d(1)+d(3))/(d(4)+d(2)))