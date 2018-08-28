% Preprocess Whole Field Polans et al. 2015 (Optica) data

% Some details on the dataset:

% We measured the right eye. Along the horizontal meridian negative values
% represent temporal retina and positive values represent nasal retina.
% Along the vertical meridian positive and negative represent respectively
% superior and inferior retina.
% 
% The data attached are data in excel using a comma to indicate decimals.
% The data are in micron calculated for a 4 mm pupil. The order of the
% Zernike is conform ANSI2004. The data is the raw data.
% 
% Along the horizontal axis the range between +/- 35 degrees is reliable.
% The extreme 5 degrees the pupil is in various subjects small and could
% generate erroneous results especially for the higher order aberrations.
% 
% For the calculation of the Zernike coefficients in an elliptic pupil a
% central circular pupil of 4 mm is used.
% 
% A chromatic off-set was taken into account on Z4 (-0.8 D).
% 
% The subjects were between 30 and 40 years and were emmetropic or mild
% myopic. The measurements were done without cycloplegia. The target was a
% red laser point displayed on a wall at 2 m distance.
% 
% The results are very consistent along the horizontal meridian compared to
% the earlier 100 subjects data. The impact of the vertical meridian is
% especially visible for Z3 where it turns around the angle kappa, and for
% Z8 where there is hardly any variation compared to the large horizontal
% variation. In general the variation along the vertical meridian is minor
% than the variation along the horizontal meridian.


pth = '/Volumes/server/Projects/PerformanceFieldsIsetBio/literature/literatureData/';
fname = 'whole_field_10_subjects.xlsx';

numSubjects = 10;
verticalAxis = -25:5:25;
horizontalAxis = -40:1:40;
zernikeCoeff = 3:20;


data = NaN(numSubjects, (length(verticalAxis)*length(horizontalAxis)), length(zernikeCoeff)+2); 

for s = 1:numSubjects
    data(s,:,:) = xlsread(fullfile(pth, fname), s);
end

% Strip axes from data
data = data(:,2:end, 3:end);

