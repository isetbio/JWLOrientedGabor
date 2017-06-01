function [OG, scenes, tseries, fname] = ogStimuli(varargin)
% OGSTIMULI - Create oriented stimuli (e.g., left and right)
%
%  [OG, scenes, tseries, fname] = ogStimuli(varargin)
%
% There is one input argument that is a struct with these parameters
%
%  ogabor    - Parameters for the oriented gabor stimuli; Default is ogaborP
%  tsamples  - Time samples (sec) ??
%  timesd    - Time standard deviation
%  display   - Display struct, default displayCreate('LCD-Apple')
%  sceneFOV  - Degrees, default 2 deg
%  distance  - Meters, viewing distance to display, default is 0.57 m
%
%  When the LCD-Apple display is [200, 200] pixels and the sceneFOV is 2 deg,
%  then 1 pixel is 0.6 arc min.
%
%  ?? See also s_EIParameters, s_EIMosaicSize ??
%
% The default display is an Apple LCD monitor.
%
% Outputs
%
%  OG       - 1x2 OIsequence array for left- and right-oriented Gabors
%  scenes   - A 1x2 cell with a uniform field scene and oriented gabor scene
%  tseries  - A vector of the temporal window modulator for the scene
%  fname    - Path and filename where stimulus is stored
%
% Example: 
%     params.cmFOV              = 2;
%     params.em                 = emCreate;
%     params.tStep              = 0.005;      % seconds
%     params.sceneFOV           = 2;
%     params.tsamples           = (-0.060:params.tStep:0.070); % seconds
%     params.timesd             = 0.100;                      % seconds
%     params.gabor              = harmonicP; % default harmonic
%     params.gabor.ang          = (pi/180)* 20;  % Gabor orientation (radians)
%     params.gabor.freq         = 6*params.sceneFOV;
%     params.gabor.GaborFlag    = .25/params.sceneFOV; % gaussian sd for gabor
%     [OG,scenes,tseries, fname] = ogStimuli(params);
% 
%     OG(1).visualize
%     ieAddObject(scenes{2}); sceneWindow;
%     vcNewGraphWin; plot(tseries)
%
% Notes
%   Luminance 2 log millilamberts
%
% JW & EK, NYU ISETBIO Team, 2017
%
% See also vaStimuli in WLVernier repository

%% Load a stored image, if it exists

pth = fullfile(ogRootPath, 'Stimulus');
if ~exist(pth, 'dir'), mkdir(pth); end
fname = fullfile(pth, 'stimulus.mat');

% if exist(fname,'file')
%     disp('Loading stimulus from file - parameters match')
%     try
%         load(fname,'OG','scenes','tseries');
%         return;
%     catch
%         disp('File found, but not the variables.  Creating.')
%     end
% else 
%     disp('Creating and saving stimulus file - no match found')
% end


%% Parse inputs
p = inputParser;

p.KeepUnmatched = true;   % Sometimes we overload p for SVM and cMosaic

p.addParameter('gabor', harmonicP,@isstruct);
p.addParameter('tsamples', (-50:100),@isvector);  % Time samples (ms)
p.addParameter('timesd', 20,@isscalar);           % Time standard deviation
p.addParameter('sceneFOV',2,@isscalar);           % Degrees
p.addParameter('distance',0.57,@isscalar);        % Meters
p.addParameter('bgColor',0.5,@isscalar);          % 0 to 1 (assumed grayscale) | currently not used

p.parse(varargin{:});

%% Derive variables from parsed inputs
oGabor     = p.Results.gabor; %

tsamples  = p.Results.tsamples;
timesd    = p.Results.timesd;
sceneFOV  = p.Results.sceneFOV;
distance  = p.Results.distance;

% Currently not being used because bgColor doesn't exist in sceneSet!
bgColor   = p.Results.bgColor;

%% Build Gaussian time series (soft window for stimulus in time)

tseries = exp(-(tsamples/timesd).^2);
tseries = ieScale(tseries,0,1);

%%  Scene parameters in general

% The size of the integration for the cone mosaic is 5 minutes for a line.
% There are two lines, so the stimulus should be more than 10 minutes.  We
% set up the cone mosaic to be about 15 minutes, and we will try shrinking
% it later.  And enlarging it.
%
% We create the scene and thus the oi to be a bit bigger to allow for eye
% movements.  We set the scene to be .35*60 ~ 21 minutes.  The oi is even a
% little bigger to allow for blur.
sparams.fov      = sceneFOV;    % degrees
sparams.distance = distance;    % Meters

% [NOT USED BY OISCREATE]
% sparams.bgColor  = bgColor;     % scaled from 0 to 1 

% Basic vernier parameters for the oiSequence.  Reverse order forces the
% allocation first so the array does not grow over the loop.
for ii = 3:-1:1
    ogparams(ii) = oGabor;
end

% Uniform field, no line, just the background color
ogparams(1).name     = 'uniform'; 
ogparams(1).contrast = 0;

% Left oriented Gabor on a zero background
ogparams(2).name     = 'leftOG';  

% Right oriented Gabor on a zero background
ogparams(3).name     = 'rightOG'; 
ogparams(3).ang      = -oGabor.ang;

% Put test params and scene params into P for use with oisCreate
P.sampleTimes       = tsamples;
P.testParameters    = ogparams([1 2]);
P.sceneParameters   = sparams;

[OG, scenes] = oisCreate('harmonic','blend', tseries, P);
% OG.visualize;
% ieAddObject(OG.oiFixed); ieAddObject(OG.oiModulated); oiWindow;
% ieAddObject(scenes{2}); sceneWindow;

% Gabor with opposite rotation
P.testParameters = ogparams([1 3]);
OG(2) = oisCreate('harmonic','blend', tseries, P);
% OG(2).visualize;
% ieAddObject(OG(2).oiFixed); ieAddObject(OG(2).oiModulated); oiWindow;


%%
save(fname,'OG','scenes','tseries','P');

% Print out the offset in degrees of arc sec 
% offsetDeg = sceneGet(scenes{1},'degrees per sample')*vparams(2).offset;
% fprintf('Offset in arc secs %.2f\n',offsetDeg*3600);

end
