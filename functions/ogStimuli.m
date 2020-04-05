function [OG, scenes, tseries, fname] = ogStimuli(varargin)
% OGSTIMULI - Create oriented stimuli (e.g., cw and ccw for 2 AFC expt)
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
%  OG       - 1x2 OIsequence array for ccw- and cw-oriented Gabors
%  scenes   - A 1x2 cell with a uniform field scene and oriented gabor scene
%  tseries  - A vector of the temporal window modulator for the scene
%  fname    - Path and filename where stimulus is stored
%
% Example: 
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
%   Luminance 2 log millilamberts ?? 
%   It appears that luminance is not used?? Check.
%
% JW & EK, NYU ISETBIO Team, 2017
%
% See also vaStimuli in WLVernier repository

%% Load a stored image, if it exists
% 
% pth = fullfile(ogRootPath, 'Stimulus');
% if ~exist(pth, 'dir'), mkdir(pth); end
% fname = fullfile(pth, 'stimulus.mat');
% 
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
p.addParameter('noStimPhase',false, @islogical);          % introduce stimulus phase shift or not

p.parse(varargin{:});

%% Derive variables from parsed inputs
oGabor     = p.Results.gabor; %

tsamples  = p.Results.tsamples;
timesd    = p.Results.timesd;
sceneFOV  = p.Results.sceneFOV;
distance  = p.Results.distance;
noStimPhase  = p.Results.noStimPhase;

% Currently not being used because bgColor doesn't exist in sceneSet!
bgColor   = p.Results.bgColor;

%% Build Gaussian time series (soft window for stimulus in time)

tseries = single(tsamples <= 0.054);

% tseries = exp(-(tsamples/timesd).^2);
% tseries = ieScale(tseries,0,1);

%%  Scene parameters in general
sparams.fov      = sceneFOV;    % degrees
sparams.distance = distance;    % Meters

% [NOT USED BY OISCREATE]
% sparams.bgColor  = bgColor;     % scaled from 0 to 1 

% Basic Gabor parameters for the oiSequence.  We make 3 for CW, CCW,
% and blank. CW and CCW are for the 2 AFC. Blank is for mixing with the
% Gabors for the temporal windowing.
ogparams(1:5) = oGabor;

% Uniform field, no Gabor, just the background color
ogparams(1).name     = 'uniform'; 
ogparams(1).contrast = 0;

% CCW oriented Gabor on a zero background. The function which makes the
% Gabors, harmonicP, treats 0 as vertical and clockwise as positive.
ogparams(2).name     = 'ccw_sin_OG';  
ogparams(2).ang      = -oGabor.ang;

ogparams(3).ang      = -oGabor.ang;

if noStimPhase == true
    ogparams(3).name     = 'ccw_sin_OG';  
else
    ogparams(3).name     = 'ccw_cos_OG';  
    ogparams(3).ph       =  oGabor.ph - pi;
end

% CW oriented Gabor on a zero background
ogparams(4).name     = 'cw_sin_OG'; 

if noStimPhase == true
    ogparams(5).name     = 'cw_sin_OG'; 
else
    ogparams(5).name     = 'cw_cos_OG'; 
    ogparams(5).ph       =  oGabor.ph - pi;
end


% Put test params and scene params into P for use with oisCreate
P.sampleTimes       = tsamples;
P.sceneParameters   = sparams;

if isfield(p.Unmatched,'oi')
    P.oi = p.Unmatched.oi;
end

% Blend uniform and ccw Gabor for temporal windowing
for ii = 1:4

    P.testParameters = ogparams([1 ii+1]);
    [OG(ii), scenes{ii}] = oisCreate('harmonic','blend', tseries, P);
    
    % OG.visualize;
    % ieAddObject(OG.oiFixed); ieAddObject(OG.oiModulated); oiWindow;
    % ieAddObject(scenes{2}); sceneWindow;

end

% save(fname,'OG','scenes','tseries','P');

% Print out the offset in degrees of arc sec 
% offsetDeg = sceneGet(scenes{1},'degrees per sample')*vparams(2).offset;
% fprintf('Offset in arc secs %.2f\n',offsetDeg*3600);

end
