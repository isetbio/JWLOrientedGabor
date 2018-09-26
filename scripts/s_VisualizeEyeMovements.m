%% s_VisualizeEyemovements.m

% Script to visualize eye movement paths and relative size compared to one
% cycle of a 1D Gabor stimulus.

%% 0. Define parameters

% Unit converter
deg2m   = 0.3 * 0.001;          % (we first choose 3 deg per mm, .001 mm per meter, but now adjusted to the default in isetbio)

% What eye to make the cone mosaic of
whichEye = 'left';
    
% Specify retinal location where stimulus is presented
cparams.eccentricity      = 4.5;             % Visual angle of stimulus center, in deg
cparams.polarAngle        = 90; % Polar angle (radians): 0 is right, pi/2 is superior, pi is left, 3*pi/2 inferior

% Cone mosaic field of view in degrees
cparams.cmFOV        = 2; % degrees


%% 1. Make cone mosaic for a given eccentricity and polar angle

% Compute x,y position in m of center of retinal patch from ecc and angle
[x, y] = pol2cart(cparams.polarAngle, cparams.eccentricity);
x = x * deg2m;  y = y * deg2m;

% Create coneMosaic for particular location and eye
cMosaic = coneMosaic('center', [x, y], 'whichEye', whichEye); % cMosaic     = coneMosaic('center', [0.0015  0]);

% Set the field of view (degrees)
cMosaic.setSizeToFOV(cparams.cmFOV);



%% 2. Create one cycle of 1D Gabor stimulus

% Get cone mosaic size
rows = cMosaic.rows;

% Define how many cones one cycle is
oneCycle = rows / 12;

% Define the Gabor stimulus
x = (1:500)/500;    % Space sampling
H = sin(x*2*pi*12); % Harmonic
G = exp(-(x-mean(x)).^2 / (2* (0.25).^2)); % Gaussian window
stim = G .* H;      % Gabor stimulus = Gaussian .* harmonic

% Define y axis limits
yl = oneCycle*[-0.5 0.5];

%% 3. Create a time vector

% Set the number of time steps and sampling rate
n = 28;
cMosaic.integrationTime = 0.002;

% Get time vector
t = (0:n-1)*cMosaic.integrationTime;

% Define x axis limits
xl      = [min(t), max(t)];


%% 4. Create eye movement paths with tremor + drift

% Genereate eyemovements
[emPaths, fixobj]     = cMosaic.emGenSequence(109*2, 'nTrials', 5, 'microsaccadeType', 'stats based', 'rSeed',6);

% Truncate paths to match stimulus duration (54 ms)
emPaths = emPaths(:,110:137,:); 

%% 5. Plot Figure 6A

figure(1); clf; set(gcf, 'NumberTitle', 'Off', 'Name', 'Figure 6A - Effect of Eye movements');
subplot(121)
imagesc(xl, .5*[-rows rows], stim'); colormap gray; hold on; axis([xl yl]);
plot(t, emPaths(:,:,1)', 'LineWidth',4); title('MS + Drift in X pos', 'FontSize',16); axis([xl yl])
set(gca, 'FontSize',16, 'TickDir', 'out');
xlabel('Time (s)', 'FontSize',16); ylabel('Position (Cones)', 'FontSize',16); box off;
axis square

subplot(122)
imagesc(xl, .5*[-rows rows], stim'); colormap gray; hold on;
plot(t, emPaths(:,:,2)', 'LineWidth',4); title('MS + Drift in Y pos', 'FontSize',16); axis([xl yl]);
set(gca, 'FontSize',16, 'TickDir', 'out');
xlabel('Time (s)', 'FontSize',16); ylabel('Position (Cones)', 'FontSize',16); box off;
axis square
