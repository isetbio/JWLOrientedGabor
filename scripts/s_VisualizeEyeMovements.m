%% s_VisualizeEyeMovements.m

% 1. Make CONE MOSAIC for a given eccentricity and polar angle
whichEye = 'left';
deg2m   = 0.3 * 0.001;          % (we first choose 3 deg per mm, .001 mm per meter, but now adjusted to the default in isetbio)
    
% Specify retinal location where stimulus is presented
cparams.eccentricity      = 4.5;             % Visual angle of stimulus center, in deg
cparams.polarAngle        = 90; % Polar angle (radians): 0 is right, pi/2 is superior, pi is left, 3*pi/2 inferior

% Cone mosaic field of view in degrees
cparams.cmFOV        = 2; % degrees

% Compute x,y position in m of center of retinal patch from ecc and angle
[x, y] = pol2cart(cparams.polarAngle, cparams.eccentricity);
x = x * deg2m;  y = y * deg2m;

% Create coneMosaic for particular location and eye
cMosaic = coneMosaic('center', [x, y], 'whichEye', whichEye); % cMosaic     = coneMosaic('center', [0.0015  0]);

% Set the field of view (degrees)
cMosaic.setSizeToFOV(cparams.cmFOV);

% Set the time step
cMosaic.integrationTime = 0.002;

% Get one cycle of a 1D gabor stimulus
rows = cMosaic.rows;
oneCycle = rows / 12;

x= (1:500)/500;
H = sin(x*2*pi*12);
G = exp(-(x-mean(x)).^2 / (2* (0.25).^2));
stim = G .* H;

% Create a time vector
n = 28; % number of timepoints
t = (0:n-1)*cMosaic.integrationTime;

%
% m =50;
yl = oneCycle*[-0.5 0.5];



%%
% Tremor + Drift

% Genereate eyemovements
% em.emFlag   = [1 1];

[emPaths, fixobj]     = cMosaic.emGenSequence(109*2, 'nTrials', 5, 'microsaccadeType', 'stats based', 'rSeed',6);

emPaths = emPaths(:,110:137,:);
t       = t(1:28); 
xl      = [min(t), max(t)];

figure(4); clf
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
% print('~/Desktop/eyemovement_panelC','-depsc')

% figure(5); clf;
% plot(squeeze(emPaths(5,:,1)), squeeze(emPaths(5,:,2))', 'LineWidth',2);  
% xlim([-6 6]); ylim([-6 6]);
% set(gca, 'FontSize',16, 'XGrid', 'on', 'YGrid', 'on');
% xlabel('X Position (cones)', 'FontSize',16); ylabel('Y Position (cones)', 'FontSize',16); box off; axis square;
% title('Drift + MS in XY pos', 'FontSize',16);

% print('~/Desktop/eyemovement_panelD','-depsc')

%% Obsolete
% figure(3);
% imagesc(xl, .5*[-rows rows], stim'); colormap gray; hold on; axis square; axis off;
% print('~/Desktop/eyemovement_panelA','-depsc')

% 
% yl = 2*oneCycle*[-.5 .5];
% xl = [min(t) max(t)];

% figure(6); imagesc(xl, .5*[-rows rows], stim'); axis([xl yl]); colormap gray; hold on; axis off;

% print('~/Desktop/eyemovement_panelB','-depsc')

% % Drift
% em.emFlag   = [0 1 0];
% emPaths     = cMosaic.emGenSequence(n*2, 'nTrials', m,'em', em);
% emPaths     = emPaths(:,n+1:end,:);
% subplot(3,2,1)
% imagesc(xl, .5*[-rows rows], stim'); colormap gray; hold on;
% plot(t,emPaths(:,:,1)', 'LineWidth',2);title('Drift in X pos', 'FontSize',16); axis([xl yl]);
% set(gca, 'FontSize',16);
% xlabel('Time (s)', 'FontSize',16); ylabel('Cones', 'FontSize',16); box off;
% 
% subplot(3,2,2)
% imagesc(xl, .5*[-rows rows], stim'); colormap gray; hold on;
% plot(t,emPaths(:,:,2)', 'LineWidth',2); title('Drift in Y pos', 'FontSize',16); axis([xl yl])
% set(gca, 'FontSize',16);
% xlabel('Time (s)', 'FontSize',16); ylabel('Cones', 'FontSize',16); box off;

% % Tremor
% % em.emFlag   = [1 0 0];
% emPaths     = cMosaic.emGenSequence(n*2, 'nTrials', m);
% emPaths     = emPaths(:,n+1:end,:);
% subplot(3,2,3)
% imagesc(xl, .5*[-rows rows], stim'); colormap gray; hold on;
% plot(t,emPaths(:,:,1)', 'LineWidth',2); title('Tremor in X pos', 'FontSize',16); axis([xl yl])
% set(gca, 'FontSize',16);
% xlabel('Time (s)', 'FontSize',16); ylabel('Cones', 'FontSize',16); box off;
% 
% subplot(3,2,4)
% imagesc(xl, .5*[-rows rows], stim'); colormap gray; hold on;
% plot(t,emPaths(:,:,2)', 'LineWidth',2);title('Tremor in Y pos', 'FontSize',16); axis([xl yl])
% set(gca, 'FontSize',16);
% xlabel('Time (s)', 'FontSize',16); ylabel('Cones', 'FontSize',16); box off;

