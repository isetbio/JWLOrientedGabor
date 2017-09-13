
cMosaic     = coneMosaic('center', [0.002  0]);
cMosaic.integrationTime = .002;
cMosaic.setSizeToFOV(2);
rows = cMosaic.rows;
oneCycle = rows / 12;

x= (1:500)/500;
H = sin(x*2*pi*12);
G = exp(-(x-mean(x)).^2 / (2* (0.125).^2));
stim = G .* H;


n = 201;
t = (0:n-1)*cMosaic.integrationTime;

em = emCreate;

figure(4); clf
m =50;
yl = 1*oneCycle*[-.5 .5];
xl = [min(t) max(t)];

% Drift
em.emFlag   = [0 1 0];
emPaths     = cMosaic.emGenSequence(n, 'nTrials', m,'em', em);
subplot(3,2,1)
imagesc(xl, .5*[-rows rows], stim'); colormap gray; hold on;
plot(t,emPaths(:,:,1)', 'LineWidth',2);title('Drift in X pos', 'FontSize',16); axis([xl yl]);
set(gca, 'FontSize',16);
xlabel('Time (s)', 'FontSize',16); ylabel('Cones', 'FontSize',16); box off;

subplot(3,2,2)
imagesc(xl, .5*[-rows rows], stim'); colormap gray; hold on;
plot(t,emPaths(:,:,2)', 'LineWidth',2); title('Drift in Y pos', 'FontSize',16); axis([xl yl])
set(gca, 'FontSize',16);
xlabel('Time (s)', 'FontSize',16); ylabel('Cones', 'FontSize',16); box off;

% Tremor
em.emFlag   = [1 0 0];
emPaths     = cMosaic.emGenSequence(n, 'nTrials', m,'em', em);
subplot(3,2,3)
imagesc(xl, .5*[-rows rows], stim'); colormap gray; hold on;
plot(t,emPaths(:,:,1)', 'LineWidth',2); title('Tremor in X pos', 'FontSize',16); axis([xl yl])
set(gca, 'FontSize',16);
xlabel('Time (s)', 'FontSize',16); ylabel('Cones', 'FontSize',16); box off;

subplot(3,2,4)
imagesc(xl, .5*[-rows rows], stim'); colormap gray; hold on;
plot(t,emPaths(:,:,2)', 'LineWidth',2);title('Tremor in Y pos', 'FontSize',16); axis([xl yl])
set(gca, 'FontSize',16);
xlabel('Time (s)', 'FontSize',16); ylabel('Cones', 'FontSize',16); box off;



% Tremor + Drift
em.emFlag   = [1 1 0];
emPaths     = cMosaic.emGenSequence(n, 'nTrials', m,'em', em);
subplot(3,2,5)
imagesc(xl, .5*[-rows rows], stim'); colormap gray; hold on;
plot(t, emPaths(:,:,1)', 'LineWidth',2); title('Tremor + Drift in X pos', 'FontSize',16); axis([xl yl])
set(gca, 'FontSize',16);
xlabel('Time (s)', 'FontSize',16); ylabel('Cones', 'FontSize',16); box off;

subplot(3,2,6)
imagesc(xl, .5*[-rows rows], stim'); colormap gray; hold on;
plot(t, emPaths(:,:,2)', 'LineWidth',2); title('Tremor + Drift in Y pos', 'FontSize',16); axis([xl yl]);
set(gca, 'FontSize',16);
xlabel('Time (s)', 'FontSize',16); ylabel('Cones', 'FontSize',16); box off;


figure;
plot(squeeze(emPaths(1,:,1)), squeeze(emPaths(1,:,2))', 'LineWidth',2); title('Tremor + Drift in XY pos', 'FontSize',16); 

xlim([-2 2]); ylim([-2 2]);
set(gca, 'FontSize',16, 'XGrid', 'on', 'YGrid', 'on');
xlabel('X Position (cones)', 'FontSize',16); ylabel('Y Position (cones)', 'FontSize',16); box off; axis square;

