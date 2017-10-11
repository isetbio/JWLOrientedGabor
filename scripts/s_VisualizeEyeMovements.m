
cMosaic     = coneMosaic('center', [0.0015  0]);
cMosaic.integrationTime = 0.0010;
cMosaic.setSizeToFOV(2);
rows = cMosaic.rows;
oneCycle = rows / 12;

x= (1:500)/500;
H = sin(x*2*pi*12);
G = exp(-(x-mean(x)).^2 / (2* (0.25).^2));
stim = G .* H;


n = 56;
t = (0:n-1)*cMosaic.integrationTime;

em = emCreate;


m =50;
yl = 4.5*oneCycle*[-.5 .5];
xl = [min(t) max(t)];

figure(3);
imagesc(xl, .5*[-rows rows], stim'); colormap gray; hold on; axis square; axis off;



yl = 1*oneCycle*[-.5 .5];
xl = [min(t) max(t)];

figure(6); imagesc(xl, .5*[-rows rows], stim'); axis([xl yl]); colormap gray; hold on; axis off;

figure(4); clf
% Drift
em.emFlag   = [0 1 0];
emPaths     = cMosaic.emGenSequence(n*2, 'nTrials', m,'em', em);
emPaths     = emPaths(:,n+1:end,:);
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
emPaths     = cMosaic.emGenSequence(n*2, 'nTrials', m,'em', em);
emPaths     = emPaths(:,n+1:end,:);
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
emPaths     = cMosaic.emGenSequence(n*2, 'nTrials', m,'em', em);
emPaths     = emPaths(:,n+1:end,:);
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

%%
figure(5);
% for ii = 1:50;  clf;
plot(squeeze(emPaths(30,:,1)), squeeze(emPaths(30,:,2))', 'LineWidth',2);  

xlim([-2 2]); ylim([-2 2]);
set(gca, 'FontSize',16, 'XGrid', 'on', 'YGrid', 'on');
xlabel('X Position (cones)', 'FontSize',16); ylabel('Y Position (cones)', 'FontSize',16); box off; axis square;

%  title(ii); waitforbuttonpress;
% end

title('Tremor + Drift in XY pos', 'FontSize',16);
