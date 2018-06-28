
cMosaic     = coneMosaic('center', [0.0015  0]);
cMosaic.integrationTime = 0.002;
cMosaic.setSizeToFOV(2);
rows = cMosaic.rows;
oneCycle = rows / 12;

x= (1:500)/500;
H = sin(x*2*pi*12);
G = exp(-(x-mean(x)).^2 / (2* (0.25).^2));
stim = G .* H;


n = 109;
t = (0:n-1)*cMosaic.integrationTime;

em = emCreate;


m =50;
yl = oneCycle*[-0.5 0.5];
xl = [min(t) max(t)];

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


%%
% Tremor + Drift
figure(4); clf
em.emFlag   = [1 1];
[emPaths, fixobj]     = cMosaic.emGenSequence(109*2, 'nTrials', 5, 'em', fixationalEM,'microsaccadeType', 'stats based');
emPaths     = emPaths(:,110:end,:);
subplot(121)
imagesc(xl, .5*[-rows rows], stim'); colormap gray; hold on; axis([xl yl]);
plot(t, emPaths(:,:,1)', 'LineWidth',2); title('MS + Drift in X pos', 'FontSize',16); axis([xl yl])
set(gca, 'FontSize',16);
xlabel('Time (s)', 'FontSize',16); ylabel('Cones', 'FontSize',16); box off;

subplot(122)
imagesc(xl, .5*[-rows rows], stim'); colormap gray; hold on;
plot(t, emPaths(:,:,2)', 'LineWidth',2); title('MS + Drift in Y pos', 'FontSize',16); axis([xl yl]);
set(gca, 'FontSize',16);
xlabel('Time (s)', 'FontSize',16); ylabel('Cones', 'FontSize',16); box off;

% print('~/Desktop/eyemovement_panelC','-depsc')

figure(5); clf;
% for ii = 1:50;  clf;
plot(squeeze(emPaths(5,:,1)), squeeze(emPaths(5,:,2))', 'LineWidth',2);  

xlim([-6 6]); ylim([-6 6]);
set(gca, 'FontSize',16, 'XGrid', 'on', 'YGrid', 'on');
xlabel('X Position (cones)', 'FontSize',16); ylabel('Y Position (cones)', 'FontSize',16); box off; axis square;

%  title(ii); waitforbuttonpress;
% end

title('Drift + MS in XY pos', 'FontSize',16);

% print('~/Desktop/eyemovement_panelD','-depsc')

