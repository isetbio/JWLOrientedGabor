
cMosaic     = coneMosaic('center', [0.000002  0]);
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
yl = 12*oneCycle*[-.5 .5];
xl = [min(t) max(t)];

% Drift
em.emFlag   = [0 1 0];
emPaths     = cMosaic.emGenSequence(n, 'nTrials', m,'em', em);
subplot(3,2,1)
imagesc(xl, .5*[-rows rows], stim'); colormap gray; hold on;
plot(t,emPaths(:,:,1)');title('X-Drift'); axis([xl yl])
xlabel('Time (s)'); ylabel('Cones')

subplot(3,2,2)
imagesc(xl, .5*[-rows rows], stim'); colormap gray; hold on;
plot(t,emPaths(:,:,2)');title('Y-Drift'); axis([xl yl])
xlabel('Time (s)'); ylabel('Cones')

% Tremor
em.emFlag   = [1 0 0];
emPaths     = cMosaic.emGenSequence(n, 'nTrials', m,'em', em);
subplot(3,2,3)
imagesc(xl, .5*[-rows rows], stim'); colormap gray; hold on;
plot(t,emPaths(:,:,1)');title('X-Tremor'); axis([xl yl])
xlabel('Time (s)'); ylabel('Cones')

subplot(3,2,4)
imagesc(xl, .5*[-rows rows], stim'); colormap gray; hold on;
plot(t,emPaths(:,:,2)');title('Y-Tremor'); axis([xl yl])
xlabel('Time (s)'); ylabel('Cones')



% Tremor + Drift
em.emFlag   = [1 1 0];
emPaths     = cMosaic.emGenSequence(n, 'nTrials', m,'em', em);
subplot(3,2,5)
imagesc(xl, .5*[-rows rows], stim'); colormap gray; hold on;
plot(t, emPaths(:,:,1)');title('X-Tremor + Drift'); axis([xl yl])
xlabel('Time (s)'); ylabel('Cones')

subplot(3,2,6)
imagesc(xl, .5*[-rows rows], stim'); colormap gray; hold on;
plot(t, emPaths(:,:,2)');title('Y-Tremor + Drift'); axis([xl yl])
xlabel('Time (s)'); ylabel('Cones')

