%% s_visualizeClassifierWeights

expName = 'defocus';
subFolderName_toSave = 'VisualizedClassifierWeights';
subFolderName_toLoad = 'run1';

% Load experiment parameters
expParams = loadExpParams(expName, false);

% Save figures?
saveFigures = false;

% Allow computation of accuracy for cone current
currentFlag    = false;

% Compute accuracy on fft component of cone absorptions
fftFlag        = true;

% Use default and 100% contrast stim
eccen   = 1; % 4.5 degreers
df      = 1; % 0 diopters
em      = 1; % Both drift and MS
sf      = expParams.spatFreq; % 4 cpd
c       = 1; % 100% contrast
                    
% Load dataset
fname = sprintf(...
    'OGconeOutputs_contrast%1.3f_pa%d_eye%s_eccen%1.2f_defocus%1.2f_noise-random_sf%1.2f.mat',...
    c, expParams.polarAngle, sprintf('%i',expParams.eyemovement(:,em)), expParams.eccentricities(eccen), expParams.defocusLevels(df), sf);

if currentFlag
    fname = ['current_' fname];
end

pth = fullfile(ogRootPath, 'data', expParams.name, subFolderName_toLoad, fname);
if ~exist(pth, 'file'), error('The file %s is not found', fname); end

tmp = load(pth);

if currentFlag
    data = getfield(tmp,'current');
else
    data = getfield(tmp,'absorptions');
    data = data(:,:,:,1:28,:); % Truncate time samples, where blank stimulus was presented.
end

fprintf('Loading and classifying %s\n', fname);
% Get the trials and samples (should be the data for all data sets though
nStimuli = size(data,5);
nTrials  = size(data,1) * nStimuli/2;
tSamples = size(data,4);
nrows    = size(data,2);
ncols    = size(data,3);
% absorptions is trials x rows x cols x time points x stimuli

%   permute to trials x stimuli x rows x cols x time points
data = permute(data, [1 5 2:4]);

%   reshape to (trials x stimuli) x rows x cols x time points
data = reshape(data, [], nrows, ncols, tSamples);

% permute to rows x cols x (trials x stimuli) x time points
data  = permute(data, [2 3 1 4]);

% Compute fourier transform the cone array outputs
if fftFlag; data  = abs(fft2(data)); end

% reshape to all trials x [rows x colums x time] for classification
data = permute(data, [3 1 2 4]);
data = reshape(data, nTrials*2, []);

% permute the trial order within each of the two classes
idx = [randperm(nTrials) randperm(nTrials)+nTrials];

data = data(idx, :);

label = [ones(nTrials, 1); -ones(nTrials, 1)];

% Fit the SVM model.
cvmdl = fitcsvm(data, label, 'Standardize', true, 'KernelFunction', 'linear', 'kFold', 10);

% predict the data not in the training set.
classLoss = kfoldLoss(cvmdl);

P(c) = (1-classLoss) * 100;

% visualize beta's
betas(em, :,:,:) = reshape(cvmdl.Trained{1}.Beta, [nrows, ncols, tSamples]);
mn_betas = squeeze(mean(betas(em,:,:,:),4));
imagesc(fftshift(mn_betas)); box off; set(gca, 'TickDir', 'out', 'FontSize',12); 
colormap gray; axis image; colorbar;

title(sprintf('FFT at input freq: %1.3f x10^6', mn_betas(8,3)*10^6));
% set(gca,'CLim', 4*10^-5*[-1 1]);

if saveFigures
    savePth = fullfile(ogRootPath, 'figs', 'subFolderName_toSave', expName);
    if ~exist('savePth', 'dir'); mkdir(savePth); end;
    print(fullfile(savePth, 'classifierWeights_averagedAcrossTime'),'-depsc')
end
                    

