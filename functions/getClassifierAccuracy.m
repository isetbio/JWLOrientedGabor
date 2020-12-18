function P = getClassifierAccuracy(data)
% Function to train and test linear SVM classifier on cone data with cross-validation:
%       P = getClassifierAccuracy(data)
%
% INPUTS: 
%   data        : 5 dimensional array (with trials x rows x cols x time
%                                       samples x stimuli) 
% OUTPUTS:
%   P           : classifier accuracy of computational observer model in
%                   percent correct for given absorption dataset


% Get dimensions of data
fprintf('(%s): Loading and classifying\n', mfilename);
% Get the trials and samples (should be the data for all data sets though
nStimuli = size(data,5);
nTrials  = size(data,1) * nStimuli/2;
tSamples = size(data,4);
nrows    = size(data,2);
ncols    = size(data,3);

%   permute to trials x stimuli x rows x cols x time points
data = permute(data, [1 5 2:4]);

%   reshape to (trials x stimuli) x rows x cols x time points
data = reshape(data, [], nrows, ncols, tSamples);

% permute to rows x cols x (trials x stimuli) x time points
data  = permute(data, [2 3 1 4]);

% Compute fourier transform the cone array outputs
data  = abs(fft2(data));

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

% Get percent accuracy
P = (1-classLoss) * 100;

% visualize beta's
% betas = reshape(cvmdl.Trained{1}.Beta, [nrows, ncols, tSamples]);
% mn_betas = squeeze(mean(betas,3));
% imagesc(fftshift(mn_betas)); box off; set(gca, 'TickDir', 'out', 'FontSize',12); 
% colormap gray; axis image; colorbar;
% 
% title(sprintf('FFT at input freq: %1.3f x10^6', mn_betas(8,3)*10^6));
% set(gca,'CLim', 1.0e-03 .*[-0.3624,0.3624]);
%     savePth = fullfile(pfRV1rootPath, 'figures');
%     print(fullfile(savePth, 'classifierWeights_averagedAcrossTime.eps'),'-depsc')

       

return