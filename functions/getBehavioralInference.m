function P = getBehavioralInference(absorptions, currents, expParams, saveFolder)
% Classify oriented gabors simulated for differnece contrast levels and
% experimental manipulations

p = inputParser;
p.addRequired('expParams', @struct);
p.addRequired('saveFolder', @isstring);
p.addParameter('seed', 1, @(x) (isstring(x) | isscalar(x)));
p.addParameter('currentFlag', false, @islogical);
p.parse(expName, saveFolder, varargin{:});

% Make sure to reset the seed
rng(expParams.seed)

% Compute accuracy on fft component of cone absorptions
expParams.fftFlag        = true;

% Predefine matrix for predictions
nrContrasts      = length(expParams.contrastLevels);
nrEyemovTypes    = size(expParams.eyemovement,2);
nrEccen          = length(expParams.eccentricities);
nrSpatFreq       = length(expParams.spatFreq);
nrDefocusLevels  = length(expParams.defocusLevels);

% Preallocate space for percent correct
P = nan(nrContrasts,1);

% define pth to save classification results
classifierSavePth = fullfile(ogRootPath, 'data', 'classification', expName, saveFolder);
if ~exist('savePth', 'dir'); mkdir(classifierSavePth); end



for eccen = 1:nrEccen
    for df = 1:nrDefocusLevels
        for em = 1:max(nrEyemovTypes)
            for sf = expParams.spatFreq
                for c = 1:nrContrasts
                    
                    %% Load dataset
                    if expParams.verbose; fprintf('(%s): Loading file %s\n', mfilename, fname); end
 
                    fname = sprintf(...
                        'OGconeOutputs_contrast%1.3f_pa%d_eye%s_eccen%1.2f_defocus%1.2f_noise-random_sf%1.2f.mat',...
                        expParams.contrastLevels(c),expParams.polarAngle,sprintf('%i',expParams.eyemovement(:,em)), expParams.eccentricities(eccen), expParams.defocusLevels(df), sf);
                    
                    if expParams.currentFlag
                        fname = ['current_' fname];
                    end
                    
                    pth = fullfile(ogRootPath, 'data', expName, saveFolder, fname);
                    if ~exist(pth, 'file'), error('The file %s is not found', fname); end
                    
                    tmp = load(pth);
                    
                    % Truncate time dimension if necessary
                    if expParams.currentFlag
                        data = getfield(tmp,'current');
                    else
                        data = getfield(tmp,'absorptions');
                        if size(data,4) > 28
                            data = data(:,:,:,1:28,:); % truncate time samples (only include stimulus on period)
                        end
                    end
                    
                    %% Permute and 2D FFT data
                    if expParams.verbose;  fprintf('(%s): Classifying data\n', mfilename); end
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
                    if expParams.fftFlag; data  = abs(fft2(data)); end
                    
                    % reshape to all trials x [rows x colums x time] for classification
                    data = permute(data, [3 1 2 4]);
                    data = reshape(data, nTrials*2, []);
                    
                    % permute the trial order within each of the two classes
                    idx = [randperm(nTrials) randperm(nTrials)+nTrials];
                    
                    data = data(idx, :);
                    
                    label = [ones(nTrials, 1); -ones(nTrials, 1)];
                    
                    % Fit the SVM model.
                    cvmdl = fitcsvm(data, label, 'Standardize', true, 'KernelFunction', 'linear', 'kFold', 10);
                    
                    % cvmdl = crossval(mdl);
                    
                    % predict the data not in the training set.
                    classLoss = kfoldLoss(cvmdl);
                    
                    % Different type of linear classifier (faster, but less
                    % accurate)
                    %             mdl = fitclinear(data', label,  'KFold', 10, 'ObservationsIn', 'columns');
                    %             classLoss = kfoldLoss(mdl);
                    
                    P(c) = (1-classLoss) * 100;
                
                if expParams.verbose; disp(P); end
                
                end % contrasts
                
                if expParams.verbose; fprintf('(%s): Saving file\n', mfilename); end
                % Save classifier accuracy
                fname = sprintf(...
                    'Classify_coneOutputs_contrast%1.3f_pa%d_eye%s_eccen%1.2f_defocus%1.2f_noise-random_sf%1.2f',...
                    expParams.contrastLevels(c), expParams.polarAngle,sprintf('%i',expParams.eyemovement(:,em)), expParams.eccentricities(eccen), expParams.defocusLevels(df), sf);
                if expParams.currentFlag; fname = ['current_' fname]; end
                parsave(fullfile(classifierSavePth, sprintf('%s.mat', fname)),'P',P)
                
                
                
            end % sf
        end % eyemov
    end % defocus
end % eccen


% Visualize results:

% Visualize beta's
% betas(df, :,:,:) = reshape(cvmdl.Trained{2}.Beta, [nrows, ncols, tSamples]);
% mn_betas = squeeze(mean(betas(df,:,:,:),4));
% figure; imagesc(fftshift(mn_betas));
% set(gca,'XTick', ncols*[0.25 .5, 0.75],'XTickLabel',ncols*[-0.5, 0, 0.5], 'TickDir','out')
% set(gca,'YTick', ncols*[0.25 .5, 0.75],'YTickLabel',ncols*[-0.5, 0, 0.5], 'TickDir','out')
% box off; colormap gray; axis image; set(gca,'FontSize',20); colorbar;
% xlabel('Frequency (cycles/pixel)'); ylabel('Frequency (cycles/pixel)'); set(gca,'FontSize',20)
%                     
% title(sprintf('Condition %s - FFT at input freq: %1.3f x10^6', (expParams.defocusLevels{df}), mn_betas(8,3)*10^6));
% set(gca,'CLim', 4*10^-5*[-1 1]);


% Visualize percent correct

% figure; clf; set(gcf,'Color','w'); hold all;
% set(gca, 'XScale','log', 'XLim', [.005 max(expParams.contrastLevels)], 'XTick', [1:7, 10:10:100]/100, ...
%     'YLim', [40 100], 'TickDir','out','TickLength',[.015 .015]);
% ylabel('Classifier Accuracy')
% xlabel('Contrast level (Michelson)')

% plot(expParams.contrastLevels, P,'o-', 'LineWidth',2); drawnow;

% Save figure?
% savefig(fullfile(classifierSavePth, sprintf('%s.fig', fname)))
% hgexport(gcf,fullfile(classifierSavePth, sprintf('%s.eps', fname)))







return




