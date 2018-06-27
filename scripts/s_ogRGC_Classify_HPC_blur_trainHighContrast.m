%% s_ogRGC_Classify_HPC_blur

% Script with first attempt to classify oriented gabors simulated at
% 7 different contrast levels, 4 polar angles.

%% Classify

% Load experiment parameters
expName = 'defocus';
subFolderName_toSave = '100trials_trainHighContrast';
subFolderName_toLoad = 'paddedStim';
expParams = loadExpParams(expName, false);

% Compute accuracy for cone current as well
currentFlag    = false;

% Compute accuracy on fft component of cone absorptions
fftFlag        = true;

% Predefine matrix for predictions
nrContrasts      = length(expParams.contrastLevels);
nrEyemovTypes    = size(expParams.eyemovement,2);
nrEccen          = length(expParams.eccentricities);
nrSpatFreq       = length(expParams.spatFreq);
nrDefocusLevels  = length(expParams.defocusLevels);


savePth = fullfile(ogRootPath, 'data', 'classification', expName, subFolderName_toSave);
if ~exist('savePth', 'dir'); mkdir(savePth); end;

% Init figure
figure; clf; set(gcf,'Color','w'); hold all;
set(gca, 'XScale','log', 'XLim', [.005 max(expParams.contrastLevels)], 'XTick', [1:7, 10:10:100]/100, ...
    'YLim', [40 100], 'TickDir','out','TickLength',[.015 .015]);
ylabel('Classifier Accuracy')
xlabel('Contrast level (Michelson)')

eccen = 1;



em = 1;
sf = expParams.spatFreq;

for df = 1:nrDefocusLevels
    P = nan(nrContrasts,1);
    
    % Train on high contrast
    [data, fname] = loadAndPermuteData(expParams, nrContrasts(end), em, eccen, df, sf, currentFlag, subFolderName_toLoad);
    nTrials = size(data,
    
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
    mdl = fitcsvm(data, label, 'Standardize', true, 'KernelFunction', 'linear');
    
    for c = 1:nrContrasts
        
        % data array = trials x rows x cols x time points x stimuli
        fprintf('Loading and classifying %s\n', fname);
        [data, fname] = loadAndPermuteData(expParams, c, em, eccen, df, sf, currentFlag, subFolderName_toLoad);
        
        % Get nr of trials (/2 for dividing the two phases)
        nTrials  = size(data,1) * nStimuli/2;
        
        % Compute fourier transform the cone array outputs
        if fftFlag; data  = abs(fft2(data)); end
        
        % reshape to all trials x [rows x colums x time] for classification
        data = permute(data, [3 1 2 4]);
        data = reshape(data, nTrials*2, []);
        
        % permute the trial order within each of the two classes
        idx = [randperm(nTrials) randperm(nTrials)+nTrials];
        
        data = data(idx, :);
        
        label = [ones(nTrials, 1); -ones(nTrials, 1)];
        
        [predictedLabel,score] = predict(mdl,data);
        
        % Fit the SVM model.
        %                     cvmdl = fitcsvm(data, label, 'Standardize', true, 'KernelFunction', 'linear', 'kFold', 10);
        
        %             cvmdl = crossval(mdl);
        
        % predict the data not in the training set.
        %                     classLoss = kfoldLoss(cvmdl);
        
        % Different type of linear classifier (faster, but less
        % accurate)
        %             mdl = fitclinear(data', label,  'KFold', 10, 'ObservationsIn', 'columns');
        %             classLoss = kfoldLoss(mdl);
        
        %                     P(c) = (1-classLoss) * 100;
        
        P(c==expParams.contrastLevels) = (sum(label==predictedLabel)/length(label))*100;
        
    end
    
    disp(P);
    
    % Save classifier accuracy
    fname = sprintf(...
        'Classify_coneOutputs_contrast%1.3f_pa%d_eye%s_eccen%1.2f_defocus%1.2f_noise-random_sf%1.2f',...
        c, expParams.polarAngle,sprintf('%i',expParams.eyemovement(:,em)), expParams.eccentricities(eccen), expParams.defocusLevels(df), sf);
    if currentFlag; fname = ['current_' fname]; end
    parsave(fullfile(savePth, sprintf('%s.mat', fname)),'P',P)
    
    
    % Visualize
    %                 plot(expParams.contrastLevels, P,'o-', 'LineWidth',2); drawnow;
end



return


% Save figure?
% savefig(fullfile(ogRootPath, 'data', 'classification', sprintf('%s.fig', fname)))
% hgexport(gcf,fullfile(ogRootPath, 'data', 'classification', sprintf('%s.eps', fname)))


% %% visualize multiple classifier accuracy's
% plot(contrastLevels,P(:,:,1,4),'Color', colors(1,:), 'LineWidth',2);
% plot(contrastLevels,P(:,:,2,4),'Color', colors(2,:), 'LineWidth',2);
% plot(contrastLevels,P(:,:,3,4),'Color', colors(3,:), 'LineWidth',2);
% plot(contrastLevels,P(:,:,4,4),'Color', colors(4,:), 'LineWidth',2);
% plot(contrastLevels,P(:,:,5,4),'Color', colors(5,:), 'LineWidth',2);
% legend(eyemovement);
% box off;
% xlabel('Contrast level (Michelson)');
% ylabel('Classifier Accuracy')
% set(gca, 'XLim', [0.008 .6], 'YLim', [0 100],'TickDir','out','TickLength',[.015 .015]);













