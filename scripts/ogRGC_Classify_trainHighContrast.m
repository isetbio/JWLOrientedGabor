function ogRGC_Classify_trainHighContrast(expName)

% Script with first attempt to classify oriented gabors simulated at
% 7 different contrast levels, 4 polar angles.

%% Classify

% Load experiment parameters
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

P = nan(nrContrasts,1);

savePth = fullfile(ogRootPath, 'data', 'classification', expName, subFolderName_toSave);
if ~exist('savePth', 'dir'); mkdir(savePth); end;

% Init figure
figure; clf; set(gcf,'Color','w'); hold all;
set(gca, 'XScale','log', 'XLim', [.005 max(expParams.contrastLevels)], 'XTick', [1:7, 10:10:100]/100, ...
    'YLim', [40 100], 'TickDir','out','TickLength',[.015 .015]);
ylabel('Classifier Accuracy')
xlabel('Contrast level (Michelson)')


for eccen = 1:nrEccen
    for df = 1:nrDefocusLevels
        for em = 1:max(nrEyemovTypes)
            for sf = expParams.spatFreq
                
                %% Train on high contrast
                [data, nTrials] = loadAndPermuteData(expParams, nrContrasts(end), em, eccen, df, sf, currentFlag, subFolderName_toLoad);
                
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
                
                %% Test on all contrast
                for c = 1:nrContrasts
                    
                    [data, nTrials] = loadAndPermuteData(expParams, c, em, eccen, df, sf, currentFlag);
                    
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
                    
 
                    P(c==expParams.contrastLevels) = (sum(label==predictedLabel)/length(label))*100;                                        
                    
                end
                
                disp(P);
                
                % Save classifier accuracy
                fname = sprintf(...
                    'Classify_coneOutputs_contrast%1.3f_pa%d_eye%s_eccen%1.2f_defocus%1.2f_noise-random_sf%1.2f',...
                    expParams.contrastLevels(c), expParams.polarAngle,sprintf('%i',expParams.eyemovement(:,em)), expParams.eccentricities(eccen), expParams.defocusLevels(df), sf);
                if currentFlag; fname = ['current_' fname]; end
                parsave(fullfile(savePth, sprintf('%s.mat', fname)),'P',P)

            end
        end
    end
end





return




