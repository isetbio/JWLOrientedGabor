%% s_ogRGC_Classify_trainHighContrast

% Script with first attempt to classify oriented gabors simulated at
% 7 different contrast levels, 4 polar angles.

%% Classify

% Load experiment parameters
expName = 'defocus';
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

savePth = fullfile(ogRootPath, 'data', 'classification', 'HPC', expName, '100trials_trainHighContrast');
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
                
                % Train on high contrast
                [data, nTrials] = loadAndPermuteData(expParams, nrContrasts(end), em, eccen, df, sf, currentFlag);
                
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

                    P(c) = (sum(label==predictedLabel)/length(label))*100;
                    
                    
                    % visualize beta's
                    %                     betas(df, :,:,:) = reshape(cvmdl.Trained{2}.Beta, [nrows, ncols, tSamples]);
                    %                     mn_betas = squeeze(mean(betas(df,:,:,:),4));
                    %                     figure; imagesc(fftshift(mn_betas));
                    %                     set(gca,'XTick', ncols*[0.25 .5, 0.75],'XTickLabel',ncols*[-0.5, 0, 0.5], 'TickDir','out')
                    %                     set(gca,'YTick', ncols*[0.25 .5, 0.75],'YTickLabel',ncols*[-0.5, 0, 0.5], 'TickDir','out')
                    %                     box off; colormap gray; axis image; set(gca,'FontSize',20); colorbar;
                    %                     xlabel('Frequency (cycles/pixel)'); ylabel('Frequency (cycles/pixel)'); set(gca,'FontSize',20)
                    %
                    %                     title(sprintf('Condition %s - FFT at input freq: %1.3f x10^6', (expParams.defocusLevels{df}), mn_betas(8,3)*10^6));
                    %                     %                     set(gca,'CLim', 4*10^-5*[-1 1]);
                    
                    
                    
                    
                end
                
                disp(P);
                
                % Save classifier accuracy
                fname = sprintf(...
                    'Classify_coneOutputs_contrast%1.3f_pa%d_eye%s_eccen%1.2f_defocus%1.2f_noise-random_sf%1.2f_trainHighContrast',...
                    expParams.contrastLevels(c), expParams.polarAngle,sprintf('%i',expParams.eyemovement(:,em)), expParams.eccentricities(eccen), expParams.defocusLevels(df), sf);
                if currentFlag; fname = ['current_' fname]; end
                parsave(fullfile(savePth, sprintf('%s.mat', fname)),'P',P)
                
                
                % Visualize
                plot(expParams.contrastLevels, P,'o-', 'LineWidth',2); drawnow;
            end
        end
    end
end




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








return




