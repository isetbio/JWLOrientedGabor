%% s_ogRGC_Classify

% Script with first attempt to classify oriented gabors simulated at
% 7 different contrast levels, 4 polar angles.

FFTflag = false;
%% Classify

contrastLevels = [0.01:0.01:0.09, 0.1:0.1:1.0];
polarAngles    = 0; % [0 90 180 270];
eyemovement    = {'110'};%{'000', '100', '010', '001'};

P = nan(length(polarAngles),length(contrastLevels),length(eyemovement));
% svmMdl = cell(1, length(contrastLevels));

for pa = polarAngles
    for c = contrastLevels
        for em = 1:length(eyemovement)
            % Load dataset
            load(fullfile(ogRootPath, 'data', sprintf('OGconeOutputs_contrast%1.2f_pa%d_eye%s.mat',c,pa,cell2mat(eyemovement(em)))));
            
            % Get the trials and samples (should be the data for all data sets though
            nTrials = size(absorptions.cw,1);
            tSamples = size(absorptions.cw,4);
            
            % If requested, fourier transform the cone array outputs
            if FFTflag
                absorptions.cwF  = abs(fft2(permute(absorptions.cw, [2 3 1 4])));
                absorptions.ccwF = abs(fft2(permute(absorptions.ccw, [2 3 1 4])));
                
                imgListCW  = trial2Matrix(permute(absorptions.cwF, [3 1 2 4]));
                imgListCCW = trial2Matrix(permute(absorptions.ccwF, [3 1 2 4]));
                
            else
                % Reformat the time series for the PCA analysis
                %
                % imgListX matrix contains the temporal response for a pixel in a
                % column. The rows represent time samples times number of trials.
                % These are the temporal responses across all trials and time
                % points. The columns represent the cells.
                % 4D array input
                imgListCW  = trial2Matrix(absorptions.cw);
                imgListCCW = trial2Matrix(absorptions.ccw);
            end
            
            % Concatenate the matrices of the two stimuli
            imgList = cat(1,imgListCW,imgListCCW);
            
            % compute the imagebases of the two stimuli
            imageBasis = ogPCA(cat(1,absorptions.cw,absorptions.ccw));
            
            % Time series of weights
            weightSeries  = imgList * imageBasis;
            
            %% Start classification training
            %
            % Put the weights from each trial into the rows of a matrix
            % Each row is another trial
            nWeights = size(weightSeries,2);
            data = zeros(2*nTrials,nWeights*tSamples);
            for ii = 1 : (2*nTrials)
                start = (ii-1)*tSamples + 1;
                thisTrial = weightSeries(start:(start+tSamples - 1),:);
                data(ii,:) = thisTrial(:)';
            end
            label = [ones(nTrials, 1); -ones(nTrials, 1)];
            
            % Select some of the data (80%) as the training set.
            train_index = zeros(nTrials, 1);
            train_index(randperm(nTrials, round(0.8*nTrials))) = 1;
            train_index = train_index > 0;
            
            % The cw and ccw trials are still matched ????
            train_index = repmat(train_index, 2, 1);
            
            % Fit the SVM model.
            mdl = fitcsvm(data(train_index, :), label(train_index), ...
                'KernelFunction', 'linear');
            
            % predict the data not in the training set.
            yp = predict(mdl, data(~train_index, :));
            classLoss = sum(label(~train_index) ~= yp) / length(yp);
            
            P(pa==polarAngles,c==contrastLevels,em) = (1-classLoss) * 100;
            
            
        end
    end
end


 disp(P);

% Visualize
labels = {'Polar Angle: 0'};%,'Polar Angle: 90','Polar Angle: 180','Polar Angle: 270'};

colors = lines(length(eyemovement));
figure(1); clf; set(gcf,'Color','w'); hold all;
plot(contrastLevels, squeeze(P),'Color', 'k', 'LineWidth',2);
set(gca, 'XScale','log', 'YLim', [0 100], 'TickDir','out','TickLength',[.015 .015]);
ylabel('Classifier Accuracy')
xlabel('Contrast level (Michelson)')

savefig('~/matlab/git/toolboxes/JWLOrientedGabor/figs/drift_contrastVSperformance_noFFT.fig')
hgexport(gcf,'~/matlab/git/toolboxes/JWLOrientedGabor/figs/drift_contrastVSperformance_noFFT.eps')


%%
% plot(contrastLevels,P(1,:),'Color', colors(1,:), 'LineWidth',2);
% plot(contrastLevels,P(2,:),'Color', colors(2,:), 'LineWidth',2);
% plot(contrastLevels,P(3,:),'Color', colors(3,:), 'LineWidth',2);
% plot(contrastLevels,P(4,:),'Color', colors(4,:), 'LineWidth',2);
% legend(labels);
% box off;
% xlabel('Contrast level (Michelson)');
% ylabel('Classifier Accuracy')
% set(gca, 'XLim', [0.4 1], 'YLim', [0 100],'TickDir','out','TickLength',[.015 .015]);

return



