%% s_ogRGC_Classify

% Script with first attempt to classify oriented gabors simulated at
% 7 different contrast levels, 4 polar angles.


%% Classify

contrastLevels = [0.4:0.1:1.0];
polarAngles    = [0 90 180 270];

P = nan(length(polarAngles),length(contrastLevels));
% svmMdl = cell(1, length(contrastLevels));

for pa = polarAngles
    for c = contrastLevels
        
        % Load dataset
        load(fullfile(ogRootPath, 'data', sprintf('OGconeOutputs_contrast%1.1f_pa%d.mat',c,pa)));
        
        % Get the trials and samples (should be the data for all data sets though
        nTrials = size(absorptions.cw,1);
        tSamples = size(absorptions.cw,4);
        
        % Reformat the time series for the PCA analysis
        %
        % imgListX matrix contains the temporal response for a pixel in a
        % column. The rows represent time samples by number of trials. These are
        % the temporal responses across all trials and time points.
        % 4D array input
        imgListCW = trial2Matrix(absorptions.cw);
        imgListCCW  = trial2Matrix(absorptions.ccw);
        
        
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
        
        % The aligned and offset trials are still matched
        train_index = repmat(train_index, 2, 1);
        
        % Fit the SVM model.
        mdl = fitcsvm(data(train_index, :), label(train_index), ...
            'KernelFunction', 'linear');
        
        % predict the data not in the training set.
        yp = predict(mdl, data(~train_index, :));
        classLoss = sum(label(~train_index) ~= yp) / length(yp);
        
        % X(bb) = barOffset(bb);
        P(pa==polarAngles,c==contrastLevels) = (1-classLoss) * 100;
        
        
    end
    
end


disp(P);

% Visualize
labels = {'Polar Angle: 0','Polar Angle: 90','Polar Angle: 180','Polar Angle: 270'};

colors = lines(4);
figure(1); clf; set(gcf,'Color','w'); hold all;
plot(contrastLevels,P(1,:),'Color', colors(1,:), 'LineWidth',2);
plot(contrastLevels,P(2,:),'Color', colors(2,:), 'LineWidth',2);
plot(contrastLevels,P(3,:),'Color', colors(3,:), 'LineWidth',2);
plot(contrastLevels,P(4,:),'Color', colors(4,:), 'LineWidth',2);
legend(labels);
box off;
xlabel('Contrast level (Michelson)');
ylabel('Classifier Accuracy')
set(gca, 'XLim', [0.4 1], 'YLim', [0 100],'TickDir','out','TickLength',[.015 .015]);

return




