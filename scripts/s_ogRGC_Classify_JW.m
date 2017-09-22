%% s_ogRGC_Classify

% Script with first attempt to classify oriented gabors simulated at
% 7 different contrast levels, 4 polar angles.

FFTflag = false;
%% Classify

contrastLevels = [0:0.01:0.09, 0.1:0.1:1.0]; %flip((0:.1:1).^2);%
polarAngles    = 0; % [0 90 180 270];
eyemovement    = {'000'};%{'000', '100', '010', '001'};
noise          = 'random';
eccen          = 6;
defocus        = 0;
P = nan(length(polarAngles),length(contrastLevels),length(eyemovement));

for pa = polarAngles
    for c = contrastLevels
        for em = 1:length(eyemovement)
            % Load dataset
            fname = sprintf(...
                'OGconeOutputs_contrast%1.2f_pa%d_eye%s_eccen%1.2f_defocus%1.2f_noise-%s.mat',...
                    c,pa,eyemovement{em},  eccen, defocus, noise);                                
            pth = fullfile(ogRootPath, 'data', fname);
            if ~exist(pth, 'file'), error('The file %s is not found', fname); end   
            load(pth);
            
            fprintf('Loading and classifying %s\n', fname);
            % Get the trials and samples (should be the data for all data sets though
            nTrials = size(absorptions,1) * size(absorptions,5)/2;
            tSamples = size(absorptions,4);
            sz = size(absorptions,2);
            
            if ~isstruct(absorptions)
                
                absorptionsPh = absorptions;
                absorptions = [];
                
                % Permute so that nTrials and diff phase conditions are last
                absorptionsPh = permute(absorptionsPh,[2 3 4 1 5]);

                % Reshape to get all trials together 
                absorptions.ccw = reshape(absorptionsPh(:,:,:,:,1:2),[sz, sz, tSamples, nTrials]);
                absorptions.cw = reshape(absorptionsPh(:,:,:,:,3:4),[sz, sz, tSamples, nTrials]);

                absorptions.ccw = absorptions.ccw(:,:,:,randperm(nTrials,nTrials));
                absorptions.cw = absorptions.cw(:,:,:,randperm(nTrials,nTrials)); 
                
                absorptions.ccw = permute(absorptions.ccw, [4 1 2 3]);
                absorptions.cw = permute(absorptions.cw, [4 1 2 3]);
            end
            
            % If requested, fourier transform the cone array outputs
            if FFTflag
                absorptions.cwF = abs(fft2(permute(absorptions.cw, [2 3 1 4])));
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
            
            %             % compute the imagebases of the two stimuli
            %             imageBasis = ogPCA(cat(1,absorptions.cw,absorptions.ccw));
            %
            %             % Time series of weights
            %             weightSeries  = imgList * imageBasis;
            weightSeries = imgList;
            
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
                                  
            % Fit the SVM model.
            mdl = fitcsvm(data, label, 'KernelFunction', 'linear');
            
            cvmdl = crossval(mdl);
            
            % predict the data not in the training set.
            classLoss = kfoldLoss(cvmdl);
            
            P(pa==polarAngles,c==contrastLevels,em) = (1-classLoss) * 100
            
            
        end
    end
end


 disp(P);

% Visualize
labels = {'Polar Angle: 0'};%,'Polar Angle: 90','Polar Angle: 180','Polar Angle: 270'};

colors = lines(length(eyemovement));
figure(1); clf; set(gcf,'Color','w'); hold all;
plot(contrastLevels, squeeze(P),'o-', 'Color', 'k', 'LineWidth',2);
set(gca, 'XScale','log', 'YLim', [0 100], 'TickDir','out','TickLength',[.015 .015]);
ylabel('Classifier Accuracy')
xlabel('Contrast level (Michelson)')

fname = sprintf(...
                'Classify_coneOutputs_contrast%1.2f_pa%d_eye%s_eccen%1.2f_defocus%1.2f_noise-%s_phasescrambled',...
                    c,pa,eyemovement{em},  eccen, defocus, noise);
savefig(fullfile(ogRootPath, 'figs', sprintf('%s.fig', fname)))
hgexport(gcf,fullfile(ogRootPath, 'figs', sprintf('%s.eps', fname)))


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




