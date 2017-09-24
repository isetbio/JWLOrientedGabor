%% s_ogRGC_Classify

% Script with first attempt to classify oriented gabors simulated at
% 7 different contrast levels, 4 polar angles.

FFTflag = true;
%% Classify

contrastLevels = flip(([1:6 10])/100);% flip(0.01:0.01:0.1); % [0:0.01:0.09, 0.1:0.1:1.0]; %flip((0:.1:1).^2);%
polarAngles    = 0; % [0 90 180 270];
eyemovement    = {'110'};%{'000', '100', '010', '001'};
noise          = 'random';
eccen          = 4.5;
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
            nStimuli = size(absorptions,5);
            nTrials  = size(absorptions,1) * nStimuli/2;
            tSamples = size(absorptions,4);
            nrows    = size(absorptions,2);
            ncols    = size(absorptions,3);            
            % absorptions is trials x rows x cols x time points x stimuli            
            %   permute to trials x stimuli x rows x cols x time points
            absorptions = permute(absorptions, [1 5 2:4]);
 
            %   reshape to (trials x stimuli) x rows x cols x time points                         
            absorptions = reshape(absorptions, [], nrows, ncols, tSamples);
            
            % permute to rows x cols x (trials x stimuli) x time points
            absorptions  = permute(absorptions, [2 3 1 4]);
                        
            % If requested, fourier transform the cone array outputs
            if FFTflag, absorptions  = abs(fft2(absorptions)); end

            
            %             s_known = absorptions(39,5,:,:) + ...
            %                 absorptions(2,36,:,:) - ...
            %                 absorptions(2,5,:,:) - ...
            %                 absorptions(39,36,:,:) ;
            %             s_known = squeeze(s_known);
            %
            %             s_known(end/2+1:end,:) = s_known(end/2+1:end,:) * -1;
            %             snr(pa==polarAngles,c==contrastLevels,em) = mean(s_known(:))/std(s_known(:));
            %             disp(snr)
                        
            % reshape to trials x everything else for classification
            absorptions = permute(absorptions, [3 1 2 4]);

            absorptions = reshape(absorptions, nTrials*2, []);

            % permute the trial order within each of the two classes
            idx = [randperm(nTrials) randperm(nTrials)+nTrials];
            
            absorptions = absorptions(idx, :);

            label = [ones(nTrials, 1); -ones(nTrials, 1)];
                                  
            % Fit the SVM model.
            mdl = fitcsvm(absorptions, label, 'KernelFunction', 'linear');
            
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
figure; clf; set(gcf,'Color','w'); hold all;
plot(contrastLevels, squeeze(P),'o-', 'Color', 'k', 'LineWidth',2);
set(gca, 'XScale','log', 'XLim', [.008 .06], 'XTick', (1:6)/100, ...
    'YLim', [40 100], 'TickDir','out','TickLength',[.015 .015]);
ylabel('Classifier Accuracy')
xlabel('Contrast level (Michelson)')

fname = sprintf(...
                'Classify_coneOutputs_contrast%1.2f_pa%d_eye%s_eccen%1.2f_defocus%1.2f_noise-%s_phasescrambled_fft0',...
                    c,pa,eyemovement{em},  eccen, defocus, noise);
save(fullfile(ogRootPath, 'figs', sprintf('%s.mat', fname)),'P')
                
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




