%% s_ogRGC_Classify

% Script with first attempt to classify oriented gabors simulated at
% 7 different contrast levels, 4 polar angles.

%% Classify

contrastLevels = 0:0.01:0.1; 
polarAngles    = 0; % [0 90 180 270];
eyemovement    = {'110'};%{'000', '100', '010', '001'};
noise          = 'random';
eccen          = 0; %2 5 10 20 40];
defocus        = 0;
spatFreq       = [0.25, 0.4, 0.65, 1, 1.6, 2.6, 4, 8, 10, 16, 26];

% Predefine matrix for predictions
P = nan(length(polarAngles),length(contrastLevels),length(eyemovement),length(spatFreq));

for pa = polarAngles
    for c = contrastLevels
        for em = 1:length(eyemovement)
            for sf = spatFreq
 
            % Load dataset
            fname = sprintf(...
                'OGconeOutputs_contrast%1.2f_pa%d_eye%s_eccen%1.2f_defocus%1.2f_noise-%s_sf%1.2f.mat',...
                    c,pa,eyemovement{em},  eccen, defocus, noise, sf);                                
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
                        
            % Compute fourier transform the cone array outputs
            absorptions  = abs(fft2(absorptions)); 

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
            
            P(pa==polarAngles,c==contrastLevels,em,sf==spatFreq) = (1-classLoss) * 100;
            
            end
        end
    end
end

disp(P);

% Visualize
labels = {'Polar Angle: 0'};%,'Polar Angle: 90','Polar Angle: 180','Polar Angle: 270'};

colors = lines(length(eyemovement));
figure; clf; set(gcf,'Color','w'); hold all;
plot(contrastLevels, squeeze(P),'o-', 'LineWidth',2);
set(gca, 'XScale','log', 'XLim', [.008 .06], 'XTick', (1:6)/100, ...
    'YLim', [40 100], 'TickDir','out','TickLength',[.015 .015]);
ylabel('Classifier Accuracy')
xlabel('Contrast level (Michelson)')

fname = sprintf(...
                'Classify_coneOutputs_contrast%1.2f_pa%d_eye%s_eccen%1.2f_defocus%1.2f_noise-%s_phasescrambled',...
                    c,pa,eyemovement{em},  eccen, defocus, noise);
save(fullfile(ogRootPath, 'data', sprintf('%s.mat', fname)),'P')
                
savefig(fullfile(ogRootPath, 'data', sprintf('%s.fig', fname)))
hgexport(gcf,fullfile(ogRootPath, 'data', sprintf('%s.eps', fname)))


%% visualize multiple classifier accuracy's
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




