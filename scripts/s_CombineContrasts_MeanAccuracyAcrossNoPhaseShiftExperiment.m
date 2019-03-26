%% s_CombineContrasts_MeanAccuracyAcrossNoPhaseShiftExperiment.m
%
% Combine low and high contrast and compute mean over no eyemovements, no
% stim phase shift

%% 0. Set general experiment parameters
expName                  = 'eyemov';
expParams                = loadExpParams(expName, false);
[xUnits, colors, labels, M] = loadWeibullPlottingParams(expName);

% Use cone current (flag = true) or cone absorptions (flag = false)
polarAngles = expParams.polarAngle;
FFTflag     = true;

% Where to find data and save figures
subFolderName = 'average_noStimPhaseShift';
dataPth     = fullfile(ogRootPath,'data','PF_data_alias','classification',expName);
savePth     = fullfile(dataPth, subFolderName);
figurePth   = fullfile(ogRootPath,'figs', expName, subFolderName);

% Number of total trials in computational observer model (50 clockwise, 50 counterclockwise)
nTotal      = 100;


% Get nr of conditions
nrEyemovTypes    = size(expParams.eyemovement,2);
lmsRatio         = expParams.cparams.spatialDensity;

cmap = copper(nrEyemovTypes);

figure(3); clf;  hold all;
xlabel('Contrast (%)'); ylabel('Accuracy (% correct)'); title('Classifier performance')
set(gca, 'TickDir', 'out', 'FontSize', 15, 'XScale','log', 'LineWidth',2); box off;
for em = 1:nrEyemovTypes
    
    fNamePre1   = sprintf('Classify_coneOutputs_contrast0.100_pa%d_eye%s_eccen%1.2f_defocus%1.2f_noise-random_sf%1.2f_lms-%1.1f%1.1f%1.1f.mat', ...
        polarAngles,sprintf('%i',expParams.eyemovement(:,em)'),expParams.eccentricities,expParams.defocusLevels,expParams.spatFreq, lmsRatio(2),lmsRatio(3),lmsRatio(4));
    
    fNamePre2   = sprintf('Classify_coneOutputs_contrast0.015_pa%d_eye%s_eccen%1.2f_defocus%1.2f_noise-random_sf%1.2f_lms-%1.1f%1.1f%1.1f.mat', ...
        polarAngles,sprintf('%i',expParams.eyemovement(:,em)'),expParams.eccentricities,expParams.defocusLevels,expParams.spatFreq, lmsRatio(2),lmsRatio(3),lmsRatio(4));
    
    
    d = dir(fullfile(dataPth, 'run*'));
    
    P =[];
    
    for ii = 1:2:size(d,1)
        
        fprintf('Load file name: %s %s\n', d(ii).name, fNamePre1);
        
        accuracy1 = load(fullfile(d(ii).folder, d(ii).name, fNamePre1));
        if em==1
            
            fprintf('Load file name: %s %s\n', d(ii+1).name, fNamePre2);

            accuracy2 = load(fullfile(d(ii+1).folder, d(ii+1).name, fNamePre2));
            accuracy.P = [squeeze(accuracy2.accuracy), squeeze(accuracy1.accuracy(5:end))];
        else
            accuracy.P = [squeeze(accuracy1.accuracy)];
        end
        
        if size(accuracy.P,1)<size(accuracy.P,2)
            accuracy.P = accuracy.P';
        end
        
        P = [P accuracy.P];
        %         fprintf('Accuracy is: %d \t',P)
    end
    
    if em ==1
        % Redefine contrast levels
        expParams.contrastLevels = [0:0.001:0.01, 0.015:0.005:0.04, 0.05:0.01:0.1];
    else
        expParams.contrastLevels = [0:0.005:0.04, 0.05:0.01:0.1];
    end
    
    for plotIdx = 1:size(P,2)
        plot(expParams.contrastLevels(2:end),P(2:end,plotIdx), 'LineWidth',1,'Color', cmap(em,:)');
        plot(10.^(-4),P(1,plotIdx), 'o', 'LineWidth',1,'Color', cmap(em,:)');
    end
    
    P_SE = std(P,[],2)./sqrt(size(P,2));
    
    P = mean(P,2);
    fNamePostMean   = sprintf('Classify_coneOutputs_contrast%1.3f_pa%d_eye%s_eccen%1.2f_defocus%1.2f_noise-random_sf%1.2f_lms-%1.1f%1.1f%1.1f_AVERAGE.mat', ...
        max(expParams.contrastLevels),polarAngles,sprintf('%i',expParams.eyemovement(:,em)'),expParams.eccentricities,expParams.defocusLevels,expParams.spatFreq,lmsRatio(2),lmsRatio(3),lmsRatio(4));
    
    fNamePostSE   = sprintf('Classify_coneOutputs_contrast%1.3f_pa%d_eye%s_eccen%1.2f_defocus%1.2f_noise-random_sf%1.2f_lms-%1.1f%1.1f%1.1f_SE.mat', ...
        max(expParams.contrastLevels),polarAngles,sprintf('%i',expParams.eyemovement(:,em)'),expParams.eccentricities,expParams.defocusLevels,expParams.spatFreq,lmsRatio(2),lmsRatio(3),lmsRatio(4));
    
    % plot mean in thick line
    plot(expParams.contrastLevels(2:end), P(2:end), 'LineWidth', 4, 'Color', cmap(em,:)');
    plot(10.^(-4), P(1), 'o', 'LineWidth', 4, 'Color', cmap(em,:)');
    
    if ~exist(savePth,'dir'), mkdir(savePth); end
    save(fullfile(dataPth, subFolderName, fNamePostMean),'P');
    save(fullfile(dataPth, subFolderName, fNamePostSE),'P_SE');
    
end

h = findobj(gca,'Type','line', '-and',{'LineWidth',4});
legend([h(6), h(4), h(2)],labels, 'Location', 'Best')
%
% legendMatrix = strings(8,3);
% legendMatrix(1,1) = labels{1};
% legendMatrix(1,2) = labels{2};
% legendMatrix(1,3) = labels{3};
%
% legend(cellstr(legendMatrix(:)))