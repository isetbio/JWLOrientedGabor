%% s_ogConeAbsorptions_Classify
%
% Script to classify cone absorption rates for clockwise and counter-
% clockwise oriented gabors simulated at different contrast levels.
%
% Absorption data are simulated with script: s_ogConeAbsorptions.m
% using the ISETBIO toolbox.
%
% By Eline Kupers & Jonathan Winawer, NYU (2018)
%
%% Set params

% Reset random number generator seed
rng;

selectTimePoints = 1:109; % for absorptions use 1:28, for cone current use all;

% Load experiment parameters
expName = 'conedensity';
for runNr = 1:5
    subFolderName_toLoad = sprintf('run%d', runNr);
    subFolderName_toSave = sprintf('run%d', runNr);
    expParams = loadExpParams(expName);
    noiseFlag = 'random'; % photon noise properties, saved in file name, could also be 'none' or 'random'
    
    % Compute accuracy for cone current as well
    currentFlag    = true;
    
    % Compute accuracy on fft component of cone absorptions
    fftFlag        = true;
    
    % Predefine matrix for predictions
    if currentFlag
        contrasts    = expParams.contrastLevelsPC;
    else
        contrasts    = expParams.contrastLevels;
    end
    nrContrasts      = length(contrasts);
    nrEyemovTypes    = size(expParams.eyemovement,2);
    nrEccen          = length(expParams.eccentricities);
    nrSpatFreq       = length(expParams.spatFreq);
    nrDefocusLevels  = length(expParams.defocusLevels);
    
    % Preallocate space
    %     P = nan(nrContrasts,1);
    
    % Folder to save data
    % serverPth = '/Volumes/server/Projects/PerformanceFieldsIsetBio';
%     serverPth = '/Volumes/server-1/Projects/PerformanceFields_RetinaV1Model';
    serverPth = '/scratch/ek99/JWLOrientedGabor/';
    savePth = fullfile(serverPth, 'data', expName, 'classification', 'current', subFolderName_toSave);
    if ~exist('savePth', 'dir'); mkdir(savePth); end
    
    % Init figure
%     figure; clf; set(gcf,'Color','w'); hold all;
%     set(gca, 'XScale','log', 'XLim', [.005 max(contrasts)], 'XTick', [1:7, 10:10:100]/100, ...
%         'YLim', [40 100], 'TickDir','out','TickLength',[.015 .015]);
%     ylabel('Classifier Accuracy')
%     xlabel('Contrast level (Michelson)')
    em = 1;
    df = 1;
    sf = expParams.spatFreq;
    for eccen = 1:nrEccen
        %         for df = 1:nrDefocusLevels
        %             for em = 1:max(nrEyemovTypes)
        %                 for sf = expParams.spatFreq
        for c = 1:nrContrasts
            
            % Load dataset
            fname = sprintf(...
                'OGconeOutputs_contrast%1.4f_pa%d_eye%s_eccen%1.2f_defocus%1.2f_noise-%s_sf%1.2f_lms-0.60.30.1.mat',...
                contrasts(c),expParams.polarAngle,sprintf('%i',expParams.eyemovement(:,em)), expParams.eccentricities(eccen), expParams.defocusLevels(df), noiseFlag, sf);
            
            if currentFlag
                fname = ['current_' fname];
            end
            pth = fullfile(serverPth, 'data', expName,'current', expName, subFolderName_toLoad, fname);
            if ~exist(pth, 'file'), error('The file %s is not found', fname); end
            tmp = load(pth);
            
            if currentFlag
                data = getfield(tmp,'current');
                % photocurrent responses are temporally delayed,
                % select same nr of time samples as for absorptions
                % that include stimulus "on" period.
                data = data(:,:,:,selectTimePoints,:);
            else
                data = getfield(tmp,'absorptions');
                % truncate time samples (only include stimulus "on" period)
                data = data(:,:,:,selectTimePoints,:);
            end
            
            
            fprintf('Loading and classifying %s\n', fname);
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
            if fftFlag; data  = abs(fft2(data)); end
            
            % reshape to all trials x [rows x colums x time] for classification
            data = permute(data, [3 1 2 4]);
            data = reshape(data, nTrials*2, []);
            
            % permute the trial order within each of the two classes
            idx = [randperm(nTrials) randperm(nTrials)+nTrials];
            
            data = data(idx, :);
            
            label = [ones(nTrials, 1); -ones(nTrials, 1)];
            
            % Fit the SVM model.
            cvmdl = fitcsvm(data, label, 'Standardize', true, 'KernelFunction', 'linear', 'kFold', 10);
            
            % predict the data not in the training set.
            classLoss = kfoldLoss(cvmdl);
            
            P = (1-classLoss) * 100;
            disp(P);
            
            % Save classifier accuracy
            fname = sprintf(...
                'Classify_coneOutputs_contrast%1.4f_pa%d_eye%s_eccen%1.2f_defocus%1.2f_noise-random_sf%1.2f',...
                contrasts(c), expParams.polarAngle,sprintf('%i',expParams.eyemovement(:,em)), expParams.eccentricities(eccen), expParams.defocusLevels(df), sf);
            if currentFlag; fname = ['current_' fname]; end
            parsave(fullfile(savePth, sprintf('%s.mat', fname)),'P',P)
            
        end
        % Visualize
        %                     plot(contrasts, P,'o-', 'LineWidth',2); drawnow;
    end
    %             end
    %         end
end


return




