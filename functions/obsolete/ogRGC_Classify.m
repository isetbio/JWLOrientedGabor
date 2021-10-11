function []=ogRGC_Classify(contrastLevels,polarAngles,eyemovement,FFTflag, phaseFlag, ...
                            usedEccentricities, zernikeDefocus, pcaFlag)

% Function to classify absorption rate outputs from s_ogRGC.m
%
% We transform absorptions rates (originally in a ntrials x rows x cols x time)
% for both clockwise (cw) and counter-clockwise (ccw) stimuli to a (space x
% time) matrix. Then absorption rates are concatenated

% We then take either the amplitudes of the


% Examples:
% Defocus
% ogRGC_Classify([],0,{'110'},0, 0, 6, [0 0.5 1 1.5 2])

% Default, no FFT
% ogRGC_Classify([], 0,{'110'},0, 0, 6, [])

% Default, with FFT
% ogRGC_Classify([], 0,{'110'},1, 0, 6, [])

% Eyemovements, no FFT
% ogRGC_Classify([], 0,{'000','110','220','330'},0, 0, 6, [],0)

% Eyemovements, with FFT
% ogRGC_Classify([], 0,{'000','110','220','330'},1, 0, 6, [],0)

% Density
% ogRGC_Classify(0.05, 0,{'110'},0, 0, [1:1:40], [],0)
% ogRGC_Classify(0.08, 0,{'110'},0, 0, [1:1:40], [],0)


% Check inputs
if ~exist('contrastLevels','var') || isempty(contrastLevels)
    contrastLevels = [0 0.01:0.01:0.09, 0.1:0.1:1.0];
end

if ~exist('polarAngles','var') || isempty(polarAngles)
    polarAngles = 0;
end

if ~exist('eyemovement','var') || isempty(eyemovement)
    eyemovement = {'000','100','010','110'};
end

if ~exist('FFTflag','var') || isempty(FFTflag)
    FFTflag = false;
end

if ~exist('phaseFlag','var') || isempty(phaseFlag)
    phaseFlag = true;
    postFix = '';
end

if ~exist('postFix','var') || isempty(postFix)
    postFix = '';
end

if ~exist('usedEccentricities','var') || isempty(usedEccentricities)
    usedEccentricities = 6; % But could be any integer
end

if ~exist('zernikeDefocus','var') || isempty(zernikeDefocus)
    zernikeDefocus = 0; % But could be [0:0.5:2];
end

if ~exist('pcaFlag','var') || isempty(pcaFlag)
    pcaFlag = false; 
end

%% Classify

% contrastLevels = [0.01:0.01:0.09, 0.1:0.1:1.0];
% polarAngles    = 0; % [0 90 180 270];
% eyemovement    = {'110'};%{'000', '100', '010', '001'};
% usedEccentricities = 6; % 2:40;

P = nan(length(polarAngles),length(contrastLevels),length(eyemovement),length(usedEccentricities),length(zernikeDefocus));

for eccen = usedEccentricities
    for pa = polarAngles
        for c = contrastLevels
            for em = 1:length(eyemovement)
                for df = zernikeDefocus
                    % Load dataset
                    load(fullfile(ogRootPath, 'data', sprintf('OGconeOutputs_contrast%1.2f_pa%d_eye%s_eccen%1.2f_defocus%1.2f.mat',c,pa,cell2mat(eyemovement(em)),eccen,df)));
                    
                    % Get the trials and samples (should be the data for all data sets though
                    nTrials = size(absorptions.cw,1);
                    tSamples = size(absorptions.cw,4);
                    
                    % If requested, fourier transform the cone array outputs
                    if FFTflag
                        % permute to put 2D spatial array first
                        absorptions.cw  = permute(absorptions.cw,  [2 3 1 4]);
                        absorptions.ccw = permute(absorptions.ccw, [2 3 1 4]);
                        
                        absorptions.cw  = fft2(absorptions.cw);
                        absorptions.ccw = fft2(absorptions.ccw);
                        
                        if phaseFlag
                            %  compute phase spectrum
                            absorptions.cw  = angle(absorptions.cw);
                            absorptions.ccw = angle(absorptions.ccw);
                            postFix = '_phase';
                        else
                            absorptions.cw  = abs(absorptions.cw);
                            absorptions.ccw = abs(absorptions.ccw);
                        end
                        
                        
                        % unpermute
                        absorptions.cw  = permute(absorptions.cw,  [3 1 2 4]);
                        absorptions.ccw = permute(absorptions.ccw, [3 1 2 4]);
                        
                    end
                    
                    if pcaFlag
                        % Reformat the time series for the PCA analysis
                        %
                        % imgListX matrix contains the temporal response for a pixel in a
                        % column. The rows represent time samples times number of trials.
                        % These are the temporal responses across all trials and time
                        % points. The columns represent the cells.
                        % 4D array input
                        imgListCW  = trial2Matrix(absorptions.cw);
                        imgListCCW = trial2Matrix(absorptions.ccw);
                        
                        % Compute the imagebases of the two stimuli
                        %                     imageBasis = ogPCA(cat(1,absorptions.cw,absorptions.ccw));
                        
                        % Concatenate the matrices of the two stimuli
                        imgList = cat(1,imgListCW,imgListCCW);
                        
                        % Time series of weights
                        weightSeries  = imgList;% * imageBasis;
                        
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
                    else
                        
                        data = cat(1,absorptions.cw,absorptions.ccw);
                        data = reshape(data,[size(data,1),size(data,2)*size(data,3)*size(data,4)]);
                    end
                    
                    label = [ones(nTrials, 1); -ones(nTrials, 1)];
                    
                    % Fit the SVM model and cross validate
                    SVMModel = fitcsvm(data,label,'KernelFunction','linear');
                    CVSVMModel = crossval(SVMModel);
                    classLoss = kfoldLoss(CVSVMModel);
                    
                    % Store performance for each polar angle, contrast level and
                    % eyemovement condition
                    P(pa==polarAngles,c==contrastLevels,em, eccen==usedEccentricities,df==zernikeDefocus) = (1-classLoss) * 100;
                    
                    
                end
            end
        end
    end
end

disp(P);
save(fullfile(ogRootPath,'figs',sprintf('contrastVSperformance_eye%s_pa%d_fft%d%s_eccen%1.2f_defocus%1.2f_pca%d.mat',cell2mat(eyemovement),polarAngles,FFTflag,postFix,max(usedEccentricities),max(zernikeDefocus),pcaFlag)),'P')

% Visualize
% labels = {'Polar Angle: 0'};%,'Polar Angle: 90','Polar Angle: 180','Polar Angle: 270'};

% colors = lines(length(eyemovement));
% figure(1); clf; set(gcf,'Color','w'); hold all;
% plot(contrastLevels, squeeze(P),'Color', 'k', 'LineWidth',2);
% set(gca, 'XScale','log', 'YLim', [0 100], 'TickDir','out','TickLength',[.015 .015]);
% ylabel('Classifier Accuracy')
% xlabel('Contrast level (Michelson)')

%savefig(fullfile(ogRootPath,'figs',sprintf('contrastVSperformance_eye%s_pa%d_fft%d%s',cell2mat(eyemovement),polarAngles,FFTflag,postFix)))
%hgexport(gcf,fullfile(ogRootPath,'figs',sprintf('contrastVSperformance_eye%s_pa%d_fft%d%s.eps',cell2mat(eyemovement),polarAngles,FFTflag,postFix)))


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




