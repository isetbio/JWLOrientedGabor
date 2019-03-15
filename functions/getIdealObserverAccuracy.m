function P = getIdealObserverAccuracy(data, fname)
% Function to train and test linear SVM classifier on cone data with cross-validation:
%       P = getClassifierAccuracy(data)
%
% INPUTS: 
%   data        : 5 dimensional array (with trials x rows x cols x time
%                                       samples x stimuli) 
% OUTPUTS:
%   P           : classifier accuracy of computational observer model in
%                   percent correct for given absorption dataset

%% 1. Transform and reshape absorption data:

% Get dimensions of data
fprintf('(%s): Loading and classifying\n', mfilename);

% Get the trials and samples (should be the data for all data sets though
nStimuli = size(data,5);
nTrials  = size(data,1) * nStimuli/2;
tSamples = size(data,4);
nrows    = size(data,2);
ncols    = size(data,3);

%   permute to trials x stimuli x rows x cols x time points
data = permute(data, [1 5 2:4]);

%   reshape to (trials x stimuli) x rows x cols x time points
data = reshape(data, [], nrows, ncols, tSamples);

% permute to rows x cols x (trials x stimuli) x time points
data  = permute(data, [2 3 1 4]);

% Compute fourier transform the cone array outputs
% data  = abs(fft2(data));

% reshape to all trials x [rows x colums x time] for classification
data = permute(data, [3 1 2 4]);
data = reshape(data, nTrials*2, []);

% permute the trial order within each of the two classes
idx = [randperm(nTrials) randperm(nTrials)+nTrials];

data = data(idx, :);

%% 2. Transform and reshape template

template = load(fullfile(ogRootPath, 'data', 'idealobserver', 'idealtemplate', fname));
template = template.absorptions;

%   permute to trials x stimuli x rows x cols x time points
template = permute(template, [1 5 2:4]);

%   reshape to (trials x stimuli) x rows x cols x time points
template = reshape(template, [], nrows, ncols, tSamples);

% permute to rows x cols x (trials x stimuli) x time points
template  = permute(template, [2 3 1 4]);

% Compute fourier transform the cone array outputs
% template  = abs(fft2(template));

% reshape to all trials x [rows x colums x time] for classification
template = permute(template, [3 1 2 4]);
% template = reshape(template, nTrials*2, []);

% Label clockwise and counterclockwise trials
label = [ones(nTrials, 1); -ones(nTrials, 1)]; % First set is CW, second set is CCW

%% 3. Define CW and CCW templates, compute log difference

% Get all trials with the same label and only take one trial (since data
% are all the same across time, without any photon noise or eyemovement or
% phase shifts)
templateCW = template(label==1,:,:,1);
templateCW = templateCW(1,:);

templateCCW = template(label==-1,:,:,1);
templateCCW = templateCCW(1,:);

% Check where templates are different, only those datapoints will be used
% for the calculation
diffIdx = find(templateCW~=templateCCW);

for ii = 1:nTeObservations       
    llCW(ii) = sum(log10(poisspdf(data(ii,diffIdx),templateCW(diffIdx))));
    llCCW(ii) = sum(log10(poisspdf(data(ii,diffIdx),templateCCW(diffIdx))));
    logLikelihoodRatio(ii) = llCW(ii)-llCCW(ii);
    
    if logLikelihoodRatio > 0
        predictedLabel(ii) = 1;
    else
        predictedLabel(ii) = -1;
    end
end



% % Compute a log difference between CW and CCW
% logTemplateDiff = log(templateCW(diffIdx))-log(templateCCW(diffIdx));
% C = sum(templateCW(diffIdx)-templateCCW(diffIdx));
% 
% %% 4. Get loglikelihood ratio of data using the template
% for ii = 1:nTrials*2
%     
%     logLikelihoodRatio(ii) = sum( data(ii,diffIdx) .* logTemplateDiff ) + C;
%     
%     
% end


% elseif strcmp(solutionType, 'svm')
% %% 5.Predict stimulus category with weights from ideal observer template
% 
% % Fit the SVM model.
% cvmdl = fitcsvm(template, label, 'Standardize', true, 'KernelFunction', 'linear');
% 
% % Predict new trials
% predictedLabel = predict(cvmdl, data);
% 
% % Transpose
% predictedLabel = predictedLabel';
% 
% end

% Get percent accuracy
P = 100*(sum([predictedLabel(1:nTrials)>0, predictedLabel((nTrials+1):end)<0])/(nTrials*2));

return