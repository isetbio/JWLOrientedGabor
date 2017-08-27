%% s_ogRGC_Classify

% Script with first attempt to classify oriented gabors simualted at 4
% different polar angles

%% Load the stored data and define parameters
load(fullfile(ogRootPath, 'data', 'OGconeOutputs'));

nTrials = size(absorptions.cw,1);
%% Classify

labels = [ones(nTrials,1);-1*ones(nTrials,1)];


thisData = [absorptions.cw;absorptions.ccw];

m = fitcsvm(thisData, labels, 'KernelFunction', 'linear');
cv = crossval(m,'kfold',5);
rocAreaCones = 1-kfoldLoss(cv);

fprintf('ROC Area for cones: %4.2f\n', rocAreaCones)



