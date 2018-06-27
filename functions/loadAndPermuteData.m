function [data, fname] = loadAndPermuteData(expParams, c, em, eccen, df, sf, currentFlag, subFolderName_toLoad)


% Load dataset
fname = sprintf(...
    'OGconeOutputs_contrast%1.3f_pa%d_eye%s_eccen%1.2f_defocus%1.2f_noise-random_sf%1.2f.mat',...
    expParams.contrastLevels(c),expParams.polarAngle,sprintf('%i',expParams.eyemovement(:,em)), expParams.eccentricities(eccen), expParams.defocusLevels(df), sf);

if currentFlag
    fname = ['current_' fname];
end

pth = fullfile(ogRootPath, 'data', expParams.name, subFolderName_toLoad, fname);
if ~exist(pth, 'file'), error('The file %s is not found', fname); end

tmp = load(pth);

if currentFlag
    data = getfield(tmp,'current');
else
    data = getfield(tmp,'absorptions');
    data = data(:,:,:,1:28,:); % Truncate time samples, where blank stimulus was presented.
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