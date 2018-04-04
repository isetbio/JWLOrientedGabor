% define experiment parameters, save in a mat file so it can be loaded
% later

expName = 'default'; % Choose from 'default','eyemov','eyemovEnhanced','coneDensity','defocus','ConeTypes'

expParams = struct;

switch expName

    case 'default'
        expParams.nTrials         = 25;               % Number of trials per stimulus condition
        expParams.contrast_levels = [0:0.01:0.1 0.2]; % Contrast levels
        expParams.eyemovement     = [1 1 0]';         % Which type of eye movments
        expParams.eccentricities  = 4.5;              % [0 2 5 10 20 40]; %4.5;
        expParams.spatFreq        = 4;                % [0.25, 0.4, 0.65, 1, 1.6, 2.6, 4, 8, 10, 16, 26];
        expParams.polarangles     = 0;
        expParams.defocuslevels   = 0;                % units??  [0 0.5 1 1.5 2]
        expParams.verbose         = true;

     case 'eyemov'
        expParams.nTrials         = 25;               % Number of trials per stimulus condition
        expParams.contrast_levels = [0:0.01:0.1 0.2]; % Contrast levels
        expParams.eyemovement     = {[0 0 0]', [1 0 0]', [0 1 0]', [1 1 0]'};         % Which type of eye movments
        expParams.eccentricities  = 4.5;              % [0 2 5 10 20 40]; %4.5;
        expParams.spatFreq        = 4;                % [0.25, 0.4, 0.65, 1, 1.6, 2.6, 4, 8, 10, 16, 26];
        expParams.polarangles     = 0;
        expParams.defocuslevels   = 0;                % units??  [0 0.5 1 1.5 2]
        expParams.verbose         = true;
     
      case 'coneDensity'
        expParams.nTrials         = 25;               % Number of trials per stimulus condition
        expParams.contrast_levels = [0:0.01:0.1 0.2]; % Contrast levels
        expParams.eyemovement     = [1 1 0]';         % Which type of eye movments
        expParams.eccentricities  = [0 2 4.5 5 10 20 40];
        expParams.spatFreq        = 4;                % [0.25, 0.4, 0.65, 1, 1.6, 2.6, 4, 8, 10, 16, 26];
        expParams.polarangles     = 0;
        expParams.defocuslevels   = 0;                % units??  [0 0.5 1 1.5 2]
        expParams.verbose         = true;
end


save(fullfile(ogRootPath, 'Stimulus', sprintf('ogDefault_%s.mat',expName)),expParams)