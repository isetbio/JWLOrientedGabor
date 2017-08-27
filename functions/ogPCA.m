function [imageBasisAbsorptions] = ogPCA(coneData)
% ogPCA - Make the principal component images for Oriented Gabor absorptions
%
%   [imageBasisAbsorptions] =ogPCA(params)
%
% The input argument is a struct that includes absorption data (could later
% also include current data)
%
% The basis images for the pattern of absorptions are returned.
%
% Algorithm:
%   Combine the absorptions into their large matrices
%   Compute the PCA across all the stimuli
%   Save the parameters and the image bases.
%
%
% EK & JW, ISETBIO NYU Team, Copyright 2017


%% At some point we want to check if the PCA has already been computed
% fname = vaFname(params);
% 
% if exist(fname,'file')
%     disp('Loading image basis from file - parameters match')
%     load(fname,'imageBasisAbsorptions','imageBasisCurrent');
%     imageBasisAbsorptions = imageBasisAbsorptions(:,1:params.nBasis); %#ok<*NODEF>
%     imageBasisCurrent     = imageBasisCurrent(:,1:params.nBasis);
%     return;
% else 
%     disp('Creating and saving image basis file - parameters do not match')
% end


%% Basic parameters

% Figure out a sensible way to set this.  We average some trials to reduce
% noise, I think.
nTrials = size(coneData,1);

% We want plenty of basis functions stored.  We may not use all of them.
nBasis = 100;

%% Convert the shape of the absorptions so we can perform the svd

tAbsorptions = trial2Matrix(coneData);

% We make bases for each type of stimulus and then concatenate them.
[~,~,V] = svd(tAbsorptions,'econ');
imageBasisAbsorptions = V(:,1:nBasis); %#ok<*NASGU>

% %% Convert the shape of the absorptions so we can perform the svd
% 
% tCurrent = trial2Matrix(current,cMosaic);
% 
% % We make bases for each type of stimulus and then concatenate them.
% [~,~,V] = svd(tCurrent,'econ');
% imageBasisCurrent = V(:,1:nBasis);

%% Save full basis set, but return the params.nBasis number

% save(fname,'imageBasisAbsorptions','imageBasisCurrent','params');
% save(fname,'imageBasisAbsorptions','params');

imageBasisAbsorptions = imageBasisAbsorptions(:,1:nBasis); %#ok<*NODEF>
% imageBasisCurrent     = imageBasisCurrent(:,1:params.nBasis);


end
