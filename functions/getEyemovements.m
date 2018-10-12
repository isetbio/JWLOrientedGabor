function [emPaths, cMosaic] = getEyemovements(OG, cMosaic, expParams, sparams, emIdx, currSeed)

% Create eye movements for individual trials. Eye movement paths are defined 
% in units of cones along the horizontal (x, first dimension), vertical (y,
% second dimension), for every time step t. This code relies on the
% emGenSequence function of isetbio, which can create drift and
% microsaccades
%
%   [emPaths, cMosaic] = getEyemovements(OG, cMosaic, expParams, sparams, emIdx, currSeed)
%
% INPUTS: 
%   OG              : Struct containing the Optical Image sequences (OIS),
%                       created by getSceneAndStimuli()
%   cMosaic         : Struct containing the isetbio cone mosaic, created by
%                       getConeMosaic()
%   expParams       : Struct containing parameters of this particular
%                       experimental condition, created by loadExpParams()
%   sparams         : Struct containing specific scene parameters, created by 
%                       getSceneAndStimuli()
%   emIdx           : Integer containing index to loop over different eye movement
%                       conditions.The function needs this index to get the 
%                       correct eyemovement condition from the emGenSequence
%                       (TODO: probably can be simplified and defined as a
%                       different input instead)
%   currSeed        : Integer containing the current random number generator
%                       seed used to produce eye movement paths  
%
% INPUTS: 
%   emPaths         : 3 dimensional array contains trials x position (x,y) x time



% Calculate number of eyemovements based on cone mosaic integration time
maxEyeMovementsNum = OG(1).maxEyeMovementsNumGivenIntegrationTime(cMosaic.integrationTime);

% Lengthen the number of eyemovements such that the starting position
% varies across trials and OIS'.
paddingFactor = 8;

% Check what eyemovements to simulate:
if all(expParams.eyemovement(:,emIdx) == [1;0])      % if only drift, no MS
    emPaths = cMosaic.emGenSequence(maxEyeMovementsNum*paddingFactor, 'nTrials', expParams.nTrials, 'microsaccadeType', 'none', 'rSeed', currSeed);
elseif all(expParams.eyemovement(:,emIdx) == [1;1])  % if drift and MS
    emPaths = cMosaic.emGenSequence(maxEyeMovementsNum*paddingFactor, 'nTrials', expParams.nTrials, 'microsaccadeType', 'stats based', 'rSeed', currSeed);
elseif all(expParams.eyemovement(:,emIdx) == [0;0]) % if none
    emPaths = zeros(expParams.nTrials, maxEyeMovementsNum*paddingFactor, 2);
elseif all(expParams.eyemovement(:,emIdx) == [2;0]) % if enhanced drift, no MS
    error('this condition is not implemented yet')
    em = fixationalEM;
    % todo: add line to adjust the velocity/speed/etc
    emPaths = cMosaic.emGenSequence(maxEyeMovementsNum*paddingFactor, 'nTrials', expParams.nTrials, 'microsaccadeType', 'none', 'rSeed', currSeed, 'emobj', em);    
elseif all(expParams.eyemovement(:,emIdx) == [1;1])  % if drift and MS
    error('this condition is not implemented yet')
    em = fixationalEM;
    % todo: add line to adjust the velocity/speed/etc
    emPaths = cMosaic.emGenSequence(maxEyeMovementsNum*paddingFactor, 'nTrials', expParams.nTrials, 'microsaccadeType', 'stats based', 'rSeed', currSeed, 'emobj', em);

end

% Truncate warm up period
emPaths = emPaths(:, end-maxEyeMovementsNum+1:end,:);

if expParams.verbose
    %plot eye movements
    figure(98); clf;
    subplot(211)
    plot(sparams.tsamples, emPaths(:,:,1)')
    xlabel('Time points')
    ylabel('Position (cones)');
    title('Eye movements X position')
    
    subplot(212)
    plot(sparams.tsamples, emPaths(:,:,2)')
    xlabel('Time points')
    ylabel('Position (cones)');
    title('Eye movements Y position')
end

return

