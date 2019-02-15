function [emPaths, cMosaic] = getEyemovements(OG, cMosaic, expParams, sparams, emIdx)

% Calculate number of eyemovements based on cone mosaic integration time
maxEyeMovementsNum = OG(1).maxEyeMovementsNumGivenIntegrationTime(cMosaic.integrationTime);

% Check what eyemovements to simulate:
if all(expParams.eyemovement(:,emIdx) == [1;0])      % if only drift, no MS
    emPaths = cMosaic.emGenSequence(maxEyeMovementsNum*2, 'nTrials', expParams.nTrials, 'microsaccadeType', 'none', 'rSeed', seed);
elseif all(expParams.eyemovement(:,emIdx) == [1;1])  % if drift and MS
    emPaths = cMosaic.emGenSequence(maxEyeMovementsNum*2, 'nTrials', expParams.nTrials, 'microsaccadeType', 'stats based');
elseif all(expParams.eyemovement(:,emIdx) == [0;0]) % if none
    emPaths = zeros(expParams.nTrials, maxEyeMovementsNum*2, 2);
end

% Truncate warm up period
emPaths = emPaths(:, end-maxEyeMovementsNum+1:end,:);

if expParams.verbose
    %plot eye movements
    figure,
    subplot(211)
    plot(sparams.tsamples, emPaths(:,:,1)')
    
    subplot(212)
    plot(sparams.tsamples, emPaths(:,:,2)')
end

return

