function [emPaths, cMosaic] = getEyemovements(OG, cMosaic, expParams, sparams)

% Calculate number of eyemovements based on cone mosaic integration time
maxEyeMovementsNum = OG(1).maxEyeMovementsNumGivenIntegrationTime(cMosaic.integrationTime);

% Check what eyemovements to simulate:
if all(expParams.eyemovement(:,emIdx) == [1;0])      % if only drift, no MS
    emPaths = cMosaic.emGenSequence(maxEyeMovementsNum*2, 'nTrials', expParams.nTrials, 'microsaccadeType', 'none', 'rSeed', expParams.seed);
elseif all(expParams.eyemovement(:,emIdx) == [1;1])  % if drift and MS
    emPaths = cMosaic.emGenSequence(maxEyeMovementsNum*2, 'nTrials', expParams.nTrials, 'microsaccadeType', 'stats based', 'rSeed', expParams.seed);
elseif all(expParams.eyemovement(:,emIdx) == [0;0]) % if none
    emPaths = zeros(expParams.nTrials, maxEyeMovementsNum*2, 2);
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

