function result = ogLogLikelihood(pred, nCorrect, nTotal, logOrNot)

% function result = logLikelihood(pred, correctNum, totalNum, logOrNot)


%% Example

exampleOn = 0;

if exampleOn
    pred = 0.1 : 0.05 : 0.99;
    nCorrect = randi(50, [1, length(pred)]);
    nTotal   = 50 .* ones(1, length(pred));
    logOrNot   = 'log';
end

%% validate inputs

checkLth = @(a, b, c) isequal(length(pred), length(nCorrect), length(nTotal));
assert(checkLth(pred, nCorrect, nTotal), 'Length of the inputs do not match.')


%% compute loglikelihood

nWrong = nTotal - nCorrect;

switch lower(logOrNot)
    case 'log'
        result = sum(log(pred).*nCorrect + log(1 - pred).* nWrong);
    case 'nolog'
        result = prod(pred.^(nCorrect).*(1 - pred).^nWrong);
end


end