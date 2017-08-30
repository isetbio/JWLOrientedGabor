function s = ogFitWeibull(vars, levels, nCorrect, nTotal)
% There are two parameters to be fitted slope and threshold

% INPUTS: 
% vars
%     vars.b = slope
%     vars.t = threshold yeilding ~80% correct

% levels     = independent variable tested on
% nCorrect   = number of correct trials per level
% nTotal     = amount of trials per level

% OUTPUTS:
% s          = -loglikelihood between data and weibull fit


figureOn = 0;

%% compute Weibull

target = ogWeibull(vars, levels);
target = target.*.99 + 0.005;

%% compute log likelihood between data and Weibull fit

% size(target)
% size(nCorrect)

s = -ogLogLikelihood(target, nCorrect', nTotal*ones(1, length(nCorrect)), 'log');

%% visualize

if figureOn
    figure (1), clf
    plot(levels, target), hold on
    plot(levels, nCorrect./nTotal), 
    drawnow
end