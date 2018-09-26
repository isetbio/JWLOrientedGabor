%% s_makeFigure8_EffectOfConeDensity

% Script to plot figure 8 (effect of cone density on computational observer
% performance

% Requires computational model to be ran with 'eccbasedcoverage' condition, 
% absorptions to be classified and averaged across iterations.


%% Visualize different cone mosaics (Figure 8A)



%% Plot effect of cone density on psychometric function (Figure 8B)
plotPsychometricFunctions('eccbasedcoverage', false)