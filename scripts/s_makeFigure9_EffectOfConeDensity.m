%% s_makeFigure9_EffectOfConeDensity

% Script to plot figure 9 (effect of cone density on computational observer
% performance

% Requires computational model to be ran with 'conedensity' condition, 
% absorptions to be classified and averaged across iterations.


%% Visualize different cone mosaics (Figure 9A)
plotExampleAbsorptionsFromConeMosaic('conedensity')


%% Plot effect of cone density on psychometric function (Figure 9B)
plotPsychometricFunctions('conedensity', false)