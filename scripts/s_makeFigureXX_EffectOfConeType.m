%% s_makeFigureXX_EffectOfConeType

% Script to plot figure xx (effect of cone type on computational observer
% performance

% Requires computational model to be ran with 'conetypes' condition, 
% absorptions to be classified and averaged across iterations.


%% Visualize different cone mosaics (Figure 8A)
plotExampleAbsorptionsFromConeMosaic('conetypes')

%% Plot effect of cone density on psychometric function (Figure 8B)
plotPsychometricFunctions('conetypes', 'run1', false)