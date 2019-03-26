%% s_makeFigure7_EffectOfConeType
%
% Script to plot figure 7 (effect of cone type on computational observer
% performance
%
% Requires computational model to be ran with 'conetypes' condition, 
% absorptions to be classified and averaged across iterations.


%% Visualize different cone mosaics with single cone type(Figure 8A - Right panel)
plotExampleAbsorptionsFromConeMosaic('conetypes')

%% Plot effect of single cone type on model performance as a psychometric function (Figure 8A - Left panel)
plotPsychometricFunctions('conetypes', 'average', false)

%% Plot effect of L-M mixture mosiacs on model performance as a psychometric function (Figure 8B and 8C)
plotPsychometricFunctions('conetypesmixed', 'average', false)