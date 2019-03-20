%% s_makeFigure8_EffectOfOptics

% Script to plot figure 8 (effect of defocus on computational observer
% performance


% Requires computational model to be ran with 'defocus' condition, 
% absorptions to be classified and averaged across iterations.


%% Visualize different levels of defocus (Figure 8A)
s_plotMTF_defocusLevels


%% Plot effect of eye movements on psychometric function (Figure 8B)
plotPsychometricFunctions('defocus', false)